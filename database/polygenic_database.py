"""
Database management for polygenic scores and version control.

Handles:
- Polygenic score storage and retrieval
- Database versioning
- Update checking and execution
"""

import sqlite3
import os
import shutil
import hashlib
import json
from datetime import datetime
from typing import List, Optional, Dict, Tuple
from contextlib import contextmanager

from models.polygenic_models import (
    PolygenicScore, PolygenicVariant, PopulationDistribution,
    DatabaseVersion, UpdateStatus, TraitCategory
)
from config import DATABASE_PATH, BASE_DIR
from utils.logging_config import get_logger

logger = get_logger(__name__)

# Database paths
PGS_DATABASE_PATH = os.path.join(BASE_DIR, 'database', 'pgs.db')
BACKUP_DIR = os.path.join(BASE_DIR, 'database', 'backups')


class PolygenicDatabaseError(Exception):
    """Exception raised for polygenic database errors."""
    pass


class PolygenicDatabase:
    """
    Database manager for polygenic scores.
    
    Provides:
    - Storage and retrieval of polygenic score definitions
    - Population distribution data
    - Version tracking and updates
    """
    
    def __init__(self, db_path: Optional[str] = None) -> None:
        """
        Initialize the database manager.
        
        Args:
            db_path: Path to the database file. Uses default if None.
        """
        self.db_path = db_path or PGS_DATABASE_PATH
        self._ensure_database()
    
    @contextmanager
    def _get_connection(self):
        """Context manager for database connections."""
        conn = None
        try:
            conn = sqlite3.connect(self.db_path)
            conn.row_factory = sqlite3.Row
            yield conn
        except sqlite3.Error as e:
            logger.error(f"Database error: {e}")
            raise PolygenicDatabaseError(f"Database error: {e}")
        finally:
            if conn:
                conn.close()
    
    def _ensure_database(self) -> None:
        """Ensure database exists with required schema."""
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            # Polygenic scores table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS polygenic_scores (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pgs_id TEXT UNIQUE NOT NULL,
                    trait_name TEXT NOT NULL,
                    trait_category TEXT NOT NULL,
                    publication_doi TEXT,
                    publication_year INTEGER,
                    study_population TEXT,
                    sample_size INTEGER,
                    num_variants INTEGER,
                    description TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            # Score variants table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS pgs_variants (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pgs_id TEXT NOT NULL,
                    rsid TEXT NOT NULL,
                    chromosome TEXT,
                    position INTEGER,
                    effect_allele TEXT NOT NULL,
                    other_allele TEXT,
                    effect_weight REAL NOT NULL,
                    effect_allele_frequency REAL,
                    FOREIGN KEY (pgs_id) REFERENCES polygenic_scores(pgs_id),
                    UNIQUE(pgs_id, rsid)
                )
            """)
            
            # Population distributions table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS population_distributions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pgs_id TEXT NOT NULL,
                    population TEXT NOT NULL,
                    mean REAL NOT NULL,
                    std REAL NOT NULL,
                    percentiles_json TEXT,
                    FOREIGN KEY (pgs_id) REFERENCES polygenic_scores(pgs_id),
                    UNIQUE(pgs_id, population)
                )
            """)
            
            # Database versions table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS database_versions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    database_name TEXT UNIQUE NOT NULL,
                    version TEXT NOT NULL,
                    release_date TEXT,
                    download_date TEXT NOT NULL,
                    source_url TEXT,
                    record_count INTEGER,
                    checksum TEXT
                )
            """)
            
            # Create indexes
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_pgs_trait ON polygenic_scores(trait_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_pgs_category ON polygenic_scores(trait_category)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_variant_rsid ON pgs_variants(rsid)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_variant_pgs ON pgs_variants(pgs_id)")
            
            conn.commit()
        
        # Insert sample data if database is empty
        self._ensure_sample_data()
    
    def _ensure_sample_data(self) -> None:
        """Ensure sample polygenic scores exist in database."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM polygenic_scores")
            if cursor.fetchone()[0] > 0:
                return
        
        # Insert sample polygenic scores
        self._insert_sample_scores()
        self._set_version("pgs_catalog", "1.0.0", "https://www.pgscatalog.org/")
    
    def _insert_sample_scores(self) -> None:
        """Insert sample polygenic scores for testing."""
        sample_scores = self._get_sample_pgs_data()
        
        for score_data in sample_scores:
            try:
                self.insert_score(score_data['score'], score_data['distribution'])
            except PolygenicDatabaseError as e:
                logger.warning(f"Error inserting sample score: {e}")
    
    def _get_sample_pgs_data(self) -> List[Dict]:
        """Generate sample polygenic score data."""
        return [
            {
                'score': PolygenicScore(
                    pgs_id="PGS000001",
                    trait_name="Type 2 Diabetes",
                    trait_category=TraitCategory.METABOLIC,
                    publication_doi="10.1038/ng.3951",
                    publication_year=2018,
                    study_population="European",
                    sample_size=898130,
                    num_variants=15,
                    description="Polygenic risk score for type 2 diabetes based on GWAS meta-analysis",
                    variants=[
                        PolygenicVariant("rs7903146", "10", 114758349, "T", "C", 0.295, 0.30),
                        PolygenicVariant("rs1801282", "3", 12393125, "C", "G", 0.143, 0.12),
                        PolygenicVariant("rs5219", "11", 17409572, "T", "C", 0.140, 0.35),
                        PolygenicVariant("rs13266634", "8", 118184783, "T", "C", 0.115, 0.70),
                        PolygenicVariant("rs10811661", "9", 22134095, "T", "C", 0.195, 0.79),
                        PolygenicVariant("rs7756992", "6", 20679709, "G", "A", 0.120, 0.31),
                        PolygenicVariant("rs9939609", "16", 53820527, "A", "T", 0.088, 0.42),
                        PolygenicVariant("rs780094", "2", 27741237, "T", "C", 0.063, 0.39),
                        PolygenicVariant("rs1260326", "2", 27730940, "T", "C", 0.055, 0.40),
                        PolygenicVariant("rs2943641", "2", 227093745, "C", "T", 0.078, 0.63),
                        PolygenicVariant("rs174547", "11", 61570783, "T", "C", 0.045, 0.33),
                        PolygenicVariant("rs560887", "2", 169763148, "C", "T", 0.052, 0.30),
                        PolygenicVariant("rs17782313", "18", 57851097, "C", "T", 0.067, 0.24),
                        PolygenicVariant("rs1558902", "16", 53803574, "A", "T", 0.082, 0.41),
                        PolygenicVariant("rs429358", "19", 45411941, "C", "T", 0.048, 0.14),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000001",
                    population="EUR",
                    mean=1.45,
                    std=0.42,
                    percentiles={5: 0.78, 10: 0.92, 25: 1.17, 50: 1.45, 75: 1.73, 90: 2.01, 95: 2.15}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000002",
                    trait_name="Coronary Artery Disease",
                    trait_category=TraitCategory.CARDIOVASCULAR,
                    publication_doi="10.1038/ng.3913",
                    publication_year=2018,
                    study_population="European",
                    sample_size=547261,
                    num_variants=12,
                    description="Polygenic risk score for coronary artery disease",
                    variants=[
                        PolygenicVariant("rs1333049", "9", 22125500, "C", "G", 0.312, 0.47),
                        PolygenicVariant("rs10757278", "9", 22124478, "G", "A", 0.285, 0.48),
                        PolygenicVariant("rs17465637", "1", 56962821, "C", "A", 0.142, 0.72),
                        PolygenicVariant("rs6511720", "19", 11202306, "T", "G", 0.198, 0.12),
                        PolygenicVariant("rs12740374", "1", 109821511, "T", "G", 0.165, 0.22),
                        PolygenicVariant("rs964184", "11", 116648917, "G", "C", 0.125, 0.18),
                        PolygenicVariant("rs4420638", "19", 45422946, "G", "A", 0.178, 0.17),
                        PolygenicVariant("rs1800961", "20", 43042364, "T", "C", 0.095, 0.05),
                        PolygenicVariant("rs662799", "11", 116663707, "A", "G", 0.108, 0.15),
                        PolygenicVariant("rs1799983", "7", 150999023, "T", "G", 0.088, 0.33),
                        PolygenicVariant("rs5186", "3", 148459988, "C", "A", 0.075, 0.28),
                        PolygenicVariant("rs1801133", "1", 11856378, "T", "C", 0.062, 0.35),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000002",
                    population="EUR",
                    mean=1.72,
                    std=0.38,
                    percentiles={5: 1.12, 10: 1.25, 25: 1.48, 50: 1.72, 75: 1.96, 90: 2.22, 95: 2.38}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000003",
                    trait_name="Breast Cancer",
                    trait_category=TraitCategory.ONCOLOGY,
                    publication_doi="10.1038/s41588-019-0393-x",
                    publication_year=2019,
                    study_population="European",
                    sample_size=228951,
                    num_variants=10,
                    description="Polygenic risk score for breast cancer",
                    variants=[
                        PolygenicVariant("rs2981582", "10", 123337335, "A", "G", 0.268, 0.38),
                        PolygenicVariant("rs889312", "5", 56031884, "C", "A", 0.132, 0.28),
                        PolygenicVariant("rs13281615", "8", 128355618, "G", "A", 0.085, 0.40),
                        PolygenicVariant("rs6504950", "17", 53056471, "G", "A", 0.052, 0.27),
                        PolygenicVariant("rs614367", "11", 69328764, "T", "C", 0.148, 0.16),
                        PolygenicVariant("rs1042522", "17", 7579472, "C", "G", 0.075, 0.72),
                        PolygenicVariant("rs6983267", "8", 128413305, "G", "A", 0.065, 0.48),
                        PolygenicVariant("rs4430796", "17", 36098040, "A", "G", 0.055, 0.48),
                        PolygenicVariant("rs10993994", "10", 51549496, "T", "C", 0.048, 0.41),
                        PolygenicVariant("rs1447295", "8", 128554220, "A", "C", 0.072, 0.11),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000003",
                    population="EUR",
                    mean=0.95,
                    std=0.28,
                    percentiles={5: 0.52, 10: 0.62, 25: 0.78, 50: 0.95, 75: 1.12, 90: 1.32, 95: 1.45}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000004",
                    trait_name="Body Mass Index",
                    trait_category=TraitCategory.METABOLIC,
                    publication_doi="10.1038/nature14177",
                    publication_year=2015,
                    study_population="European",
                    sample_size=339224,
                    num_variants=8,
                    description="Polygenic risk score for BMI/obesity",
                    variants=[
                        PolygenicVariant("rs9939609", "16", 53820527, "A", "T", 0.392, 0.42),
                        PolygenicVariant("rs1558902", "16", 53803574, "A", "T", 0.385, 0.41),
                        PolygenicVariant("rs17782313", "18", 57851097, "C", "T", 0.125, 0.24),
                        PolygenicVariant("rs6265", "11", 27679916, "T", "C", 0.082, 0.20),
                        PolygenicVariant("rs1800497", "11", 113400106, "T", "C", 0.055, 0.19),
                        PolygenicVariant("rs571312", "18", 57839769, "A", "C", 0.095, 0.24),
                        PolygenicVariant("rs2943641", "2", 227093745, "C", "T", 0.068, 0.63),
                        PolygenicVariant("rs10938397", "4", 45186139, "G", "A", 0.078, 0.44),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000004",
                    population="EUR",
                    mean=1.12,
                    std=0.35,
                    percentiles={5: 0.58, 10: 0.70, 25: 0.89, 50: 1.12, 75: 1.35, 90: 1.58, 95: 1.72}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000005",
                    trait_name="Alzheimer's Disease",
                    trait_category=TraitCategory.NEUROPSYCHIATRIC,
                    publication_doi="10.1038/ng.2802",
                    publication_year=2013,
                    study_population="European",
                    sample_size=74046,
                    num_variants=6,
                    description="Polygenic risk score for Alzheimer's disease",
                    variants=[
                        PolygenicVariant("rs429358", "19", 45411941, "C", "T", 1.245, 0.14),
                        PolygenicVariant("rs7412", "19", 45412079, "T", "C", -0.892, 0.08),
                        PolygenicVariant("rs6265", "11", 27679916, "T", "C", 0.085, 0.20),
                        PolygenicVariant("rs1006737", "3", 53127857, "A", "G", 0.062, 0.33),
                        PolygenicVariant("rs4680", "22", 19963748, "G", "A", 0.045, 0.48),
                        PolygenicVariant("rs1800497", "11", 113400106, "T", "C", 0.038, 0.19),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000005",
                    population="EUR",
                    mean=0.42,
                    std=0.58,
                    percentiles={5: -0.55, 10: -0.35, 25: 0.05, 50: 0.42, 75: 0.78, 90: 1.15, 95: 1.42}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000006",
                    trait_name="Prostate Cancer",
                    trait_category=TraitCategory.ONCOLOGY,
                    publication_doi="10.1038/ng.3669",
                    publication_year=2016,
                    study_population="European",
                    sample_size=140306,
                    num_variants=8,
                    description="Polygenic risk score for prostate cancer",
                    variants=[
                        PolygenicVariant("rs4430796", "17", 36098040, "A", "G", 0.215, 0.48),
                        PolygenicVariant("rs1447295", "8", 128554220, "A", "C", 0.385, 0.11),
                        PolygenicVariant("rs10993994", "10", 51549496, "T", "C", 0.248, 0.41),
                        PolygenicVariant("rs6983267", "8", 128413305, "G", "A", 0.125, 0.48),
                        PolygenicVariant("rs1042522", "17", 7579472, "C", "G", 0.068, 0.72),
                        PolygenicVariant("rs2981582", "10", 123337335, "A", "G", 0.045, 0.38),
                        PolygenicVariant("rs889312", "5", 56031884, "C", "A", 0.035, 0.28),
                        PolygenicVariant("rs13281615", "8", 128355618, "G", "A", 0.055, 0.40),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000006",
                    population="EUR",
                    mean=1.08,
                    std=0.32,
                    percentiles={5: 0.58, 10: 0.68, 25: 0.88, 50: 1.08, 75: 1.28, 90: 1.52, 95: 1.65}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000007",
                    trait_name="Rheumatoid Arthritis",
                    trait_category=TraitCategory.IMMUNE,
                    publication_doi="10.1038/ng.2462",
                    publication_year=2012,
                    study_population="European",
                    sample_size=100000,
                    num_variants=7,
                    description="Polygenic risk score for rheumatoid arthritis",
                    variants=[
                        PolygenicVariant("rs2476601", "1", 114377568, "A", "G", 0.542, 0.09),
                        PolygenicVariant("rs3087243", "2", 204738919, "G", "A", 0.125, 0.42),
                        PolygenicVariant("rs6897932", "5", 35910332, "C", "T", 0.088, 0.23),
                        PolygenicVariant("rs11209026", "1", 67705958, "A", "G", 0.175, 0.07),
                        PolygenicVariant("rs2201841", "1", 67649663, "T", "C", 0.095, 0.08),
                        PolygenicVariant("rs10781499", "9", 117552851, "A", "G", 0.072, 0.10),
                        PolygenicVariant("rs3135388", "6", 32681631, "A", "G", 0.285, 0.13),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000007",
                    population="EUR",
                    mean=0.65,
                    std=0.28,
                    percentiles={5: 0.22, 10: 0.32, 25: 0.48, 50: 0.65, 75: 0.82, 90: 1.02, 95: 1.15}
                )
            },
            {
                'score': PolygenicScore(
                    pgs_id="PGS000008",
                    trait_name="Height",
                    trait_category=TraitCategory.PHYSICAL,
                    publication_doi="10.1038/ng.3097",
                    publication_year=2014,
                    study_population="European",
                    sample_size=253288,
                    num_variants=9,
                    description="Polygenic score for adult height",
                    variants=[
                        PolygenicVariant("rs143384", "20", 34025756, "A", "G", 0.082, 0.58),
                        PolygenicVariant("rs1042725", "12", 66339827, "T", "C", 0.055, 0.49),
                        PolygenicVariant("rs6060369", "20", 6604619, "T", "C", 0.042, 0.35),
                        PolygenicVariant("rs2284746", "2", 55308273, "G", "C", 0.065, 0.60),
                        PolygenicVariant("rs6060373", "20", 6608684, "T", "G", 0.038, 0.33),
                        PolygenicVariant("rs12913832", "15", 28365618, "G", "A", 0.025, 0.78),
                        PolygenicVariant("rs1426654", "15", 48426484, "A", "G", 0.018, 0.99),
                        PolygenicVariant("rs16891982", "5", 33951693, "C", "G", 0.022, 0.95),
                        PolygenicVariant("rs4988235", "2", 136608646, "T", "C", 0.015, 0.50),
                    ]
                ),
                'distribution': PopulationDistribution(
                    pgs_id="PGS000008",
                    population="EUR",
                    mean=0.38,
                    std=0.12,
                    percentiles={5: 0.18, 10: 0.23, 25: 0.30, 50: 0.38, 75: 0.46, 90: 0.54, 95: 0.60}
                )
            },
        ]
    
    def insert_score(
        self,
        score: PolygenicScore,
        distribution: Optional[PopulationDistribution] = None
    ) -> None:
        """
        Insert a polygenic score into the database.
        
        Args:
            score: Polygenic score to insert.
            distribution: Optional population distribution.
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            # Insert score metadata
            cursor.execute("""
                INSERT OR REPLACE INTO polygenic_scores 
                (pgs_id, trait_name, trait_category, publication_doi, publication_year,
                 study_population, sample_size, num_variants, description)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                score.pgs_id, score.trait_name, score.trait_category.value,
                score.publication_doi, score.publication_year, score.study_population,
                score.sample_size, score.num_variants, score.description
            ))
            
            # Insert variants
            for variant in score.variants:
                cursor.execute("""
                    INSERT OR REPLACE INTO pgs_variants
                    (pgs_id, rsid, chromosome, position, effect_allele, other_allele,
                     effect_weight, effect_allele_frequency)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    score.pgs_id, variant.rsid, variant.chromosome, variant.position,
                    variant.effect_allele, variant.other_allele, variant.effect_weight,
                    variant.effect_allele_frequency
                ))
            
            # Insert distribution if provided
            if distribution:
                cursor.execute("""
                    INSERT OR REPLACE INTO population_distributions
                    (pgs_id, population, mean, std, percentiles_json)
                    VALUES (?, ?, ?, ?, ?)
                """, (
                    distribution.pgs_id, distribution.population,
                    distribution.mean, distribution.std,
                    json.dumps(distribution.percentiles)
                ))
            
            conn.commit()
    
    def get_all_scores(self) -> List[PolygenicScore]:
        """
        Get all polygenic scores from the database.
        
        Returns:
            List of PolygenicScore objects (without variants loaded).
        """
        scores = []
        
        with self._get_connection() as conn:
            cursor = conn.cursor()
            # Use columns that exist in the actual database schema
            # Map 'ancestry' to 'study_population', use publication_title as description
            cursor.execute("""
                SELECT pgs_id, trait_name, trait_category, publication_doi,
                       publication_year, ancestry, sample_size, num_variants, publication_title
                FROM polygenic_scores
                WHERE download_status = 'complete'
                ORDER BY trait_name
            """)
            
            for row in cursor.fetchall():
                try:
                    category = TraitCategory(row['trait_category']) if row['trait_category'] else TraitCategory.OTHER
                except ValueError:
                    category = TraitCategory.OTHER
                
                score = PolygenicScore(
                    pgs_id=row['pgs_id'],
                    trait_name=row['trait_name'],
                    trait_category=category,
                    publication_doi=row['publication_doi'],
                    publication_year=row['publication_year'],
                    study_population=row['ancestry'],  # Map ancestry to study_population
                    sample_size=row['sample_size'],
                    num_variants=row['num_variants'],
                    description=row['publication_title'] or '',  # Use publication_title as description
                    variants=[]  # Loaded separately for performance
                )
                scores.append(score)
        
        return scores
    
    def get_score_with_variants(self, pgs_id: str) -> Optional[PolygenicScore]:
        """
        Get a polygenic score with all its variants loaded.
        
        Args:
            pgs_id: Polygenic score identifier.
            
        Returns:
            PolygenicScore with variants, or None if not found.
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            
            # Get score metadata (using actual database schema)
            cursor.execute("""
                SELECT pgs_id, trait_name, trait_category, publication_doi,
                       publication_year, ancestry, sample_size, num_variants, publication_title
                FROM polygenic_scores WHERE pgs_id = ?
            """, (pgs_id,))
            
            row = cursor.fetchone()
            if not row:
                return None
            
            try:
                category = TraitCategory(row['trait_category']) if row['trait_category'] else TraitCategory.OTHER
            except ValueError:
                category = TraitCategory.OTHER
            
            # Get variants (using actual database schema)
            cursor.execute("""
                SELECT rsid, chromosome, position, effect_allele, other_allele,
                       effect_weight, allele_frequency
                FROM pgs_variants WHERE pgs_id = ?
            """, (pgs_id,))
            
            variants = []
            for v_row in cursor.fetchall():
                variant = PolygenicVariant(
                    rsid=v_row['rsid'] or "",
                    chromosome=v_row['chromosome'] or "",
                    position=v_row['position'] or 0,
                    effect_allele=v_row['effect_allele'] or "",
                    other_allele=v_row['other_allele'] or "",
                    effect_weight=v_row['effect_weight'],
                    effect_allele_frequency=v_row['allele_frequency']  # Map column name
                )
                variants.append(variant)
            
            return PolygenicScore(
                pgs_id=row['pgs_id'],
                trait_name=row['trait_name'],
                trait_category=category,
                publication_doi=row['publication_doi'],
                publication_year=row['publication_year'],
                study_population=row['ancestry'],  # Map ancestry to study_population
                sample_size=row['sample_size'],
                num_variants=row['num_variants'],
                description=row['publication_title'] or '',  # Use publication_title as description
                variants=variants
            )
    
    def get_population_distribution(
        self,
        pgs_id: str,
        population: str = "EUR"
    ) -> Optional[PopulationDistribution]:
        """
        Get population distribution for a score.
        
        Args:
            pgs_id: Polygenic score identifier.
            population: Population code (default: EUR).
            
        Returns:
            PopulationDistribution or None if not found.
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT pgs_id, population, mean, std, percentiles_json
                FROM population_distributions
                WHERE pgs_id = ? AND population = ?
            """, (pgs_id, population))
            
            row = cursor.fetchone()
            if not row:
                return None
            
            percentiles = {}
            if row['percentiles_json']:
                try:
                    percentiles = {int(k): v for k, v in json.loads(row['percentiles_json']).items()}
                except (json.JSONDecodeError, ValueError):
                    pass
            
            return PopulationDistribution(
                pgs_id=row['pgs_id'],
                population=row['population'],
                mean=row['mean'],
                std=row['std'],
                percentiles=percentiles
            )
    
    def get_all_distributions(self) -> Dict[str, PopulationDistribution]:
        """
        Get all population distributions.
        
        Returns:
            Dictionary mapping pgs_id to PopulationDistribution.
        """
        distributions = {}
        
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT pgs_id, population, mean, std, percentiles_json
                FROM population_distributions
            """)
            
            for row in cursor.fetchall():
                percentiles = {}
                if row['percentiles_json']:
                    try:
                        percentiles = {int(k): v for k, v in json.loads(row['percentiles_json']).items()}
                    except (json.JSONDecodeError, ValueError):
                        pass
                
                dist = PopulationDistribution(
                    pgs_id=row['pgs_id'],
                    population=row['population'],
                    mean=row['mean'],
                    std=row['std'],
                    percentiles=percentiles
                )
                distributions[row['pgs_id']] = dist
        
        return distributions
    
    def _set_version(self, db_name: str, version: str, source_url: str) -> None:
        """Set version info for a database."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT OR REPLACE INTO database_versions
                (database_name, version, download_date, source_url, record_count)
                VALUES (?, ?, ?, ?, (SELECT COUNT(*) FROM polygenic_scores))
            """, (db_name, version, datetime.now().isoformat(), source_url))
            conn.commit()
    
    def get_version(self, db_name: str) -> Optional[DatabaseVersion]:
        """
        Get version info for a database.
        
        Args:
            db_name: Database name.
            
        Returns:
            DatabaseVersion or None.
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT database_name, version, release_date, download_date,
                       source_url, record_count, checksum
                FROM database_versions WHERE database_name = ?
            """, (db_name,))
            
            row = cursor.fetchone()
            if not row:
                return None
            
            release_date = None
            if row['release_date']:
                try:
                    release_date = datetime.fromisoformat(row['release_date'])
                except ValueError:
                    pass
            
            download_date = datetime.fromisoformat(row['download_date'])
            
            return DatabaseVersion(
                database_name=row['database_name'],
                version=row['version'],
                release_date=release_date,
                download_date=download_date,
                source_url=row['source_url'],
                record_count=row['record_count'],
                checksum=row['checksum']
            )
    
    def get_score_count(self) -> int:
        """Get total number of complete polygenic scores."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM polygenic_scores WHERE download_status = 'complete'")
            return cursor.fetchone()[0]


class DatabaseVersionManager:
    """
    Manages database version checking and updates.
    """
    
    def __init__(self) -> None:
        """Initialize the version manager."""
        self.gwas_db_path = DATABASE_PATH
        self.pgs_db = PolygenicDatabase()
        os.makedirs(BACKUP_DIR, exist_ok=True)
    
    def get_gwas_version(self) -> Optional[DatabaseVersion]:
        """Get GWAS database version info."""
        if not os.path.exists(self.gwas_db_path):
            return None
        
        try:
            conn = sqlite3.connect(self.gwas_db_path)
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT name FROM sqlite_master 
                WHERE type='table' AND name='database_versions'
            """)
            if not cursor.fetchone():
                # Legacy database without version tracking
                cursor.execute("SELECT COUNT(*) FROM gwas_variants")
                count = cursor.fetchone()[0]
                conn.close()
                return DatabaseVersion(
                    database_name="gwas_catalog",
                    version="1.0.0",
                    release_date=None,
                    download_date=datetime.fromtimestamp(os.path.getmtime(self.gwas_db_path)),
                    source_url="https://www.ebi.ac.uk/gwas/",
                    record_count=count
                )
            
            cursor.execute("""
                SELECT * FROM database_versions WHERE database_name = 'gwas_catalog'
            """)
            row = cursor.fetchone()
            conn.close()
            
            if not row:
                return None
            
            return DatabaseVersion(
                database_name=row['database_name'],
                version=row['version'],
                release_date=datetime.fromisoformat(row['release_date']) if row['release_date'] else None,
                download_date=datetime.fromisoformat(row['download_date']),
                source_url=row['source_url'],
                record_count=row['record_count'],
                checksum=row.get('checksum')
            )
        except sqlite3.Error as e:
            logger.error(f"Error getting GWAS version: {e}")
            return None
    
    def get_pgs_version(self) -> Optional[DatabaseVersion]:
        """Get PGS database version info."""
        return self.pgs_db.get_version("pgs_catalog")
    
    def create_backup(self, db_path: str) -> Optional[str]:
        """
        Create a backup of a database.
        
        Args:
            db_path: Path to database to backup.
            
        Returns:
            Path to backup file, or None if failed.
        """
        if not os.path.exists(db_path):
            return None
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        db_name = os.path.basename(db_path)
        backup_path = os.path.join(BACKUP_DIR, f"{db_name}.{timestamp}.bak")
        
        try:
            shutil.copy2(db_path, backup_path)
            logger.info(f"Created backup: {backup_path}")
            return backup_path
        except IOError as e:
            logger.error(f"Failed to create backup: {e}")
            return None
    
    def restore_backup(self, backup_path: str, target_path: str) -> bool:
        """
        Restore a database from backup.
        
        Args:
            backup_path: Path to backup file.
            target_path: Path to restore to.
            
        Returns:
            True if successful.
        """
        try:
            shutil.copy2(backup_path, target_path)
            logger.info(f"Restored backup from {backup_path}")
            return True
        except IOError as e:
            logger.error(f"Failed to restore backup: {e}")
            return False
    
    def list_backups(self) -> List[Tuple[str, datetime]]:
        """
        List available backups.
        
        Returns:
            List of (backup_path, timestamp) tuples.
        """
        backups = []
        if not os.path.exists(BACKUP_DIR):
            return backups
        
        for filename in os.listdir(BACKUP_DIR):
            if filename.endswith('.bak'):
                filepath = os.path.join(BACKUP_DIR, filename)
                mtime = datetime.fromtimestamp(os.path.getmtime(filepath))
                backups.append((filepath, mtime))
        
        return sorted(backups, key=lambda x: x[1], reverse=True)


def get_gwas_database_stats() -> Dict[str, int]:
    """Get statistics about the GWAS database."""
    if not os.path.exists(DATABASE_PATH):
        return {'variants': 0, 'traits': 0, 'genes': 0}
    
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM gwas_variants")
        variants = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(DISTINCT reported_trait) FROM gwas_variants")
        traits = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(DISTINCT mapped_gene) FROM gwas_variants")
        genes = cursor.fetchone()[0]
        
        conn.close()
        return {'variants': variants, 'traits': traits, 'genes': genes}
    except sqlite3.Error:
        return {'variants': 0, 'traits': 0, 'genes': 0}
