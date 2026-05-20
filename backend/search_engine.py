"""
Search engine for matching user SNPs against GWAS database.
"""

import sqlite3
from collections import defaultdict
from typing import List, Optional, Dict, Any, Callable, Tuple
from contextlib import contextmanager

from models.data_models import SNPRecord, GWASMatch, FilterCriteria
from backend.scoring import calculate_impact_score
from config import DATABASE_PATH, DEFAULT_ALLELE_FREQUENCY
from utils.logging_config import get_logger

logger = get_logger(__name__)


class DatabaseError(Exception):
    """Exception raised for database-related errors."""
    pass


class SearchEngine:
    """
    Search engine for matching user SNPs against the GWAS Catalog database.
    
    Provides functionality for:
    - Matching user SNPs to GWAS entries
    - Full-text search in traits and genes
    - Filtering and sorting results
    """
    
    def __init__(self, db_path: Optional[str] = None) -> None:
        """
        Initialize the search engine.
        
        Args:
            db_path: Path to the SQLite database. Uses config default if None.
        """
        self.db_path = db_path or DATABASE_PATH
        self._cache: Dict[str, Any] = {}
    
    @contextmanager
    def _get_connection(self):
        """
        Context manager for database connections.
        
        Yields:
            sqlite3.Connection: Database connection.
            
        Raises:
            DatabaseError: If connection fails.
        """
        try:
            conn = sqlite3.connect(self.db_path)
            conn.row_factory = sqlite3.Row
            yield conn
        except sqlite3.Error as e:
            logger.error(f"Database connection error: {e}")
            raise DatabaseError(f"Failed to connect to database: {e}")
        finally:
            conn.close()
    
    def verify_database(self) -> bool:
        """
        Verify that the database exists and has the required tables.
        
        Returns:
            bool: True if database is valid.
        """
        import os
        if not os.path.exists(self.db_path):
            logger.error(f"Database file not found: {self.db_path}")
            return False
        
        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='gwas_variants'"
                )
                if not cursor.fetchone():
                    logger.error("Missing table: gwas_variants")
                    return False
                
                cursor.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='allele_frequencies'"
                )
                if not cursor.fetchone():
                    logger.error("Missing table: allele_frequencies")
                    return False
                
                return True
        except DatabaseError:
            return False
    
    def match_user_snps(
        self,
        user_snps: List[SNPRecord],
        progress_callback: Optional[Callable[[int, int, str], None]] = None,
    ) -> List[GWASMatch]:
        """
        Match a list of user SNPs against the GWAS database.
        
        Args:
            user_snps: List of SNP records from user's 23andMe file.
            
        Returns:
            List[GWASMatch]: List of GWAS matches found.
            
        Raises:
            DatabaseError: If database query fails.
        """
        matches: List[GWASMatch] = []
        rsid_to_genotype: Dict[str, str] = {}
        positional_to_genotype: Dict[Tuple[str, int], str] = {}

        for snp in user_snps:
            if snp.rsid:
                rsid_to_genotype[snp.rsid] = snp.genotype
            else:
                positional_to_genotype[(snp.chromosome, snp.position)] = snp.genotype

        rsids = list(rsid_to_genotype.keys())
        positions = list(positional_to_genotype.keys())

        logger.info(
            "Searching GWAS database for %s rsIDs and %s chrom/pos fallback keys",
            len(rsids),
            len(positions),
        )

        BATCH_SIZE = 500
        total_batches = ((len(rsids) + BATCH_SIZE - 1) // BATCH_SIZE) + ((len(positions) + BATCH_SIZE - 1) // BATCH_SIZE)
        completed_batches = 0

        def _report_progress(message: str) -> None:
            if progress_callback and total_batches > 0:
                percent = int((completed_batches / total_batches) * 100)
                progress_callback(percent, total_batches, message)
            logger.info(
                "GWAS matching progress: %s/%s batches - %s",
                completed_batches,
                total_batches,
                message,
            )
        
        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                
                for i in range(0, len(rsids), BATCH_SIZE):
                    batch_rsids = rsids[i:i + BATCH_SIZE]
                    completed_batches += 1
                    _report_progress(f"Matching RSIDs {i + 1:,}-{min(i + BATCH_SIZE, len(rsids)):,}...")
                    
                    placeholders = ','.join('?' * len(batch_rsids))
                    query = f"""
                        SELECT 
                            g.variant_id,
                            g.chromosome,
                            g.position,
                            g.risk_allele,
                            g.reported_trait,
                            g.mapped_gene,
                            g.p_value,
                            g.odds_ratio,
                            g.sample_size,
                            g.category,
                            COALESCE(a.af_overall, ?) as af_overall
                        FROM gwas_variants g
                        LEFT JOIN allele_frequencies a ON g.variant_id = a.variant_id
                        WHERE g.variant_id IN ({placeholders})
                    """
                    
                    params = [DEFAULT_ALLELE_FREQUENCY] + batch_rsids
                    cursor.execute(query, params)
                    
                    for row in cursor.fetchall():
                        rsid = row['variant_id']
                        user_genotype = rsid_to_genotype.get(rsid, '')
                        if not user_genotype:
                            continue
                        
                        p_value = row['p_value']
                        af = row['af_overall']
                        
                        try:
                            impact_score = calculate_impact_score(p_value, af)
                        except ValueError:
                            impact_score = 5.0
                        
                        match = GWASMatch(
                            rsid=rsid,
                            chromosome=row['chromosome'],
                            position=row['position'],
                            user_genotype=user_genotype,
                            gene=row['mapped_gene'],
                            trait=row['reported_trait'],
                            risk_allele=row['risk_allele'],
                            p_value=p_value,
                            odds_ratio=row['odds_ratio'],
                            sample_size=row['sample_size'],
                            category=row['category'],
                            allele_frequency=af,
                            impact_score=impact_score
                        )
                        matches.append(match)

                position_seen: set[Tuple[str, int, str]] = set()
                for i in range(0, len(positions), BATCH_SIZE):
                    batch_positions = positions[i:i + BATCH_SIZE]
                    completed_batches += 1
                    _report_progress(f"Matching chrom/pos {i + 1:,}-{min(i + BATCH_SIZE, len(positions)):,}...")

                    if not batch_positions:
                        continue

                    position_clauses = ' OR '.join('(g.chromosome = ? AND g.position = ?)' for _ in batch_positions)
                    query = f"""
                        SELECT
                            g.variant_id,
                            g.chromosome,
                            g.position,
                            g.risk_allele,
                            g.reported_trait,
                            g.mapped_gene,
                            g.p_value,
                            g.odds_ratio,
                            g.sample_size,
                            g.category,
                            COALESCE(a.af_overall, ?) as af_overall
                        FROM gwas_variants g
                        LEFT JOIN allele_frequencies a ON g.variant_id = a.variant_id
                        WHERE {position_clauses}
                    """

                    params: List[Any] = [DEFAULT_ALLELE_FREQUENCY]
                    for chromosome, position in batch_positions:
                        params.extend([chromosome, position])

                    cursor.execute(query, params)

                    for row in cursor.fetchall():
                        lookup_key = (row['variant_id'], row['chromosome'], row['position'])
                        if lookup_key in position_seen:
                            continue
                        position_seen.add(lookup_key)

                        user_genotype = positional_to_genotype.get((row['chromosome'], row['position']), '')
                        if not user_genotype:
                            continue

                        p_value = row['p_value']
                        af = row['af_overall']

                        try:
                            impact_score = calculate_impact_score(p_value, af)
                        except ValueError:
                            impact_score = 5.0

                        match = GWASMatch(
                            rsid=row['variant_id'],
                            chromosome=row['chromosome'],
                            position=row['position'],
                            user_genotype=user_genotype,
                            gene=row['mapped_gene'],
                            trait=row['reported_trait'],
                            risk_allele=row['risk_allele'],
                            p_value=p_value,
                            odds_ratio=row['odds_ratio'],
                            sample_size=row['sample_size'],
                            category=row['category'],
                            allele_frequency=af,
                            impact_score=impact_score
                        )
                        matches.append(match)

                if progress_callback:
                    progress_callback(100, total_batches or 1, f"Matched {len(matches):,} GWAS entries")
                
                logger.info(f"Found {len(matches)} GWAS matches")
                return matches
                
        except sqlite3.Error as e:
            logger.error(f"Database query error: {e}")
            raise DatabaseError(f"Failed to query database: {e}")
    
    def search_text(self, search_term: str, matches: List[GWASMatch]) -> List[GWASMatch]:
        """
        Filter matches by text search in traits and genes.
        
        Uses simple case-insensitive substring matching with fuzzy tolerance.
        
        Args:
            search_term: Text to search for.
            matches: List of GWAS matches to filter.
            
        Returns:
            List[GWASMatch]: Filtered matches.
        """
        if not search_term:
            return matches
        
        search_lower = search_term.lower()
        
        return [
            m for m in matches
            if (search_lower in m.trait.lower() or
                (m.gene and search_lower in m.gene.lower()) or
                search_lower in m.rsid.lower())
        ]
    
    def search_fts(self, search_term: str) -> List[str]:
        """
        Perform full-text search in the database.
        
        Args:
            search_term: Search query.
            
        Returns:
            List[str]: List of matching variant IDs.
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                
                cursor.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='gwas_fts'"
                )
                if cursor.fetchone():
                    query = """
                        SELECT variant_id FROM gwas_fts
                        WHERE gwas_fts MATCH ?
                    """
                    cursor.execute(query, (search_term + '*',))
                else:
                    query = """
                        SELECT variant_id FROM gwas_variants
                        WHERE reported_trait LIKE ? OR mapped_gene LIKE ?
                    """
                    like_term = f'%{search_term}%'
                    cursor.execute(query, (like_term, like_term))
                
                return [row['variant_id'] for row in cursor.fetchall()]
                
        except sqlite3.Error as e:
            logger.warning(f"FTS search error: {e}")
            return []
    
    def filter_matches(
        self,
        matches: List[GWASMatch],
        criteria: FilterCriteria
    ) -> List[GWASMatch]:
        """
        Apply filter criteria to a list of matches.
        
        Args:
            matches: List of GWAS matches.
            criteria: Filter criteria to apply.
            
        Returns:
            List[GWASMatch]: Filtered and sorted matches.
        """
        return criteria.apply_to_matches(matches)
    
    def get_categories(self) -> List[str]:
        """
        Get all unique trait categories from the database.
        
        Returns:
            List[str]: List of category names.
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT DISTINCT category FROM gwas_variants ORDER BY category")
                categories = [row['category'] for row in cursor.fetchall() if row['category']]
                return ['ALL'] + categories
        except sqlite3.Error as e:
            logger.warning(f"Error fetching categories: {e}")
            return ['ALL']
    
    def get_database_stats(self) -> dict:
        """
        Get statistics about the GWAS database.
        
        Returns:
            dict: Database statistics.
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                
                cursor.execute("SELECT COUNT(*) FROM gwas_variants")
                variant_count = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(DISTINCT reported_trait) FROM gwas_variants")
                trait_count = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(DISTINCT mapped_gene) FROM gwas_variants")
                gene_count = cursor.fetchone()[0]
                
                return {
                    'variants': variant_count,
                    'traits': trait_count,
                    'genes': gene_count
                }
        except sqlite3.Error as e:
            logger.error(f"Error getting database stats: {e}")
            return {'variants': 0, 'traits': 0, 'genes': 0}
    
    def clear_cache(self) -> None:
        """Clear the search cache."""
        self._cache.clear()
