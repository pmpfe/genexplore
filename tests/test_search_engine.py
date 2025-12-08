"""
Unit tests for the search engine.
"""

import pytest
import tempfile
import os
import sys
import sqlite3

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.search_engine import SearchEngine, DatabaseError
from models.data_models import SNPRecord, GWASMatch, FilterCriteria


def create_test_database(db_path: str) -> None:
    """Create a test database with sample data."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE gwas_variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id TEXT NOT NULL,
            chromosome TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref_allele TEXT,
            alt_allele TEXT,
            risk_allele TEXT,
            reported_trait TEXT NOT NULL,
            mapped_gene TEXT,
            p_value REAL NOT NULL,
            odds_ratio REAL,
            sample_size INTEGER,
            category TEXT,
            pubmed_id TEXT,
            UNIQUE(variant_id, reported_trait)
        )
    """)
    
    cursor.execute("""
        CREATE TABLE allele_frequencies (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id TEXT NOT NULL UNIQUE,
            af_overall REAL,
            af_eur REAL,
            af_afr REAL,
            af_eas REAL,
            af_amr REAL
        )
    """)
    
    # Insert test data
    test_variants = [
        ("rs6983267", "8", 128413305, "A", "G", "G", "Colorectal cancer", "MYC", 5.2e-11, 1.27, 201518, "Oncology", "17618284"),
        ("rs7903146", "10", 114758349, "C", "T", "T", "Type 2 diabetes", "TCF7L2", 2.3e-36, 1.37, 120000, "Metabolic", "17293876"),
        ("rs7903146", "10", 114758349, "C", "T", "T", "Fasting glucose", "TCF7L2", 1.1e-20, 1.25, 80000, "Metabolic", "20081858"),
        ("rs1333049", "9", 22125500, "C", "G", "C", "Coronary artery disease", "CDKN2A", 4.8e-14, 1.36, 256000, "Cardiovascular", "17554300"),
        ("rs12913832", "15", 28365618, "A", "G", "G", "Eye color", "HERC2", 1.0e-300, 25.0, 10000, "Physical Trait", "18252222"),
        ("rs1006737", "3", 53127857, "A", "G", "A", "Bipolar disorder", "CACNA1C", 7.0e-8, 1.18, 50000, "Neuropsychiatric", "18711365"),
    ]
    
    cursor.executemany("""
        INSERT INTO gwas_variants 
        (variant_id, chromosome, position, ref_allele, alt_allele, risk_allele,
         reported_trait, mapped_gene, p_value, odds_ratio, sample_size, category, pubmed_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, test_variants)
    
    test_afs = [
        ("rs6983267", 0.48, 0.50, 0.42, 0.38, 0.45),
        ("rs7903146", 0.30, 0.35, 0.28, 0.05, 0.25),
        ("rs1333049", 0.47, 0.50, 0.40, 0.35, 0.45),
        ("rs12913832", 0.78, 0.80, 0.02, 0.05, 0.45),
        ("rs1006737", 0.33, 0.35, 0.28, 0.30, 0.32),
    ]
    
    cursor.executemany("""
        INSERT INTO allele_frequencies 
        (variant_id, af_overall, af_eur, af_afr, af_eas, af_amr)
        VALUES (?, ?, ?, ?, ?, ?)
    """, test_afs)
    
    conn.commit()
    conn.close()


class TestSearchEngine:
    """Tests for the SearchEngine class."""
    
    @pytest.fixture
    def temp_db(self):
        """Create a temporary test database."""
        fd, path = tempfile.mkstemp(suffix='.db')
        os.close(fd)
        create_test_database(path)
        yield path
        os.unlink(path)
    
    @pytest.fixture
    def search_engine(self, temp_db):
        """Create a search engine with test database."""
        return SearchEngine(db_path=temp_db)
    
    def test_verify_database(self, search_engine):
        """Test database verification."""
        assert search_engine.verify_database() is True
    
    def test_verify_database_missing(self):
        """Test verification with missing database."""
        engine = SearchEngine(db_path='/nonexistent/path.db')
        assert engine.verify_database() is False
    
    def test_match_user_snps_found(self, search_engine):
        """Test matching user SNPs that exist in database."""
        user_snps = [
            SNPRecord(rsid="rs6983267", chromosome="8", position=128413305, genotype="GG"),
            SNPRecord(rsid="rs7903146", chromosome="10", position=114758349, genotype="TT"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        assert len(matches) > 0
        
        rsids = [m.rsid for m in matches]
        assert "rs6983267" in rsids
        assert "rs7903146" in rsids
        
        # rs7903146 has multiple traits
        rs7903146_matches = [m for m in matches if m.rsid == "rs7903146"]
        assert len(rs7903146_matches) == 2
    
    def test_match_user_snps_not_found(self, search_engine):
        """Test matching user SNPs that don't exist in database."""
        user_snps = [
            SNPRecord(rsid="rs9999999", chromosome="1", position=100000, genotype="AA"),
            SNPRecord(rsid="rs8888888", chromosome="2", position=200000, genotype="CC"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        assert len(matches) == 0
    
    def test_match_mixed_snps(self, search_engine):
        """Test matching with some SNPs found and some not."""
        user_snps = [
            SNPRecord(rsid="rs6983267", chromosome="8", position=128413305, genotype="GG"),
            SNPRecord(rsid="rs9999999", chromosome="1", position=100000, genotype="AA"),
            SNPRecord(rsid="rs1333049", chromosome="9", position=22125500, genotype="CC"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        rsids = [m.rsid for m in matches]
        assert "rs6983267" in rsids
        assert "rs1333049" in rsids
        assert "rs9999999" not in rsids
    
    def test_match_includes_user_genotype(self, search_engine):
        """Test that matches include the user's genotype."""
        user_snps = [
            SNPRecord(rsid="rs6983267", chromosome="8", position=128413305, genotype="AG"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        assert len(matches) == 1
        assert matches[0].user_genotype == "AG"
    
    def test_match_calculates_impact_score(self, search_engine):
        """Test that matches have calculated impact scores."""
        user_snps = [
            SNPRecord(rsid="rs6983267", chromosome="8", position=128413305, genotype="GG"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        assert len(matches) == 1
        assert 0 <= matches[0].impact_score <= 10
    
    def test_get_database_stats(self, search_engine):
        """Test getting database statistics."""
        stats = search_engine.get_database_stats()
        
        assert 'variants' in stats
        assert 'traits' in stats
        assert 'genes' in stats
        assert stats['variants'] > 0
    
    def test_search_text(self, search_engine):
        """Test text search functionality."""
        user_snps = [
            SNPRecord(rsid="rs6983267", chromosome="8", position=128413305, genotype="GG"),
            SNPRecord(rsid="rs7903146", chromosome="10", position=114758349, genotype="TT"),
            SNPRecord(rsid="rs1333049", chromosome="9", position=22125500, genotype="CC"),
        ]
        
        matches = search_engine.match_user_snps(user_snps)
        
        # Search for diabetes
        filtered = search_engine.search_text("diabetes", matches)
        assert len(filtered) > 0
        assert all("diabetes" in m.trait.lower() for m in filtered)
        
        # Search for gene
        filtered = search_engine.search_text("TCF7L2", matches)
        assert len(filtered) > 0
        assert all(m.gene == "TCF7L2" for m in filtered)


class TestFilterCriteria:
    """Tests for the FilterCriteria class."""
    
    @pytest.fixture
    def sample_matches(self):
        """Create sample GWASMatch objects for testing."""
        return [
            GWASMatch(
                rsid="rs001", chromosome="1", position=100, user_genotype="AA",
                gene="GENE1", trait="Type 2 diabetes", risk_allele="A",
                p_value=1e-20, odds_ratio=1.5, sample_size=10000,
                category="Metabolic", allele_frequency=0.3, impact_score=9.0
            ),
            GWASMatch(
                rsid="rs002", chromosome="2", position=200, user_genotype="GG",
                gene="GENE2", trait="Heart disease", risk_allele="G",
                p_value=1e-10, odds_ratio=1.3, sample_size=20000,
                category="Cardiovascular", allele_frequency=0.5, impact_score=7.5
            ),
            GWASMatch(
                rsid="rs003", chromosome="3", position=300, user_genotype="TT",
                gene="GENE3", trait="Height", risk_allele="T",
                p_value=1e-5, odds_ratio=1.1, sample_size=50000,
                category="Physical Trait", allele_frequency=0.7, impact_score=5.0
            ),
            GWASMatch(
                rsid="rs004", chromosome="4", position=400, user_genotype="CC",
                gene="GENE4", trait="Obesity", risk_allele="C",
                p_value=0.001, odds_ratio=1.2, sample_size=15000,
                category="Metabolic", allele_frequency=0.4, impact_score=3.5
            ),
        ]
    
    def test_filter_by_score(self, sample_matches):
        """Test filtering by minimum impact score."""
        criteria = FilterCriteria(min_score=6.0)
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert len(filtered) == 2
        assert all(m.impact_score >= 6.0 for m in filtered)
    
    def test_filter_by_pvalue(self, sample_matches):
        """Test filtering by maximum p-value."""
        criteria = FilterCriteria(max_pvalue=1e-8)
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert len(filtered) == 2
        assert all(m.p_value <= 1e-8 for m in filtered)
    
    def test_filter_by_category(self, sample_matches):
        """Test filtering by trait category."""
        criteria = FilterCriteria(category="Metabolic")
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert len(filtered) == 2
        assert all(m.category == "Metabolic" for m in filtered)
    
    def test_filter_all_category(self, sample_matches):
        """Test that 'ALL' category returns all matches."""
        criteria = FilterCriteria(category="ALL")
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert len(filtered) == 4
    
    def test_search_text_fuzzy(self, sample_matches):
        """Test fuzzy text search in traits and genes."""
        # Search in trait
        criteria = FilterCriteria(search_text="diabetes")
        filtered = criteria.apply_to_matches(sample_matches)
        assert len(filtered) == 1
        assert filtered[0].trait == "Type 2 diabetes"
        
        # Search in gene
        criteria = FilterCriteria(search_text="GENE2")
        filtered = criteria.apply_to_matches(sample_matches)
        assert len(filtered) == 1
        assert filtered[0].gene == "GENE2"
        
        # Case insensitive search
        criteria = FilterCriteria(search_text="HEART")
        filtered = criteria.apply_to_matches(sample_matches)
        assert len(filtered) == 1
    
    def test_sort_by_score_descending(self, sample_matches):
        """Test sorting by impact score descending."""
        criteria = FilterCriteria(sort_by="score", sort_ascending=False)
        filtered = criteria.apply_to_matches(sample_matches)
        
        scores = [m.impact_score for m in filtered]
        assert scores == sorted(scores, reverse=True)
    
    def test_sort_by_score_ascending(self, sample_matches):
        """Test sorting by impact score ascending."""
        criteria = FilterCriteria(sort_by="score", sort_ascending=True)
        filtered = criteria.apply_to_matches(sample_matches)
        
        scores = [m.impact_score for m in filtered]
        assert scores == sorted(scores)
    
    def test_sort_by_pvalue(self, sample_matches):
        """Test sorting by p-value."""
        criteria = FilterCriteria(sort_by="pvalue", sort_ascending=True)
        filtered = criteria.apply_to_matches(sample_matches)
        
        pvalues = [m.p_value for m in filtered]
        assert pvalues == sorted(pvalues)
    
    def test_sort_by_trait(self, sample_matches):
        """Test sorting by trait name."""
        criteria = FilterCriteria(sort_by="trait", sort_ascending=True)
        filtered = criteria.apply_to_matches(sample_matches)
        
        traits = [m.trait.lower() for m in filtered]
        assert traits == sorted(traits)
    
    def test_combined_filters(self, sample_matches):
        """Test applying multiple filters together."""
        criteria = FilterCriteria(
            min_score=4.0,
            max_pvalue=1e-5,
            category="ALL",
            search_text="",
            sort_by="score",
            sort_ascending=False
        )
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert all(m.impact_score >= 4.0 for m in filtered)
        assert all(m.p_value <= 1e-5 for m in filtered)
    
    def test_filter_returns_empty_list(self, sample_matches):
        """Test that impossible filters return empty list."""
        criteria = FilterCriteria(min_score=9.5)  # Higher than any sample match
        filtered = criteria.apply_to_matches(sample_matches)
        
        assert len(filtered) == 0
    
    def test_invalid_criteria_raises_error(self):
        """Test that invalid criteria values raise errors."""
        with pytest.raises(ValueError):
            FilterCriteria(min_score=15.0)
        
        with pytest.raises(ValueError):
            FilterCriteria(max_pvalue=2.0)
        
        with pytest.raises(ValueError):
            FilterCriteria(sort_by="invalid")


class TestGWASMatch:
    """Tests for the GWASMatch dataclass."""
    
    def test_has_risk_allele(self):
        """Test risk allele detection."""
        match = GWASMatch(
            rsid="rs001", chromosome="1", position=100, user_genotype="AG",
            gene="GENE1", trait="Test", risk_allele="A",
            p_value=1e-10, odds_ratio=1.5, sample_size=10000,
            category="Test", allele_frequency=0.5, impact_score=8.0
        )
        
        assert match.has_risk_allele() is True
        
        match2 = GWASMatch(
            rsid="rs001", chromosome="1", position=100, user_genotype="GG",
            gene="GENE1", trait="Test", risk_allele="A",
            p_value=1e-10, odds_ratio=1.5, sample_size=10000,
            category="Test", allele_frequency=0.5, impact_score=8.0
        )
        
        assert match2.has_risk_allele() is False
    
    def test_risk_allele_count(self):
        """Test counting risk allele copies."""
        # Homozygous risk
        match1 = GWASMatch(
            rsid="rs001", chromosome="1", position=100, user_genotype="AA",
            gene="GENE1", trait="Test", risk_allele="A",
            p_value=1e-10, odds_ratio=1.5, sample_size=10000,
            category="Test", allele_frequency=0.5, impact_score=8.0
        )
        assert match1.risk_allele_count() == 2
        
        # Heterozygous
        match2 = GWASMatch(
            rsid="rs001", chromosome="1", position=100, user_genotype="AG",
            gene="GENE1", trait="Test", risk_allele="A",
            p_value=1e-10, odds_ratio=1.5, sample_size=10000,
            category="Test", allele_frequency=0.5, impact_score=8.0
        )
        assert match2.risk_allele_count() == 1
        
        # No risk allele
        match3 = GWASMatch(
            rsid="rs001", chromosome="1", position=100, user_genotype="GG",
            gene="GENE1", trait="Test", risk_allele="A",
            p_value=1e-10, odds_ratio=1.5, sample_size=10000,
            category="Test", allele_frequency=0.5, impact_score=8.0
        )
        assert match3.risk_allele_count() == 0
    
    def test_sorting_by_impact_score(self):
        """Test that GWASMatch objects can be sorted by impact score."""
        matches = [
            GWASMatch(
                rsid="rs001", chromosome="1", position=100, user_genotype="AA",
                gene="G1", trait="T1", risk_allele="A",
                p_value=1e-10, odds_ratio=1.5, sample_size=10000,
                category="C1", allele_frequency=0.5, impact_score=5.0
            ),
            GWASMatch(
                rsid="rs002", chromosome="2", position=200, user_genotype="GG",
                gene="G2", trait="T2", risk_allele="G",
                p_value=1e-20, odds_ratio=1.8, sample_size=20000,
                category="C2", allele_frequency=0.3, impact_score=9.0
            ),
            GWASMatch(
                rsid="rs003", chromosome="3", position=300, user_genotype="TT",
                gene="G3", trait="T3", risk_allele="T",
                p_value=0.01, odds_ratio=1.1, sample_size=5000,
                category="C3", allele_frequency=0.7, impact_score=3.0
            ),
        ]
        
        sorted_matches = sorted(matches)
        
        # Default sort is descending by impact_score
        assert sorted_matches[0].impact_score == 9.0
        assert sorted_matches[1].impact_score == 5.0
        assert sorted_matches[2].impact_score == 3.0
