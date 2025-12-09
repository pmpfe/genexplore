"""
Tests for polygenic scoring functionality.
"""

import pytest
from typing import List

from models.polygenic_models import (
    PolygenicScore, PolygenicVariant, PolygenicResult,
    PopulationDistribution, RiskCategory, TraitCategory
)
from models.data_models import SNPRecord
from backend.polygenic_scoring import (
    PolygenicScorer, get_risk_interpretation, format_score_summary
)


class TestRiskCategory:
    """Tests for RiskCategory enum."""
    
    def test_from_percentile_high(self):
        """Test high risk category assignment."""
        assert RiskCategory.from_percentile(95) == RiskCategory.HIGH
        assert RiskCategory.from_percentile(80) == RiskCategory.HIGH
        assert RiskCategory.from_percentile(85) == RiskCategory.HIGH
    
    def test_from_percentile_low(self):
        """Test low risk category assignment."""
        assert RiskCategory.from_percentile(5) == RiskCategory.LOW
        assert RiskCategory.from_percentile(19) == RiskCategory.LOW
        assert RiskCategory.from_percentile(10) == RiskCategory.LOW
    
    def test_from_percentile_intermediate(self):
        """Test intermediate risk category assignment."""
        assert RiskCategory.from_percentile(50) == RiskCategory.INTERMEDIATE
        assert RiskCategory.from_percentile(20) == RiskCategory.INTERMEDIATE
        assert RiskCategory.from_percentile(79) == RiskCategory.INTERMEDIATE


class TestPopulationDistribution:
    """Tests for PopulationDistribution class."""
    
    def test_score_to_percentile_with_percentiles(self):
        """Test percentile calculation using known percentiles."""
        dist = PopulationDistribution(
            pgs_id="TEST001",
            population="EUR",
            mean=1.0,
            std=0.5,
            percentiles={5: 0.2, 25: 0.7, 50: 1.0, 75: 1.3, 95: 1.8}
        )
        
        # Test exact matches
        assert dist.score_to_percentile(1.0) == 50.0
        
        # Test interpolation
        percentile = dist.score_to_percentile(0.85)
        assert 25 < percentile < 50
    
    def test_score_to_percentile_z_score_fallback(self):
        """Test percentile calculation using z-score when no percentiles."""
        dist = PopulationDistribution(
            pgs_id="TEST001",
            population="EUR",
            mean=1.0,
            std=0.5,
            percentiles={}
        )
        
        # Mean should be around 50th percentile
        percentile = dist.score_to_percentile(1.0)
        assert 45 < percentile < 55
        
        # High score should be high percentile
        high_percentile = dist.score_to_percentile(2.0)
        assert high_percentile > 90
        
        # Low score should be low percentile
        low_percentile = dist.score_to_percentile(0.0)
        assert low_percentile < 10


class TestPolygenicVariant:
    """Tests for PolygenicVariant class."""
    
    def test_variant_creation(self):
        """Test basic variant creation."""
        variant = PolygenicVariant(
            rsid="rs7903146",
            chromosome="10",
            position=114758349,
            effect_allele="T",
            other_allele="C",
            effect_weight=0.295,
            effect_allele_frequency=0.30
        )
        
        assert variant.rsid == "rs7903146"
        assert variant.effect_weight == 0.295
        assert variant.effect_allele == "T"


class TestPolygenicScore:
    """Tests for PolygenicScore class."""
    
    def test_score_creation(self):
        """Test basic score creation."""
        variants = [
            PolygenicVariant("rs123", "1", 100, "A", "G", 0.1),
            PolygenicVariant("rs456", "2", 200, "T", "C", 0.2),
        ]
        
        score = PolygenicScore(
            pgs_id="PGS000001",
            trait_name="Test Trait",
            trait_category=TraitCategory.METABOLIC,
            publication_doi="10.1000/test",
            publication_year=2020,
            study_population="EUR",
            sample_size=100000,
            num_variants=2,
            variants=variants
        )
        
        assert score.pgs_id == "PGS000001"
        assert len(score.variants) == 2
        assert score.trait_category == TraitCategory.METABOLIC


class TestPolygenicResult:
    """Tests for PolygenicResult class."""
    
    def test_is_low_coverage(self):
        """Test low coverage detection."""
        result_low = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.0,
            normalized_score=0.0,
            percentile=50.0,
            risk_category=RiskCategory.INTERMEDIATE,
            variants_found=70,
            variants_total=100,
            coverage_percent=70.0,
            population_reference="EUR"
        )
        assert result_low.is_low_coverage()
        
        result_ok = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.0,
            normalized_score=0.0,
            percentile=50.0,
            risk_category=RiskCategory.INTERMEDIATE,
            variants_found=85,
            variants_total=100,
            coverage_percent=85.0,
            population_reference="EUR"
        )
        assert not result_ok.is_low_coverage()
    
    def test_get_top_contributors(self):
        """Test getting top contributing variants."""
        result = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.0,
            normalized_score=0.0,
            percentile=50.0,
            risk_category=RiskCategory.INTERMEDIATE,
            variants_found=5,
            variants_total=10,
            coverage_percent=50.0,
            population_reference="EUR",
            variant_contributions=[
                ("rs1", 0.1),
                ("rs2", -0.3),
                ("rs3", 0.2),
                ("rs4", 0.05),
                ("rs5", -0.15),
            ]
        )
        
        top3 = result.get_top_contributors(3)
        assert len(top3) == 3
        # Should be sorted by absolute value descending
        assert top3[0] == ("rs2", -0.3)
        assert top3[1] == ("rs3", 0.2)
        assert top3[2] == ("rs5", -0.15)


class TestPolygenicScorer:
    """Tests for PolygenicScorer class."""
    
    @pytest.fixture
    def sample_genotypes(self) -> List[SNPRecord]:
        """Create sample genotype data."""
        return [
            SNPRecord("rs7903146", "10", 114758349, "TT"),
            SNPRecord("rs1801282", "3", 12393125, "CG"),
            SNPRecord("rs5219", "11", 17409572, "CC"),
        ]
    
    @pytest.fixture
    def sample_score(self) -> PolygenicScore:
        """Create sample polygenic score."""
        return PolygenicScore(
            pgs_id="PGS_TEST",
            trait_name="Test Diabetes",
            trait_category=TraitCategory.METABOLIC,
            publication_doi=None,
            publication_year=2020,
            study_population="EUR",
            sample_size=100000,
            num_variants=3,
            variants=[
                PolygenicVariant("rs7903146", "10", 114758349, "T", "C", 0.3, 0.30),
                PolygenicVariant("rs1801282", "3", 12393125, "C", "G", 0.15, 0.12),
                PolygenicVariant("rs5219", "11", 17409572, "T", "C", 0.1, 0.35),
            ]
        )
    
    @pytest.fixture
    def sample_distribution(self) -> PopulationDistribution:
        """Create sample population distribution."""
        return PopulationDistribution(
            pgs_id="PGS_TEST",
            population="EUR",
            mean=0.5,
            std=0.2,
            percentiles={5: 0.17, 25: 0.37, 50: 0.5, 75: 0.63, 95: 0.83}
        )
    
    def test_load_genotypes(self, sample_genotypes):
        """Test loading genotypes into scorer."""
        scorer = PolygenicScorer()
        scorer.load_genotypes(sample_genotypes)
        
        assert scorer.genotype_count == 3
    
    def test_compute_score_no_genotypes(self, sample_score):
        """Test error when computing without loaded genotypes."""
        scorer = PolygenicScorer()
        
        with pytest.raises(Exception):
            scorer.compute_score(sample_score)
    
    def test_compute_score(self, sample_genotypes, sample_score, sample_distribution):
        """Test basic score computation."""
        scorer = PolygenicScorer()
        scorer.load_genotypes(sample_genotypes)
        
        result = scorer.compute_score(sample_score, sample_distribution)
        
        assert result.pgs_id == "PGS_TEST"
        assert result.trait_name == "Test Diabetes"
        assert result.variants_total == 3
        assert result.variants_found > 0
        assert 0 <= result.percentile <= 100
        assert result.risk_category in [RiskCategory.LOW, RiskCategory.INTERMEDIATE, RiskCategory.HIGH]
    
    def test_compute_score_with_contributions(self, sample_genotypes, sample_score, sample_distribution):
        """Test score computation with contribution tracking."""
        scorer = PolygenicScorer()
        scorer.load_genotypes(sample_genotypes)
        
        result = scorer.compute_score(sample_score, sample_distribution, track_contributions=True)
        
        assert result.variant_contributions is not None
        assert len(result.variant_contributions) > 0
    
    def test_count_effect_alleles(self):
        """Test effect allele counting."""
        scorer = PolygenicScorer()
        
        # Homozygous for effect allele
        assert scorer._count_effect_alleles("TT", "T", "C") == 2
        
        # Heterozygous
        assert scorer._count_effect_alleles("TC", "T", "C") == 1
        
        # Homozygous for other allele
        assert scorer._count_effect_alleles("CC", "T", "C") == 0
        
        # Complementary strand
        assert scorer._count_effect_alleles("AA", "T", "C") == 2  # A complements T
    
    def test_clear_cache(self, sample_genotypes):
        """Test cache clearing."""
        scorer = PolygenicScorer()
        scorer.load_genotypes(sample_genotypes)
        assert scorer.genotype_count == 3
        
        scorer.clear_cache()
        assert scorer.genotype_count == 0


class TestInterpretationFunctions:
    """Tests for interpretation helper functions."""
    
    def test_get_risk_interpretation_high(self):
        """Test interpretation for high risk."""
        result = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test Disease",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.0,
            normalized_score=2.0,
            percentile=90.0,
            risk_category=RiskCategory.HIGH,
            variants_found=100,
            variants_total=100,
            coverage_percent=100.0,
            population_reference="EUR"
        )
        
        interpretation = get_risk_interpretation(result)
        
        assert "90th percentile" in interpretation
        assert "High" in interpretation
        assert "higher" in interpretation.lower()
    
    def test_get_risk_interpretation_low_coverage(self):
        """Test interpretation includes low coverage warning."""
        result = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test Disease",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.0,
            normalized_score=0.0,
            percentile=50.0,
            risk_category=RiskCategory.INTERMEDIATE,
            variants_found=50,
            variants_total=100,
            coverage_percent=50.0,
            population_reference="EUR"
        )
        
        interpretation = get_risk_interpretation(result)
        
        assert "coverage" in interpretation.lower() or "warning" in interpretation.lower()
    
    def test_format_score_summary(self):
        """Test score summary formatting."""
        result = PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test Disease",
            trait_category=TraitCategory.METABOLIC,
            raw_score=1.234567,
            normalized_score=0.5,
            percentile=75.0,
            risk_category=RiskCategory.INTERMEDIATE,
            variants_found=90,
            variants_total=100,
            coverage_percent=90.0,
            population_reference="EUR"
        )
        
        summary = format_score_summary(result)
        
        assert summary['pgs_id'] == "PGS000001"
        assert summary['percentile'] == "75.0%"
        assert summary['risk'] == "Intermediate"
        assert summary['coverage'] == "90.0%"
        assert "1.2346" in summary['raw_score']
