"""
Unit tests for the impact score calculation.
"""

import pytest
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.scoring import (
    calculate_impact_score,
    calculate_score_batch,
    get_score_interpretation
)
from config import DEFAULT_IMPACT_SCORE, DEFAULT_ALLELE_FREQUENCY


class TestCalculateImpactScore:
    """Tests for the calculate_impact_score function."""
    
    def test_impact_score_valid_inputs(self):
        """Test impact score calculation with valid inputs."""
        # Highly significant p-value with common variant
        score1 = calculate_impact_score(p_value=1e-10, allele_frequency=0.5)
        assert 0 <= score1 <= 10
        assert score1 > 5  # Should be high impact
        
        # Less significant p-value
        score2 = calculate_impact_score(p_value=0.01, allele_frequency=0.5)
        assert 0 <= score2 <= 10
        
        # More significant should have higher score
        assert score1 > score2
    
    def test_impact_score_edge_cases(self):
        """Test impact score with edge case values."""
        # Very small p-value
        score1 = calculate_impact_score(p_value=1e-100, allele_frequency=0.5)
        assert 0 <= score1 <= 10
        
        # P-value at boundary
        score2 = calculate_impact_score(p_value=1.0, allele_frequency=0.5)
        assert 0 <= score2 <= 10
        
        # AF at boundaries
        score3 = calculate_impact_score(p_value=1e-5, allele_frequency=0.0)
        assert 0 <= score3 <= 10
        
        score4 = calculate_impact_score(p_value=1e-5, allele_frequency=1.0)
        assert 0 <= score4 <= 10
        
        # Rare variant should score higher than common
        assert score3 > score4
    
    def test_impact_score_clamped_to_10(self):
        """Test that impact score is clamped to [0, 10]."""
        # Very extreme values
        score = calculate_impact_score(p_value=1e-200, allele_frequency=0.01)
        assert score <= 10.0
        assert score >= 0.0
        
        # Moderate values
        score2 = calculate_impact_score(p_value=0.5, allele_frequency=0.9)
        assert score2 <= 10.0
        assert score2 >= 0.0
    
    def test_impact_score_with_missing_af(self):
        """Test impact score when allele frequency is None."""
        score = calculate_impact_score(p_value=1e-8, allele_frequency=None)
        assert 0 <= score <= 10
        
        # Should use default AF (0.5)
        score_explicit = calculate_impact_score(
            p_value=1e-8, 
            allele_frequency=DEFAULT_ALLELE_FREQUENCY
        )
        assert score == score_explicit
    
    def test_impact_score_vectorized(self):
        """Test impact score calculation with numpy arrays."""
        p_values = np.array([1e-10, 1e-5, 0.01])
        
        scores = calculate_impact_score(p_value=p_values, allele_frequency=0.5)
        
        assert isinstance(scores, np.ndarray)
        assert len(scores) == 3
        assert all(0 <= s <= 10 for s in scores)
        
        # More significant should have higher score
        assert scores[0] > scores[1] > scores[2]
    
    def test_impact_score_vectorized_with_af_array(self):
        """Test vectorized calculation with AF array."""
        p_values = np.array([1e-8, 1e-8, 1e-8])
        afs = np.array([0.1, 0.5, 0.9])
        
        scores = calculate_impact_score(p_value=p_values, allele_frequency=afs)
        
        assert isinstance(scores, np.ndarray)
        assert len(scores) == 3
        
        # Rare variants (low AF) should have higher scores
        assert scores[0] > scores[1] > scores[2]
    
    def test_invalid_p_value_raises_error(self):
        """Test that invalid p-values raise ValueError."""
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=0, allele_frequency=0.5)
        
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=-0.01, allele_frequency=0.5)
        
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=1.5, allele_frequency=0.5)
    
    def test_invalid_af_raises_error(self):
        """Test that invalid allele frequencies raise ValueError."""
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=1e-5, allele_frequency=-0.1)
        
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=1e-5, allele_frequency=1.5)
    
    def test_invalid_array_values_raise_error(self):
        """Test that invalid array values raise ValueError."""
        p_values = np.array([1e-10, 0, 0.01])
        
        with pytest.raises(ValueError):
            calculate_impact_score(p_value=p_values, allele_frequency=0.5)


class TestCalculateScoreBatch:
    """Tests for the calculate_score_batch function."""
    
    def test_batch_calculation(self):
        """Test batch score calculation."""
        p_values = [1e-10, 1e-5, 0.01]
        afs = [0.5, 0.3, 0.7]
        
        scores = calculate_score_batch(p_values, afs)
        
        assert len(scores) == 3
        assert all(0 <= s <= 10 for s in scores)
    
    def test_batch_with_none_af(self):
        """Test batch calculation with None allele frequencies."""
        p_values = [1e-10, 1e-5, 0.01]
        afs = [0.5, None, 0.7]
        
        scores = calculate_score_batch(p_values, afs)
        
        assert len(scores) == 3
        assert all(0 <= s <= 10 for s in scores)
    
    def test_batch_with_invalid_values(self):
        """Test batch calculation handles invalid values gracefully."""
        p_values = [1e-10, 0, 0.01]  # Zero p-value is invalid
        afs = [0.5, 0.5, 0.5]
        
        scores = calculate_score_batch(p_values, afs)
        
        assert len(scores) == 3
        # Invalid value should get default score
        assert scores[1] == DEFAULT_IMPACT_SCORE


class TestGetScoreInterpretation:
    """Tests for the get_score_interpretation function."""
    
    def test_very_high_impact(self):
        """Test very high impact interpretation."""
        assert "Very High" in get_score_interpretation(10.0)
        assert "Very High" in get_score_interpretation(8.5)
        assert "Very High" in get_score_interpretation(8.0)
    
    def test_high_impact(self):
        """Test high impact interpretation."""
        assert "High Impact" in get_score_interpretation(7.5)
        assert "High Impact" in get_score_interpretation(6.0)
    
    def test_moderate_impact(self):
        """Test moderate impact interpretation."""
        assert "Moderate" in get_score_interpretation(5.5)
        assert "Moderate" in get_score_interpretation(4.0)
    
    def test_low_impact(self):
        """Test low impact interpretation."""
        assert "Low" in get_score_interpretation(3.5)
        assert "Low" in get_score_interpretation(2.0)
    
    def test_minimal_impact(self):
        """Test minimal impact interpretation."""
        assert "Minimal" in get_score_interpretation(1.5)
        assert "Minimal" in get_score_interpretation(0.0)


class TestScoreConsistency:
    """Tests for score calculation consistency."""
    
    def test_p_value_inversely_related(self):
        """Test that lower p-values give higher scores."""
        scores = []
        p_values = [0.1, 0.01, 1e-5, 1e-10, 1e-20]
        
        for p in p_values:
            scores.append(calculate_impact_score(p, 0.5))
        
        # Each score should be higher than the previous
        for i in range(1, len(scores)):
            assert scores[i] >= scores[i-1], f"Score at p={p_values[i]} should be >= score at p={p_values[i-1]}"
    
    def test_af_inversely_related(self):
        """Test that lower AFs (rarer variants) give higher scores."""
        scores = []
        afs = [0.9, 0.5, 0.1, 0.01]
        
        for af in afs:
            scores.append(calculate_impact_score(1e-8, af))
        
        # Each score should be higher than the previous (rarer = higher)
        for i in range(1, len(scores)):
            assert scores[i] >= scores[i-1], f"Score at AF={afs[i]} should be >= score at AF={afs[i-1]}"
    
    def test_reproducibility(self):
        """Test that the same inputs give the same output."""
        for _ in range(10):
            score = calculate_impact_score(p_value=5.2e-11, allele_frequency=0.48)
            assert abs(score - calculate_impact_score(5.2e-11, 0.48)) < 1e-10
