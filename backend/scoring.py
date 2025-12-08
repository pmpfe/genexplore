"""
Impact score calculation for genetic variants.
"""

import math
from typing import Union, Optional
import numpy as np

from config import DEFAULT_ALLELE_FREQUENCY, DEFAULT_IMPACT_SCORE
from utils.logging_config import get_logger

logger = get_logger(__name__)


def calculate_impact_score(
    p_value: Union[float, np.ndarray],
    allele_frequency: Optional[Union[float, np.ndarray]] = None
) -> Union[float, np.ndarray]:
    """
    Calculate the impact score for a genetic variant.
    
    The impact score combines p-value significance and allele frequency rarity
    into a single 0-10 scale where higher scores indicate more impactful variants.
    
    Formula:
        impact_score = score_p_value + score_allele_frequency
        
        score_p_value = min(10, -log10(p_value) / 10) * 9
        score_allele_frequency = (1 - af) * 3
        
        Final score is clamped to [0, 10]
    
    Args:
        p_value: GWAS p-value (must be in (0, 1]). Can be scalar or numpy array.
        allele_frequency: Population allele frequency [0, 1]. Defaults to 0.5 if None.
        
    Returns:
        Union[float, np.ndarray]: Impact score in range [0.0, 10.0].
        
    Raises:
        ValueError: If p_value is <= 0 or allele_frequency is outside [0, 1].
    """
    is_array = isinstance(p_value, np.ndarray)
    
    if is_array:
        if np.any(p_value <= 0):
            raise ValueError("p_value must be > 0")
        if np.any(p_value > 1):
            raise ValueError("p_value must be <= 1")
    else:
        if p_value <= 0:
            raise ValueError(f"p_value must be > 0, got: {p_value}")
        if p_value > 1:
            raise ValueError(f"p_value must be <= 1, got: {p_value}")
    
    if allele_frequency is None:
        allele_frequency = DEFAULT_ALLELE_FREQUENCY
    
    if isinstance(allele_frequency, np.ndarray):
        if np.any(allele_frequency < 0) or np.any(allele_frequency > 1):
            raise ValueError("allele_frequency must be in [0, 1]")
    else:
        if not 0 <= allele_frequency <= 1:
            raise ValueError(f"allele_frequency must be in [0, 1], got: {allele_frequency}")
    
    try:
        if is_array:
            neg_log_p = -np.log10(p_value)
            score_p_value = np.minimum(neg_log_p / 10, 1.0) * 7.0
        else:
            neg_log_p = -math.log10(p_value)
            score_p_value = min(neg_log_p / 10, 1.0) * 7.0
        
        score_af = (1 - allele_frequency) * 3.0
        
        total_score = score_p_value + score_af
        
        if is_array:
            return np.clip(total_score, 0.0, 10.0)
        else:
            return max(0.0, min(10.0, total_score))
    
    except (ValueError, OverflowError) as e:
        logger.warning(f"Error calculating impact score: {e}. Using default.")
        if is_array:
            return np.full_like(p_value, DEFAULT_IMPACT_SCORE)
        return DEFAULT_IMPACT_SCORE


def calculate_score_batch(
    p_values: list,
    allele_frequencies: list
) -> list:
    """
    Calculate impact scores for a batch of variants.
    
    Args:
        p_values: List of p-values.
        allele_frequencies: List of allele frequencies (can contain None).
        
    Returns:
        list: List of impact scores.
    """
    scores = []
    
    for p_val, af in zip(p_values, allele_frequencies):
        try:
            if af is None:
                af = DEFAULT_ALLELE_FREQUENCY
            score = calculate_impact_score(p_val, af)
            scores.append(score)
        except ValueError as e:
            logger.warning(f"Invalid values for scoring: p={p_val}, af={af}: {e}")
            scores.append(DEFAULT_IMPACT_SCORE)
    
    return scores


def get_score_interpretation(score: float) -> str:
    """
    Get a human-readable interpretation of an impact score.
    
    Args:
        score: Impact score (0-10).
        
    Returns:
        str: Interpretation text.
    """
    if score >= 8.0:
        return "Very High Impact"
    elif score >= 6.0:
        return "High Impact"
    elif score >= 4.0:
        return "Moderate Impact"
    elif score >= 2.0:
        return "Low Impact"
    else:
        return "Minimal Impact"
