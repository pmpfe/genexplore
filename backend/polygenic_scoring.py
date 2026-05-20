"""
Polygenic Risk Score calculation engine.

Provides computation of polygenic risk scores from user genotype data.
"""

import time
from typing import List, Dict, Optional, Tuple, Callable
from dataclasses import dataclass
import math

from models.polygenic_models import (
    PolygenicScore, PolygenicVariant, PolygenicResult,
    PopulationDistribution, RiskCategory, TraitCategory
)
from models.data_models import SNPRecord
from utils.logging_config import get_logger

logger = get_logger(__name__)


class PolygenicScoringError(Exception):
    """Exception raised for polygenic scoring errors."""
    pass


class PolygenicScorer:
    """
    Engine for computing polygenic risk scores.
    
    Handles:
    - Matching user genotypes to score variants
    - Computing weighted sums
    - Normalizing scores against population distributions
    - Handling missing variants gracefully
    """
    
    def __init__(self) -> None:
        """Initialize the scorer."""
        self._genotype_cache: Dict[str, str] = {}
        # Number of original SNPRecord entries loaded (distinct records)
        self._genotype_record_count: int = 0
        self._progress_callback: Optional[Callable[[int, int, str], None]] = None
    
    def set_progress_callback(
        self, 
        callback: Callable[[int, int, str], None]
    ) -> None:
        """
        Set callback for progress updates.
        
        Args:
            callback: Function taking (current, total, message) parameters.
        """
        self._progress_callback = callback
    
    def _report_progress(self, current: int, total: int, message: str) -> None:
        """Report progress if callback is set."""
        if self._progress_callback:
            self._progress_callback(current, total, message)
    
    def load_genotypes(self, snp_records: List[SNPRecord]) -> None:
        """
        Load user genotypes into cache for efficient lookup.
        
        Args:
            snp_records: List of SNP records from user's file.
        """
        # Build a genotype lookup cache using multiple keys per SNP so that
        # different variant identifier styles (rsid, chr:pos, chr:pos:ref>alt)
        # all resolve to the same genotype. This fixes matching for VCF files
        # which may expose variant keys different from the simple locus key
        # used by polygenic definitions.
        self._genotype_cache = {}
        # Track number of actual SNP records provided (tests expect this)
        try:
            self._genotype_record_count = len(snp_records)
        except Exception:
            self._genotype_record_count = 0
        for snp in snp_records:
            genotype = snp.genotype

            # 1) Store by rsid when present (e.g. 'rs12345')
            if getattr(snp, 'rsid', None):
                try:
                    if snp.rsid:
                        self._genotype_cache[snp.rsid] = genotype
                except Exception:
                    pass

            # 2) Store by simple locus key 'chrom:position' (used by PolygenicVariant.locus_key)
            try:
                locus = f"{snp.chromosome}:{snp.position}"
                self._genotype_cache[locus] = genotype
            except Exception:
                pass

            # 3) Store by any parser-provided variant_key (e.g. 'chr:pos:ref>alt')
            try:
                vk = snp.source_metadata.get('variant_key') if snp.source_metadata else None
                if vk:
                    self._genotype_cache[str(vk)] = genotype
            except Exception:
                pass

        logger.info(f"Loaded {len(self._genotype_cache)} genotypes into cache")
    
    def compute_score(
        self,
        pgs: PolygenicScore,
        population_dist: Optional[PopulationDistribution] = None,
        track_contributions: bool = False
    ) -> PolygenicResult:
        """
        Compute a polygenic score for the loaded genotypes.
        
        Args:
            pgs: Polygenic score definition with variants and weights.
            population_dist: Optional population distribution for normalization.
            track_contributions: Whether to track per-variant contributions.
            
        Returns:
            PolygenicResult: Computed result with score and metadata.
            
        Raises:
            PolygenicScoringError: If no genotypes are loaded.
        """
        if not self._genotype_cache:
            raise PolygenicScoringError("No genotypes loaded. Call load_genotypes first.")
        
        start_time = time.perf_counter()
        
        raw_score = 0.0
        variants_found = 0
        contributions = [] if track_contributions else None
        
        for variant in pgs.variants:
            genotype = self._genotype_cache.get(variant.rsid)
            if genotype is None:
                genotype = self._genotype_cache.get(variant.locus_key)
            
            if genotype is None:
                continue
            
            # Count effect alleles in genotype
            effect_count = self._count_effect_alleles(
                genotype, 
                variant.effect_allele,
                variant.other_allele
            )
            
            if effect_count is not None:
                contribution = self._calculate_variant_contribution(variant, effect_count)
                if contribution is None:
                    continue
                raw_score += contribution
                variants_found += 1
                
                if track_contributions:
                    contributions.append((variant.rsid, contribution))
        
        # Compute coverage
        coverage_percent = (variants_found / pgs.num_variants * 100) if pgs.num_variants > 0 else 0.0
        
        # Normalize and compute percentile
        if population_dist:
            normalized_score = self._normalize_score(raw_score, population_dist)
            percentile = population_dist.score_to_percentile(raw_score)
            population_reference = population_dist.population
        else:
            # Estimate distribution from variant weights and allele frequencies
            estimated_dist = self._estimate_population_distribution(pgs)
            if estimated_dist:
                normalized_score = self._normalize_score(raw_score, estimated_dist)
                percentile = estimated_dist.score_to_percentile(raw_score)
                population_reference = "Estimated"
            else:
                normalized_score = raw_score
                percentile = 50.0
                population_reference = "Unknown"
        
        # Determine risk category
        risk_category = RiskCategory.from_percentile(percentile)
        
        computation_time = (time.perf_counter() - start_time) * 1000  # ms
        
        return PolygenicResult(
            pgs_id=pgs.pgs_id,
            trait_name=pgs.trait_name,
            trait_category=pgs.trait_category,
            raw_score=raw_score,
            normalized_score=normalized_score,
            percentile=percentile,
            risk_category=risk_category,
            variants_found=variants_found,
            variants_total=pgs.num_variants,
            coverage_percent=coverage_percent,
            population_reference=population_reference,
            variant_contributions=contributions,
            computation_time_ms=computation_time
        )
    
    def compute_all_scores(
        self,
        scores: List[PolygenicScore],
        population_dists: Dict[str, PopulationDistribution]
    ) -> List[PolygenicResult]:
        """
        Compute all polygenic scores for loaded genotypes.
        
        Args:
            scores: List of polygenic score definitions.
            population_dists: Dictionary mapping pgs_id to population distribution.
            
        Returns:
            List[PolygenicResult]: Results for all scores.
        """
        results = []
        total = len(scores)
        
        for i, pgs in enumerate(scores):
            self._report_progress(i + 1, total, f"Computing {pgs.trait_name}...")
            
            pop_dist = population_dists.get(pgs.pgs_id)
            result = self.compute_score(pgs, pop_dist)
            results.append(result)
        
        logger.info(f"Computed {len(results)} polygenic scores")
        return results
    
    def _count_effect_alleles(
        self,
        genotype: str,
        effect_allele: str,
        other_allele: str
    ) -> Optional[int]:
        """
        Count the number of effect alleles in a genotype.
        
        Args:
            genotype: User's genotype (e.g., "AG", "AA").
            effect_allele: The effect allele to count.
            other_allele: The alternative allele.
            
        Returns:
            Optional[int]: Count of effect alleles (0, 1, or 2), or None if ambiguous.
        """
        if len(genotype) != 2:
            return None
        
        allele1, allele2 = genotype[0], genotype[1]
        valid_alleles = {effect_allele.upper(), other_allele.upper()}
        
        # Check if genotype alleles match expected alleles
        if allele1.upper() not in valid_alleles or allele2.upper() not in valid_alleles:
            # Genotype doesn't match expected alleles - could be strand issue
            # Try complement
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            allele1_comp = complement.get(allele1.upper(), allele1)
            allele2_comp = complement.get(allele2.upper(), allele2)
            
            if allele1_comp in valid_alleles and allele2_comp in valid_alleles:
                allele1, allele2 = allele1_comp, allele2_comp
            else:
                return None
        
        # Count effect alleles
        count = 0
        if allele1.upper() == effect_allele.upper():
            count += 1
        if allele2.upper() == effect_allele.upper():
            count += 1
        
        return count

    def _calculate_variant_contribution(
        self,
        variant: PolygenicVariant,
        effect_count: int,
    ) -> Optional[float]:
        """Calculate a variant contribution using either linear or dosage weights."""
        if effect_count not in (0, 1, 2):
            return None

        if variant.has_dosage_weights:
            dosage_weights = (
                variant.dosage_0_weight,
                variant.dosage_1_weight,
                variant.dosage_2_weight,
            )
            if any(weight is None for weight in dosage_weights):
                return None
            return float(dosage_weights[effect_count])

        return effect_count * variant.effect_weight
    
    def _normalize_score(
        self,
        raw_score: float,
        population_dist: PopulationDistribution
    ) -> float:
        """
        Normalize a raw score using population distribution.
        
        Args:
            raw_score: Raw polygenic score.
            population_dist: Population distribution for normalization.
            
        Returns:
            float: Z-score normalized value.
        """
        if population_dist.std <= 0:
            return 0.0
        return (raw_score - population_dist.mean) / population_dist.std
    
    def _estimate_population_distribution(
        self,
        pgs: PolygenicScore
    ) -> Optional[PopulationDistribution]:
        """
        Estimate population distribution from variant weights and allele frequencies.
        
        Uses the theoretical expectation and variance of polygenic scores under
        Hardy-Weinberg equilibrium. For each variant with effect weight β and
        effect allele frequency p, the expected contribution is 2pβ and the
        variance is 2p(1-p)β².
        
        Args:
            pgs: Polygenic score with variants.
            
        Returns:
            PopulationDistribution with estimated mean and std, or None if cannot estimate.
        """
        if not pgs.variants:
            return None
        
        total_mean = 0.0
        total_variance = 0.0
        variants_with_freq = 0
        
        for variant in pgs.variants:
            # Use effect allele frequency if available, otherwise assume 0.5
            p = variant.effect_allele_frequency if variant.effect_allele_frequency else 0.5
            if variant.has_dosage_weights:
                w0 = variant.dosage_0_weight or 0.0
                w1 = variant.dosage_1_weight or 0.0
                w2 = variant.dosage_2_weight or 0.0
                prob0 = (1 - p) ** 2
                prob1 = 2 * p * (1 - p)
                prob2 = p ** 2
                expected = prob0 * w0 + prob1 * w1 + prob2 * w2
                expected_square = prob0 * (w0 ** 2) + prob1 * (w1 ** 2) + prob2 * (w2 ** 2)
                total_mean += expected
                total_variance += max(0.0, expected_square - expected ** 2)
            else:
                beta = variant.effect_weight
                
                # Expected value: E[score] = 2 * p * β (diploid, two alleles)
                total_mean += 2 * p * beta
                
                # Variance: Var[score] = 2 * p * (1-p) * β² (binomial variance for allele count)
                total_variance += 2 * p * (1 - p) * beta * beta
            
            if variant.effect_allele_frequency:
                variants_with_freq += 1
        
        # Standard deviation
        total_std = math.sqrt(total_variance) if total_variance > 0 else 0.01
        
        # Log coverage of frequency data
        freq_coverage = variants_with_freq / len(pgs.variants) * 100 if pgs.variants else 0
        if freq_coverage < 50:
            logger.debug(f"Low frequency coverage ({freq_coverage:.0f}%) for {pgs.pgs_id}, estimates may be less accurate")
        
        return PopulationDistribution(
            pgs_id=pgs.pgs_id,
            population="Estimated",
            mean=total_mean,
            std=total_std,
            percentiles={}  # Will use z-score approximation
        )
    
    def clear_cache(self) -> None:
        """Clear the genotype cache."""
        self._genotype_cache.clear()
        self._genotype_record_count = 0
    
    @property
    def genotype_count(self) -> int:
        """Get the number of cached genotypes."""
        # Return number of original SNPRecord entries loaded (not internal lookup keys)
        return int(self._genotype_record_count)


def get_risk_interpretation(result: PolygenicResult) -> str:
    """
    Get a human-readable interpretation of a polygenic result.
    
    Args:
        result: Polygenic score result.
        
    Returns:
        str: Interpretation text.
    """
    category = result.risk_category.value
    percentile = result.percentile
    
    if result.is_low_coverage():
        coverage_warning = (
            f"\n\n⚠️ Low coverage warning: Only {result.coverage_percent:.0f}% of variants "
            f"({result.variants_found}/{result.variants_total}) were found in your genotype. "
            "This may reduce the accuracy of this score."
        )
    else:
        coverage_warning = ""
    
    if result.risk_category == RiskCategory.HIGH:
        base_text = (
            f"Your polygenic score for {result.trait_name} places you in the "
            f"{percentile:.0f}th percentile, which is classified as {category} risk. "
            f"This means your genetic predisposition is higher than {percentile:.0f}% of the "
            f"reference population ({result.population_reference})."
        )
    elif result.risk_category == RiskCategory.LOW:
        base_text = (
            f"Your polygenic score for {result.trait_name} places you in the "
            f"{percentile:.0f}th percentile, which is classified as {category} risk. "
            f"This means your genetic predisposition is lower than {100 - percentile:.0f}% of the "
            f"reference population ({result.population_reference})."
        )
    else:
        base_text = (
            f"Your polygenic score for {result.trait_name} places you in the "
            f"{percentile:.0f}th percentile, which is classified as {category} risk. "
            f"This is within the average range for the reference population ({result.population_reference})."
        )
    
    return base_text + coverage_warning


def format_score_summary(result: PolygenicResult) -> Dict[str, str]:
    """
    Format a polygenic result into display-ready strings.
    
    Args:
        result: Polygenic score result.
        
    Returns:
        Dict with formatted values for display.
    """
    return {
        'pgs_id': result.pgs_id,
        'trait': result.trait_name,
        'category': result.trait_category.value,
        'percentile': f"{result.percentile:.1f}%",
        'risk': result.risk_category.value,
        'coverage': f"{result.coverage_percent:.0f}%",
        'raw_score': f"{result.raw_score:.4f}",
        'z_score': f"{result.normalized_score:.2f}",
        'variants': f"{result.variants_found}/{result.variants_total}",
        'population': result.population_reference
    }
