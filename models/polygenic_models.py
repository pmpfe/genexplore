"""
Data models for Polygenic Risk Score analysis.

Contains dataclasses for polygenic scores, variants, and computation results.
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict
from enum import Enum
from datetime import datetime


class RiskCategory(Enum):
    """Risk level categories for polygenic scores."""
    LOW = "Low"
    INTERMEDIATE = "Intermediate"
    HIGH = "High"
    
    @classmethod
    def from_percentile(cls, percentile: float) -> 'RiskCategory':
        """
        Determine risk category from percentile.
        
        Args:
            percentile: Percentile value (0-100).
            
        Returns:
            RiskCategory: Corresponding risk category.
        """
        if percentile >= 80:
            return cls.HIGH
        elif percentile >= 20:
            return cls.INTERMEDIATE
        else:
            return cls.LOW


class TraitCategory(Enum):
    """Categories for polygenic traits."""
    METABOLIC = "Metabolic"
    CARDIOVASCULAR = "Cardiovascular"
    NEUROPSYCHIATRIC = "Neuropsychiatric"
    ONCOLOGY = "Oncology"
    IMMUNE = "Immune"
    PHYSICAL = "Physical Trait"
    OTHER = "Other"


@dataclass
class PolygenicVariant:
    """
    Represents a single variant within a polygenic score.
    
    Attributes:
        rsid: SNP identifier (e.g., "rs6983267")
        chromosome: Chromosome location
        position: Genomic position
        effect_allele: Allele associated with increased trait/risk
        other_allele: Alternative allele
        effect_weight: Weight/beta coefficient for this variant
        effect_allele_frequency: Population frequency of effect allele
    """
    rsid: str
    chromosome: str
    position: int
    effect_allele: str
    other_allele: str
    effect_weight: float
    effect_allele_frequency: Optional[float] = None
    
    def __repr__(self) -> str:
        return f"PolygenicVariant({self.rsid}, weight={self.effect_weight:.4f})"


@dataclass
class PolygenicScore:
    """
    Represents a Polygenic Risk Score definition.
    
    Attributes:
        pgs_id: Unique identifier for the score (e.g., "PGS000001")
        trait_name: Name of the trait/disease
        trait_category: Category classification
        publication_doi: DOI of the source publication
        publication_year: Year of publication
        study_population: Population used to develop the score
        sample_size: Number of participants in the study
        num_variants: Total number of variants in the score
        variants: List of variants with weights
        description: Optional description of the score
    """
    pgs_id: str
    trait_name: str
    trait_category: TraitCategory
    publication_doi: Optional[str]
    publication_year: Optional[int]
    study_population: str
    sample_size: int
    num_variants: int
    variants: List[PolygenicVariant] = field(default_factory=list)
    description: Optional[str] = None
    
    def __repr__(self) -> str:
        return f"PolygenicScore({self.pgs_id}, {self.trait_name}, {self.num_variants} variants)"


@dataclass
class PopulationDistribution:
    """
    Reference population distribution for score interpretation.
    
    Attributes:
        pgs_id: Associated polygenic score ID
        population: Population name (e.g., "EUR", "AFR")
        mean: Population mean score
        std: Population standard deviation
        percentiles: Dictionary mapping percentile values to scores
    """
    pgs_id: str
    population: str
    mean: float
    std: float
    percentiles: Dict[int, float] = field(default_factory=dict)
    
    def score_to_percentile(self, score: float) -> float:
        """
        Convert a raw score to a percentile.
        
        Uses linear interpolation between known percentile values.
        
        Args:
            score: Raw polygenic score.
            
        Returns:
            float: Estimated percentile (0-100).
        """
        if not self.percentiles:
            # Use z-score approximation if no percentiles available
            if self.std > 0:
                z = (score - self.mean) / self.std
                # Approximate percentile from z-score using normal distribution
                import math
                percentile = 50 * (1 + math.erf(z / math.sqrt(2)))
                return max(0, min(100, percentile))
            return 50.0
        
        # Find surrounding percentiles
        sorted_pcts = sorted(self.percentiles.items())
        
        if score <= sorted_pcts[0][1]:
            return float(sorted_pcts[0][0])
        if score >= sorted_pcts[-1][1]:
            return float(sorted_pcts[-1][0])
        
        for i in range(len(sorted_pcts) - 1):
            pct1, val1 = sorted_pcts[i]
            pct2, val2 = sorted_pcts[i + 1]
            if val1 <= score <= val2:
                # Linear interpolation
                ratio = (score - val1) / (val2 - val1) if val2 != val1 else 0
                return pct1 + ratio * (pct2 - pct1)
        
        return 50.0


@dataclass
class PolygenicResult:
    """
    Result of computing a polygenic score for a user.
    
    Attributes:
        pgs_id: Polygenic score identifier
        trait_name: Name of the trait
        trait_category: Category of the trait
        raw_score: Computed raw score
        normalized_score: Score normalized to population
        percentile: User's percentile in reference population
        risk_category: Risk classification
        variants_found: Number of variants matched in user's genotype
        variants_total: Total variants in the score
        coverage_percent: Percentage of variants found
        population_reference: Reference population used
        variant_contributions: Optional list of (rsid, contribution) tuples
    """
    pgs_id: str
    trait_name: str
    trait_category: TraitCategory
    raw_score: float
    normalized_score: float
    percentile: float
    risk_category: RiskCategory
    variants_found: int
    variants_total: int
    coverage_percent: float
    population_reference: str
    variant_contributions: Optional[List[tuple]] = None
    computation_time_ms: Optional[float] = None
    
    def __repr__(self) -> str:
        return (f"PolygenicResult({self.pgs_id}, {self.trait_name}, "
                f"percentile={self.percentile:.1f}, {self.risk_category.value})")
    
    def is_low_coverage(self) -> bool:
        """Check if variant coverage is below acceptable threshold (80%)."""
        return self.coverage_percent < 80.0
    
    def get_top_contributors(self, n: int = 10) -> List[tuple]:
        """
        Get top N variant contributions to the score.
        
        Args:
            n: Number of top contributors to return.
            
        Returns:
            List of (rsid, contribution) tuples sorted by absolute contribution.
        """
        if not self.variant_contributions:
            return []
        sorted_contribs = sorted(
            self.variant_contributions,
            key=lambda x: abs(x[1]),
            reverse=True
        )
        return sorted_contribs[:n]


@dataclass
class DatabaseVersion:
    """
    Tracks version information for a database.
    
    Attributes:
        database_name: Name of the database (e.g., "gwas_catalog", "pgs_catalog")
        version: Version string
        release_date: Date of the release
        download_date: Date when downloaded
        source_url: URL where data was obtained
        record_count: Number of records in the database
        checksum: MD5/SHA256 checksum for validation
    """
    database_name: str
    version: str
    release_date: Optional[datetime]
    download_date: datetime
    source_url: str
    record_count: int
    checksum: Optional[str] = None
    
    def __repr__(self) -> str:
        return f"DatabaseVersion({self.database_name}, v{self.version})"


@dataclass
class UpdateStatus:
    """
    Status of a database update operation.
    
    Attributes:
        database_name: Name of the database being updated
        current_version: Current installed version
        available_version: Available version for update
        update_available: Whether an update is available
        update_in_progress: Whether update is currently running
        progress_percent: Progress of update (0-100)
        status_message: Current status message
        error_message: Error message if update failed
    """
    database_name: str
    current_version: Optional[str] = None
    available_version: Optional[str] = None
    update_available: bool = False
    update_in_progress: bool = False
    progress_percent: float = 0.0
    status_message: str = ""
    error_message: Optional[str] = None
    
    def is_error(self) -> bool:
        """Check if update is in error state."""
        return self.error_message is not None
