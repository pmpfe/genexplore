"""
Data models for the Genetic Analysis Application.

Contains dataclasses for SNP records, GWAS matches, and filter criteria.
"""

from dataclasses import dataclass, field
from typing import Optional, List
import re

from config import (
    VALID_CHROMOSOMES, 
    RSID_PATTERN, 
    GENOTYPE_PATTERN,
    TRAIT_CATEGORIES
)


@dataclass
class SNPRecord:
    """
    Represents a SNP read from a 23andMe raw data file.
    
    Attributes:
        rsid: SNP identifier (e.g., "rs6983267")
        chromosome: Chromosome number or letter (1-22, X, Y, MT)
        position: Genomic position
        genotype: User's genotype at this position (e.g., "GG", "AT")
        
    Raises:
        ValueError: If any field fails validation.
    """
    rsid: str
    chromosome: str
    position: int
    genotype: str
    
    def __post_init__(self) -> None:
        """Validate all fields after initialization."""
        if not RSID_PATTERN.match(self.rsid):
            raise ValueError(f"Invalid RSID format: {self.rsid}")
        
        if self.chromosome not in VALID_CHROMOSOMES:
            raise ValueError(f"Invalid chromosome: {self.chromosome}")
        
        if self.position <= 0:
            raise ValueError(f"Position must be > 0, got: {self.position}")
        
        if not GENOTYPE_PATTERN.match(self.genotype):
            raise ValueError(f"Invalid genotype: {self.genotype}")
    
    def __repr__(self) -> str:
        return f"SNPRecord({self.rsid}, chr{self.chromosome}:{self.position}, {self.genotype})"


@dataclass
class GWASMatch:
    """
    Represents a match between a user SNP and a GWAS Catalog entry.
    
    Attributes:
        rsid: SNP identifier
        chromosome: Chromosome location
        position: Genomic position
        user_genotype: User's genotype at this position
        gene: Mapped gene symbol (may be None)
        trait: Reported trait/disease from GWAS
        risk_allele: Allele associated with higher risk
        p_value: GWAS p-value
        odds_ratio: Effect size (may be None)
        sample_size: Number of subjects in study (may be None)
        category: Trait category
        allele_frequency: Population allele frequency
        impact_score: Calculated impact score (0-10)
    """
    rsid: str
    chromosome: str
    position: int
    user_genotype: str
    gene: Optional[str]
    trait: str
    risk_allele: str
    p_value: float
    odds_ratio: Optional[float]
    sample_size: Optional[int]
    category: str
    allele_frequency: float
    impact_score: float
    
    def __lt__(self, other: 'GWASMatch') -> bool:
        """Enable sorting by impact score (descending)."""
        return self.impact_score > other.impact_score
    
    def __gt__(self, other: 'GWASMatch') -> bool:
        """Enable sorting by impact score (descending)."""
        return self.impact_score < other.impact_score
    
    def __repr__(self) -> str:
        return (f"GWASMatch({self.rsid}, {self.trait[:30]}..., "
                f"score={self.impact_score:.2f})")
    
    def has_risk_allele(self) -> bool:
        """
        Check if user carries the risk allele.
        
        Returns:
            bool: True if user genotype contains the risk allele.
        """
        return self.risk_allele in self.user_genotype
    
    def risk_allele_count(self) -> int:
        """
        Count how many copies of the risk allele the user has.
        
        Returns:
            int: 0, 1, or 2 copies of the risk allele.
        """
        return self.user_genotype.count(self.risk_allele)


@dataclass
class FilterCriteria:
    """
    Stores the current state of UI filters for results.
    
    Attributes:
        min_score: Minimum impact score threshold (0-10)
        max_pvalue: Maximum p-value threshold (0-1)
        category: Trait category filter ('ALL' for no filter)
        search_text: Free text search in traits and genes
        sort_by: Sort column ('score', 'pvalue', 'trait')
        sort_ascending: Sort direction
    """
    min_score: float = 0.0
    max_pvalue: float = 1.0
    category: str = 'ALL'
    search_text: str = ''
    sort_by: str = 'score'
    sort_ascending: bool = False
    
    def __post_init__(self) -> None:
        """Validate filter criteria values."""
        if not 0.0 <= self.min_score <= 10.0:
            raise ValueError(f"min_score must be in [0, 10], got: {self.min_score}")
        
        if not 0.0 <= self.max_pvalue <= 1.0:
            raise ValueError(f"max_pvalue must be in [0, 1], got: {self.max_pvalue}")
        
        valid_sort_columns = {'score', 'pvalue', 'trait', 'gene', 'rsid'}
        if self.sort_by not in valid_sort_columns:
            raise ValueError(f"sort_by must be one of {valid_sort_columns}")
    
    def apply_to_matches(self, matches: List[GWASMatch]) -> List[GWASMatch]:
        """
        Apply all filter criteria to a list of matches.
        
        Args:
            matches: List of GWASMatch objects to filter.
            
        Returns:
            List[GWASMatch]: Filtered and sorted list.
        """
        filtered = matches
        
        # Filter by minimum score
        filtered = [m for m in filtered if m.impact_score >= self.min_score]
        
        # Filter by maximum p-value
        filtered = [m for m in filtered if m.p_value <= self.max_pvalue]
        
        # Filter by category
        if self.category and self.category != 'ALL':
            filtered = [m for m in filtered if m.category == self.category]
        
        # Filter by search text (case-insensitive)
        if self.search_text:
            search_lower = self.search_text.lower()
            filtered = [
                m for m in filtered 
                if (search_lower in m.trait.lower() or 
                    (m.gene and search_lower in m.gene.lower()) or
                    search_lower in m.rsid.lower())
            ]
        
        # Sort results
        if self.sort_by == 'score':
            filtered.sort(key=lambda m: m.impact_score, reverse=not self.sort_ascending)
        elif self.sort_by == 'pvalue':
            filtered.sort(key=lambda m: m.p_value, reverse=not self.sort_ascending)
        elif self.sort_by == 'trait':
            filtered.sort(key=lambda m: m.trait.lower(), reverse=not self.sort_ascending)
        elif self.sort_by == 'gene':
            filtered.sort(key=lambda m: (m.gene or '').lower(), reverse=not self.sort_ascending)
        elif self.sort_by == 'rsid':
            filtered.sort(key=lambda m: m.rsid, reverse=not self.sort_ascending)
        
        return filtered
