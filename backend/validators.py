"""
Input validation functions for the Genetic Analysis Application.
"""

import re
from typing import Tuple, Optional

from config import VALID_CHROMOSOMES, RSID_PATTERN, GENOTYPE_PATTERN, GENOTYPE_UNDETERMINED
from utils.logging_config import get_logger

logger = get_logger(__name__)


def validate_rsid(rsid: str) -> bool:
    """
    Validate that a string is a valid RSID format.
    
    Args:
        rsid: String to validate.
        
    Returns:
        bool: True if valid RSID format (rs followed by digits).
    """
    return bool(RSID_PATTERN.match(rsid))


def validate_chromosome(chromosome: str) -> bool:
    """
    Validate that a chromosome value is valid.
    
    Args:
        chromosome: Chromosome string to validate.
        
    Returns:
        bool: True if valid chromosome (1-22, X, Y, MT).
    """
    return chromosome in VALID_CHROMOSOMES


def validate_position(position: str) -> Tuple[bool, Optional[int]]:
    """
    Validate and parse a genomic position.
    
    Args:
        position: Position string to validate.
        
    Returns:
        Tuple[bool, Optional[int]]: (is_valid, parsed_value or None)
    """
    try:
        pos_int = int(position)
        if pos_int > 0:
            return True, pos_int
        return False, None
    except ValueError:
        return False, None


def validate_genotype(genotype: str) -> Tuple[bool, bool]:
    """
    Validate a genotype string.
    
    Args:
        genotype: Genotype string to validate.
        
    Returns:
        Tuple[bool, bool]: (is_valid_format, is_determined)
            - is_valid_format: True if format is correct
            - is_determined: False if genotype is "--" (undetermined)
    """
    genotype = genotype.strip()
    
    if genotype == GENOTYPE_UNDETERMINED:
        return True, False
    
    if GENOTYPE_PATTERN.match(genotype):
        return True, True
    
    return False, False


def validate_23andme_line(line: str, line_number: int) -> Tuple[bool, Optional[dict], str]:
    """
    Validate a single line from a 23andMe raw data file.
    
    Args:
        line: Line content to validate.
        line_number: Line number for error reporting.
        
    Returns:
        Tuple[bool, Optional[dict], str]: 
            - success: True if line is valid and should be included
            - data: Parsed data dict or None
            - message: Error/info message if applicable
    """
    line = line.strip()
    
    if not line or line.startswith('#'):
        return False, None, ''
    
    parts = line.split()
    
    if len(parts) < 4:
        msg = f"Line {line_number}: Missing fields (expected 4, got {len(parts)})"
        logger.warning(msg)
        return False, None, msg
    
    rsid, chromosome, position, genotype = parts[0], parts[1], parts[2], parts[3]
    
    if not validate_rsid(rsid):
        msg = f"Line {line_number}: Invalid RSID format '{rsid}'"
        logger.warning(msg)
        return False, None, msg
    
    if not validate_chromosome(chromosome):
        msg = f"Line {line_number}: Invalid chromosome '{chromosome}'"
        logger.warning(msg)
        return False, None, msg
    
    is_valid_pos, pos_value = validate_position(position)
    if not is_valid_pos:
        msg = f"Line {line_number}: Invalid position '{position}'"
        logger.warning(msg)
        return False, None, msg
    
    is_valid_geno, is_determined = validate_genotype(genotype)
    if not is_valid_geno:
        msg = f"Line {line_number}: Invalid genotype format '{genotype}'"
        logger.warning(msg)
        return False, None, msg
    
    if not is_determined:
        msg = f"Line {line_number}: Undetermined genotype '{genotype}'"
        logger.debug(msg)
        return False, None, msg
    
    return True, {
        'rsid': rsid,
        'chromosome': chromosome,
        'position': pos_value,
        'genotype': genotype.upper()
    }, ''


def validate_p_value(p_value: float) -> bool:
    """
    Validate a p-value is within acceptable range.
    
    Args:
        p_value: P-value to validate.
        
    Returns:
        bool: True if p_value is in (0, 1].
    """
    return 0 < p_value <= 1


def validate_allele_frequency(af: Optional[float]) -> bool:
    """
    Validate an allele frequency is within acceptable range.
    
    Args:
        af: Allele frequency to validate.
        
    Returns:
        bool: True if af is None or in [0, 1].
    """
    if af is None:
        return True
    return 0 <= af <= 1
