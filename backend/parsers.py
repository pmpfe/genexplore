"""
Parser for 23andMe raw data files.
"""

from typing import List, Tuple
import os

from models.data_models import SNPRecord
from backend.validators import validate_23andme_line
from utils.logging_config import get_logger
from utils.file_utils import validate_file_exists

logger = get_logger(__name__)


class ParseError(Exception):
    """Exception raised when parsing fails completely."""
    pass


class Parser23andMe:
    """
    Parser for 23andMe raw genetic data files.
    
    Parses tab-separated files with format:
    RSID  CHROMOSOME  POSITION  GENOTYPE
    """
    
    def __init__(self) -> None:
        """Initialize the parser."""
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.skipped_undetermined: int = 0
        self.total_lines: int = 0
        self.valid_lines: int = 0
    
    def parse_file(self, filepath: str) -> List[SNPRecord]:
        """
        Parse a 23andMe raw data file and extract SNP records.
        
        Args:
            filepath: Path to the 23andMe raw data file.
            
        Returns:
            List[SNPRecord]: List of valid SNP records.
            
        Raises:
            ParseError: If file cannot be read or is completely invalid.
        """
        self.errors = []
        self.warnings = []
        self.skipped_undetermined = 0
        self.total_lines = 0
        self.valid_lines = 0
        
        if not validate_file_exists(filepath):
            raise ParseError(f"File not found or not readable: {filepath}")
        
        records: List[SNPRecord] = []
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, start=1):
                    self.total_lines += 1
                    
                    success, data, message = validate_23andme_line(line, line_num)
                    
                    if not success:
                        if message:
                            if 'Undetermined' in message:
                                self.skipped_undetermined += 1
                            elif message:
                                self.warnings.append(message)
                        continue
                    
                    try:
                        record = SNPRecord(
                            rsid=data['rsid'],
                            chromosome=data['chromosome'],
                            position=data['position'],
                            genotype=data['genotype']
                        )
                        records.append(record)
                        self.valid_lines += 1
                    except ValueError as e:
                        msg = f"Line {line_num}: {str(e)}"
                        self.warnings.append(msg)
                        logger.warning(msg)
        
        except UnicodeDecodeError as e:
            raise ParseError(f"File encoding error: {str(e)}")
        except IOError as e:
            raise ParseError(f"Error reading file: {str(e)}")
        
        if not records:
            raise ParseError("No valid SNP records found in file")
        
        logger.info(
            f"Parsed {filepath}: {self.valid_lines} valid SNPs, "
            f"{self.skipped_undetermined} undetermined, "
            f"{len(self.warnings)} warnings"
        )
        
        return records
    
    def get_parse_stats(self) -> dict:
        """
        Get statistics from the last parse operation.
        
        Returns:
            dict: Parse statistics including counts and warnings.
        """
        return {
            'total_lines': self.total_lines,
            'valid_snps': self.valid_lines,
            'skipped_undetermined': self.skipped_undetermined,
            'warnings_count': len(self.warnings),
            'warnings': self.warnings[:10],  # First 10 warnings
            'errors': self.errors
        }


def parse_23andme_file(filepath: str) -> Tuple[List[SNPRecord], dict]:
    """
    Convenience function to parse a 23andMe file.
    
    Args:
        filepath: Path to the 23andMe raw data file.
        
    Returns:
        Tuple[List[SNPRecord], dict]: (SNP records, parse statistics)
        
    Raises:
        ParseError: If parsing fails.
    """
    parser = Parser23andMe()
    records = parser.parse_file(filepath)
    stats = parser.get_parse_stats()
    return records, stats
