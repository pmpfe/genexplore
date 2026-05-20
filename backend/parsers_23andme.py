"""Parser for 23andMe raw data files."""

from collections import Counter
from typing import Callable, List, Optional, Tuple
import gzip
import io
import os

from backend.parser_common import ParseError
from backend.parser_interface import GeneticFileParser
from backend.validators import validate_23andme_line
from models.data_models import SNPRecord
from utils.file_utils import validate_file_exists
from utils.logging_config import get_logger

logger = get_logger(__name__)


class Parser23andMe(GeneticFileParser):
    """
    Parser for 23andMe raw genetic data files.

    Parses tab-separated files with format:
    RSID  CHROMOSOME  POSITION  GENOTYPE
    """

    format_name = '23andme'

    def __init__(self) -> None:
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.warning_details: List[dict] = []
        self.skipped_undetermined: int = 0
        self.invalid_lines: int = 0
        self.comment_lines: int = 0
        self.blank_lines: int = 0
        self.data_lines: int = 0
        self.total_lines: int = 0
        self.valid_lines: int = 0
        self.file_size_bytes: int = 0
        self.chromosome_counts = Counter()
        self.genotype_counts = Counter()
        self.allele_counts = Counter()
        self.homozygous_snps: int = 0
        self.heterozygous_snps: int = 0
        self.progress_callback: Optional[Callable[[int, int], None]] = None

    @classmethod
    def can_parse(cls, filepath: str) -> bool:
        lower_path = str(filepath).lower()
        return lower_path.endswith(('.txt', '.csv', '.txt.gz', '.csv.gz'))

    def set_progress_callback(self, callback) -> None:
        """
        Set a callback function for progress updates.

        The callback should accept (lines_processed, total_lines) parameters.
        """
        self.progress_callback = callback

    def inspect_file(self, filepath: str, progress_callback=None):
        """Inspect a 23andMe file using the existing parsing path."""
        if progress_callback is not None:
            self.set_progress_callback(lambda processed, total: progress_callback(processed, f"{processed:,}/{total:,} lines"))
        records = self.parse_file(filepath)
        return records, self.get_parse_stats()

    def parse_file(self, filepath: str) -> List[SNPRecord]:
        """
        Parse a 23andMe raw data file and extract SNP records.
        """
        self.errors = []
        self.warnings = []
        self.warning_details = []
        self.skipped_undetermined = 0
        self.invalid_lines = 0
        self.comment_lines = 0
        self.blank_lines = 0
        self.data_lines = 0
        self.total_lines = 0
        self.valid_lines = 0
        self.chromosome_counts = Counter()
        self.genotype_counts = Counter()
        self.allele_counts = Counter()
        self.homozygous_snps = 0
        self.heterozygous_snps = 0

        if not validate_file_exists(filepath):
            raise ParseError(f"File not found or not readable: {filepath}")

        self.file_size_bytes = os.path.getsize(filepath)

        try:
            # Support gzipped 23andMe files transparently
            if str(filepath).lower().endswith('.gz'):
                with gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as handle:
                    total_lines = sum(1 for _ in handle)
            else:
                with open(filepath, 'r', encoding='utf-8') as handle:
                    total_lines = sum(1 for _ in handle)
        except OSError as exc:
            raise ParseError(f"Error reading file: {str(exc)}")

        records: List[SNPRecord] = []

        try:
            if str(filepath).lower().endswith('.gz'):
                fh = gzip.open(filepath, 'rt', encoding='utf-8', errors='replace')
            else:
                fh = open(filepath, 'r', encoding='utf-8')

            with fh as handle:
                for line_num, line in enumerate(handle, start=1):
                    self.total_lines += 1

                    stripped_line = line.strip()
                    if not stripped_line:
                        self.blank_lines += 1
                        continue

                    if stripped_line.startswith('#'):
                        self.comment_lines += 1
                        continue

                    self.data_lines += 1

                    if self.progress_callback:
                        self.progress_callback(self.total_lines, total_lines)

                    success, data, message = validate_23andme_line(line, line_num)

                    if not success:
                        if message:
                            if 'Undetermined' in message:
                                self.skipped_undetermined += 1
                            else:
                                self.invalid_lines += 1
                                self.warnings.append(message)
                                self.warning_details.append(
                                    {
                                        'severity': 'warning',
                                        'code': 'invalid_23andme_line',
                                        'line_number': line_num,
                                        'message': message,
                                    }
                                )
                        continue

                    try:
                        genotype = data['genotype'].upper()
                        record = SNPRecord(
                            rsid=data['rsid'],
                            chromosome=data['chromosome'],
                            position=data['position'],
                            genotype=genotype,
                            source_format='23andMe raw',
                            source_metadata={
                                'raw_line_number': line_num,
                                'variant_key': data['rsid'],
                            },
                        )
                        records.append(record)
                        self.valid_lines += 1
                        self.chromosome_counts[record.chromosome] += 1
                        self.genotype_counts[genotype] += 1

                        if len(genotype) == 2 and genotype[0] == genotype[1]:
                            self.homozygous_snps += 1
                        else:
                            self.heterozygous_snps += 1

                        for allele in genotype:
                            if allele in 'ACGT':
                                self.allele_counts[allele] += 1
                    except ValueError as exc:
                        msg = f"Line {line_num}: {str(exc)}"
                        self.invalid_lines += 1
                        self.warnings.append(msg)
                        self.warning_details.append(
                            {
                                'severity': 'warning',
                                'code': 'invalid_23andme_record',
                                'line_number': line_num,
                                'message': msg,
                            }
                        )
                        logger.warning(msg)

        except UnicodeDecodeError as exc:
            raise ParseError(f"File encoding error: {str(exc)}")
        except OSError as exc:
            raise ParseError(f"Error reading file: {str(exc)}")

        if not records:
            raise ParseError("No valid SNP records found in file")

        logger.info(
            f"Parsed {filepath}: {self.valid_lines} valid SNPs, "
            f"{self.skipped_undetermined} undetermined, "
            f"{self.invalid_lines} invalid, "
            f"{len(self.warnings)} warnings"
        )

        return records

    def get_parse_stats(self) -> dict:
        """Get statistics from the last parse operation."""
        return {
            'file_size_bytes': self.file_size_bytes,
            'total_lines': self.total_lines,
            'comment_lines': self.comment_lines,
            'blank_lines': self.blank_lines,
            'data_lines': self.data_lines,
            'valid_snps': self.valid_lines,
            'invalid_lines': self.invalid_lines,
            'skipped_undetermined': self.skipped_undetermined,
            'warnings_count': len(self.warnings),
            'warnings': self.warnings[:10],
            'warning_details': self.warning_details[:10],
            'errors': self.errors,
            'chromosome_counts': dict(self.chromosome_counts),
            'genotype_counts': dict(self.genotype_counts),
            'allele_counts': dict(self.allele_counts),
            'homozygous_snps': self.homozygous_snps,
            'heterozygous_snps': self.heterozygous_snps,
        }


def parse_23andme_file(filepath: str) -> Tuple[List[SNPRecord], dict]:
    """Parse a 23andMe file and return records plus parse statistics."""
    parser = Parser23andMe()
    records = parser.parse_file(filepath)
    stats = parser.get_parse_stats()
    return records, stats
