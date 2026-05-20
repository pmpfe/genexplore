"""Parser selection helpers for supported genetic file formats."""

from typing import Type
import gzip
import io

from backend.parser_common import ParseError
from backend.parser_interface import GeneticFileParser
from backend.parsers_23andme import Parser23andMe
from backend.parsers_vcf import VCFStatsParser
from utils.file_utils import validate_file_exists


class FileFormatDetector:
    """Detect the most likely parser for a genetic input file."""

    @staticmethod
    def detect(filepath: str) -> str:
        if not validate_file_exists(filepath):
            raise ParseError(f'File not found or not readable: {filepath}')

        if VCFStatsParser.can_parse(filepath):
            return 'vcf'

        # Attempt to read the first meaningful lines safely, supporting gzipped files.
        is_gzip = False
        try:
            with open(filepath, 'rb') as fh:
                magic = fh.read(2)
                is_gzip = (magic == b'\x1f\x8b')
        except OSError as exc:
            raise ParseError(f'Error reading file: {str(exc)}')

        try:
            if is_gzip:
                with gzip.open(filepath, 'rt', encoding='utf-8', errors='replace') as handle:
                    for line in handle:
                        stripped = line.strip()
                        if not stripped:
                            continue
                        if stripped.startswith('##fileformat=VCF') or stripped.startswith('#CHROM'):
                            return 'vcf'
                        if stripped.startswith('#'):
                            continue
                        parts = stripped.split()
                        if len(parts) >= 4:
                            return '23andme'
                        break
            else:
                with open(filepath, 'r', encoding='utf-8', errors='replace') as handle:
                    for line in handle:
                        stripped = line.strip()
                        if not stripped:
                            continue
                        if stripped.startswith('##fileformat=VCF') or stripped.startswith('#CHROM'):
                            return 'vcf'
                        if stripped.startswith('#'):
                            continue
                        parts = stripped.split()
                        if len(parts) >= 4:
                            return '23andme'
                        break
        except OSError as exc:
            raise ParseError(f'Error reading file: {str(exc)}')

        if Parser23andMe.can_parse(filepath):
            return '23andme'

        raise ParseError(f'Unsupported genetic file format: {filepath}')


class ParserFactory:
    """Factory for creating the parser that best matches a file."""

    _parsers: dict[str, Type[GeneticFileParser]] = {
        'vcf': VCFStatsParser,
        '23andme': Parser23andMe,
    }

    @classmethod
    def get_parser(cls, filepath: str) -> GeneticFileParser:
        format_name = FileFormatDetector.detect(filepath)
        parser_class = cls._parsers.get(format_name)
        if parser_class is None:
            raise ParseError(f'No parser available for file: {filepath}')
        return parser_class()