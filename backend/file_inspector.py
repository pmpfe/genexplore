"""
File inspection helpers for genetic data files.
"""

from dataclasses import dataclass, field
import os
from typing import Any, Dict, List, Optional

from backend.parser_common import ParseError
from backend.parser_factory import ParserFactory
from backend.parsers_23andme import Parser23andMe
from backend.parsers_vcf import VCFStatsParser
from models.data_models import SNPRecord
from utils.file_utils import get_file_size, validate_file_exists
from utils.logging_config import get_logger

logger = get_logger(__name__)


@dataclass
class FileInspectionResult:
    """Summary of a genetic file inspection."""

    filepath: str
    detected_format: str
    file_size_bytes: Optional[int]
    is_compressed: bool = False
    stats: Dict[str, Any] = field(default_factory=dict)
    records: Optional[List[SNPRecord]] = None
    analysis_supported: bool = True

    @property
    def filename(self) -> str:
        return os.path.basename(self.filepath)

    @property
    def file_size_human(self) -> str:
        if self.file_size_bytes is None:
            return "Unknown"
        size = float(self.file_size_bytes)
        units = ["B", "KB", "MB", "GB", "TB"]
        unit_index = 0
        while size >= 1024 and unit_index < len(units) - 1:
            size /= 1024.0
            unit_index += 1
        if unit_index == 0:
            return f"{int(size)} {units[unit_index]}"
        return f"{size:.2f} {units[unit_index]}"


def detect_genetic_file_format(filepath: str) -> str:
    """Detect whether a genetic file is 23andMe raw data or VCF."""
    if not validate_file_exists(filepath):
        raise ParseError(f"File not found or not readable: {filepath}")

    return ParserFactory.get_parser(filepath).format_name.lower()


def inspect_genetic_file(filepath: str, progress_callback=None) -> FileInspectionResult:
    """Inspect a genetic file and collect file and genotype statistics."""
    parser = ParserFactory.get_parser(filepath)
    file_format = parser.format_name.lower()
    file_size_bytes = get_file_size(filepath)
    lower_path = str(filepath).lower()
    is_compressed = lower_path.endswith(".gz")

    if file_format == "vcf":
        vcf_result = parser.inspect_file(filepath, progress_callback=progress_callback)
        stats = dict(vcf_result.stats)
        stats.setdefault("analysis_records", len(vcf_result.records))
        stats["genetic_records"] = stats.get("analysis_records", len(vcf_result.records))
        return FileInspectionResult(
            filepath=filepath,
            detected_format="VCF",
            file_size_bytes=file_size_bytes,
            is_compressed=is_compressed,
            stats=stats,
            records=vcf_result.records,
            analysis_supported=True,
        )

    if isinstance(parser, Parser23andMe) and progress_callback is not None:
        def _progress_bridge(lines_processed: int, total_lines: int) -> None:
            if total_lines > 0:
                percent = int((lines_processed / total_lines) * 100)
            else:
                percent = 0
            progress_callback(percent, f"{lines_processed:,}/{total_lines:,} lines")

        parser.set_progress_callback(_progress_bridge)

    records = parser.parse_file(filepath)
    stats = parser.get_parse_stats()
    stats["file_size_bytes"] = file_size_bytes
    stats.setdefault("analysis_records", stats.get("valid_snps", len(records)))
    stats["genetic_records"] = stats.get("analysis_records", len(records))

    return FileInspectionResult(
        filepath=filepath,
        detected_format="23andMe raw",
        file_size_bytes=file_size_bytes,
        is_compressed=is_compressed,
        stats=stats,
        records=records,
        analysis_supported=True,
    )
