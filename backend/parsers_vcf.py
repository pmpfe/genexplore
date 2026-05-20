"""Streaming parser and inspector for VCF and VCF.gz files."""

from collections import Counter
from dataclasses import dataclass, field
import gzip
import io
from typing import Any, Callable, Dict, List, Optional

from backend.parser_common import ParseError, ParserIssue
from backend.parser_interface import GeneticFileParser
from config import RSID_PATTERN, VALID_CHROMOSOMES
from models.data_models import SNPRecord
from utils.file_utils import get_file_size, validate_file_exists
from utils.logging_config import get_logger

logger = get_logger(__name__)


@dataclass
class VCFInspectionResult:
    """Aggregated statistics and normalized records for a VCF file."""

    filepath: str
    file_size_bytes: Optional[int]
    is_compressed: bool
    stats: Dict[str, Any] = field(default_factory=dict)
    records: List[SNPRecord] = field(default_factory=list)

    @property
    def file_size_human(self) -> str:
        return _format_file_size(self.file_size_bytes)


class VCFStatsParser(GeneticFileParser):
    """Streaming parser that gathers structure and genotype statistics from VCF files."""

    format_name = 'VCF'

    def __init__(self) -> None:
        self.reset()

    @classmethod
    def can_parse(cls, filepath: str) -> bool:
        lower_path = str(filepath).lower()
        return lower_path.endswith(('.vcf', '.vcf.gz'))

    def reset(self) -> None:
        self.total_lines = 0
        self.blank_lines = 0
        self.metadata_lines = 0
        self.header_lines = 0
        self.data_lines = 0
        self.sample_count = 0
        self.samples: List[str] = []
        self.analysis_sample: Optional[str] = None
        self.column_count = 0
        self.vcf_version: Optional[str] = None
        self.chromosome_counts = Counter()
        self.variant_types = Counter()
        self.genotype_counts = Counter()
        self.filter_counts = Counter()
        self.info_key_counts = Counter()
        self.called_genotypes = 0
        self.missing_genotypes = 0
        self.haploid_genotypes = 0
        self.homozygous_reference = 0
        self.homozygous_alternate = 0
        self.heterozygous = 0
        self.complex_genotypes = 0
        self.phased_genotypes = 0
        self.unphased_genotypes = 0
        self.missing_id_count = 0
        self.qual_known_count = 0
        self.snv_count = 0
        self.indel_count = 0
        self.multiallelic_count = 0
        self.transition_snps = 0
        self.transversion_snps = 0
        self.analysis_records = 0
        self.skipped_missing_rsid = 0
        self.skipped_missing_sample = 0
        self.skipped_invalid_chromosome = 0
        self.skipped_non_snp = 0
        self.skipped_invalid_genotype = 0
        self.warning_details: List[ParserIssue] = []
        self.warnings: List[str] = []

    def inspect_file(
        self,
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> VCFInspectionResult:
        """Inspect a VCF or VCF.gz file and return aggregated stats plus records."""
        self.reset()
        records = self._scan_file(filepath, progress_callback=progress_callback)
        file_size_bytes = get_file_size(filepath)
        is_compressed = self._is_gzip_file(filepath)

        stats = self._build_stats(file_size_bytes)
        logger.info(
            "Inspected VCF %s: %s analyzable records, %s samples, %s warnings",
            filepath,
            self.analysis_records,
            self.sample_count,
            len(self.warnings),
        )

        return VCFInspectionResult(
            filepath=filepath,
            file_size_bytes=file_size_bytes,
            is_compressed=is_compressed,
            stats=stats,
            records=records,
        )

    def parse_file(
        self,
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> List[SNPRecord]:
        """Return normalized records from a VCF file."""
        return self.inspect_file(filepath, progress_callback=progress_callback).records

    def _scan_file(
        self,
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> List[SNPRecord]:
        if not validate_file_exists(filepath):
            raise ParseError(f"File not found or not readable: {filepath}")

        is_compressed = self._is_gzip_file(filepath)
        file_size_bytes = get_file_size(filepath)
        records: List[SNPRecord] = []

        def process_stream(handle: io.TextIOBase, raw_handle) -> None:
            for line_number, line in enumerate(handle, start=1):
                self.total_lines += 1
                if progress_callback and self.total_lines % 10000 == 0 and file_size_bytes:
                    current_bytes = raw_handle.tell()
                    progress_percent = min(99, int((current_bytes / file_size_bytes) * 100))
                    progress_callback(progress_percent, f"Inspected {self.total_lines:,} lines")

                stripped = line.strip()
                if not stripped:
                    self.blank_lines += 1
                    continue

                if stripped.startswith('##'):
                    self.metadata_lines += 1
                    if stripped.startswith('##fileformat='):
                        self.vcf_version = stripped.split('=', 1)[1].strip()
                    continue

                if stripped.startswith('#'):
                    self.header_lines += 1
                    columns = stripped.lstrip('#').split('\t')
                    self.column_count = len(columns)
                    if len(columns) >= 9:
                        self.samples = columns[9:]
                        self.sample_count = len(self.samples)
                        self.analysis_sample = self.samples[0] if self.samples else None
                    continue

                self.data_lines += 1
                fields = stripped.split('\t')
                if len(fields) < 8:
                    self._warn(
                        code='invalid_vcf_record',
                        line_number=line_number,
                        message='invalid VCF record (expected at least 8 columns)',
                    )
                    continue

                chrom, position, identifier, ref, alt, qual, filt, info = fields[:8]
                format_field = fields[8] if len(fields) > 8 else ''
                sample_fields = fields[9:] if len(fields) > 9 else []

                variant_type = self._record_variant_type(ref, alt)
                self._record_metadata(chrom, identifier, qual, filt, info)

                if variant_type != 'snp':
                    self.skipped_non_snp += 1
                    continue

                if not sample_fields:
                    self.skipped_missing_sample += 1
                    self._warn(
                        code='missing_sample',
                        line_number=line_number,
                        message='no sample genotype columns found',
                        details={'chromosome': chrom, 'position': position},
                    )
                    continue

                raw_gt = sample_fields[0].split(':', 1)[0].strip()
                canonical_genotype = self._normalize_vcf_genotype(raw_gt, ref, alt)
                if canonical_genotype is None:
                    self.skipped_invalid_genotype += 1
                    continue

                self._record_genotype(canonical_genotype, raw_gt, ref)
                self.called_genotypes += 1

                normalized_chromosome = self._normalize_chromosome(chrom)
                if normalized_chromosome is None:
                    self.skipped_invalid_chromosome += 1
                    self._warn(
                        code='invalid_chromosome',
                        line_number=line_number,
                        message=f'invalid chromosome value: {chrom}',
                    )
                    continue

                rsid = self._extract_rsid(identifier)
                variant_key = self._build_variant_key(normalized_chromosome, position, ref, alt, rsid)

                try:
                    position_value = int(position)
                    record = SNPRecord(
                        rsid=rsid,
                        chromosome=normalized_chromosome,
                        position=position_value,
                        genotype=canonical_genotype,
                        source_format='VCF',
                        source_metadata={
                            'variant_key': variant_key,
                            'ref': ref,
                            'alt': alt,
                            'qual': qual,
                            'filter': filt,
                            'info': info,
                            'format': format_field,
                            'raw_gt': raw_gt,
                            'sample_name': self.analysis_sample,
                            'vcf_identifier': identifier,
                        },
                    )
                except ValueError as exc:
                    self.skipped_invalid_genotype += 1
                    self._warn(
                        code='invalid_record',
                        line_number=line_number,
                        message=str(exc),
                    )
                    continue

                records.append(record)
                self.analysis_records += 1
                self.chromosome_counts[record.chromosome] += 1
                self.genotype_counts[record.genotype] += 1

        try:
            with open(filepath, 'rb') as raw_handle:
                if is_compressed:
                    with gzip.GzipFile(fileobj=raw_handle, mode='rb') as binary_handle:
                        with io.TextIOWrapper(binary_handle, encoding='utf-8', errors='replace') as handle:
                            process_stream(handle, raw_handle)
                else:
                    with io.TextIOWrapper(raw_handle, encoding='utf-8', errors='replace') as handle:
                        process_stream(handle, raw_handle)

            if progress_callback and file_size_bytes:
                progress_callback(100, f"Inspected {self.total_lines:,} lines")

        except UnicodeDecodeError as exc:
            raise ParseError(f"File encoding error: {str(exc)}")
        except OSError as exc:
            raise ParseError(f"Error reading file: {str(exc)}")

        if not records:
            raise ParseError('No analyzable SNP records found in VCF file')

        return records

    def _build_stats(self, file_size_bytes: Optional[int]) -> Dict[str, Any]:
        return {
            'file_size_bytes': file_size_bytes,
            'total_lines': self.total_lines,
            'blank_lines': self.blank_lines,
            'metadata_lines': self.metadata_lines,
            'header_lines': self.header_lines,
            'data_lines': self.data_lines,
            'record_count': self.data_lines,
            'analysis_records': self.analysis_records,
            'sample_count': self.sample_count,
            'samples': self.samples,
            'analysis_sample': self.analysis_sample,
            'column_count': self.column_count,
            'vcf_version': self.vcf_version,
            'chromosome_counts': dict(self.chromosome_counts),
            'variant_types': dict(self.variant_types),
            'genotype_counts': dict(self.genotype_counts),
            'filter_counts': dict(self.filter_counts),
            'info_key_counts': dict(self.info_key_counts),
            'called_genotypes': self.called_genotypes,
            'missing_genotypes': self.missing_genotypes,
            'haploid_genotypes': self.haploid_genotypes,
            'homozygous_reference_genotypes': self.homozygous_reference,
            'homozygous_alternate_genotypes': self.homozygous_alternate,
            'heterozygous_genotypes': self.heterozygous,
            'complex_genotypes': self.complex_genotypes,
            'phased_genotypes': self.phased_genotypes,
            'unphased_genotypes': self.unphased_genotypes,
            'missing_id_count': self.missing_id_count,
            'qual_known_count': self.qual_known_count,
            'snv_count': self.snv_count,
            'indel_count': self.indel_count,
            'multiallelic_count': self.multiallelic_count,
            'transition_snps': self.transition_snps,
            'transversion_snps': self.transversion_snps,
            'skipped_missing_rsid': self.skipped_missing_rsid,
            'skipped_missing_sample': self.skipped_missing_sample,
            'skipped_invalid_chromosome': self.skipped_invalid_chromosome,
            'skipped_non_snp': self.skipped_non_snp,
            'skipped_invalid_genotype': self.skipped_invalid_genotype,
            'warnings_count': len(self.warnings),
            'warnings': self.warnings[:10],
            'warning_details': [
                {
                    'severity': issue.severity,
                    'code': issue.code,
                    'message': issue.message,
                    'line_number': issue.line_number,
                    'details': issue.details,
                }
                for issue in self.warning_details[:10]
            ],
        }

    def _record_variant_type(self, ref: str, alt_field: str) -> str:
        variant_type = self._classify_variant(ref, alt_field)
        self.variant_types[variant_type] += 1

        if variant_type == 'snp':
            self.snv_count += 1
            if self._is_transition(ref, alt_field):
                self.transition_snps += 1
            else:
                self.transversion_snps += 1
        elif variant_type == 'indel':
            self.indel_count += 1
        elif variant_type == 'multiallelic':
            self.multiallelic_count += 1

        return variant_type

    def _record_metadata(self, chrom: str, identifier: str, qual: str, filt: str, info: str) -> None:
        if identifier == '.':
            self.missing_id_count += 1

        if qual != '.':
            self.qual_known_count += 1
        self.filter_counts[filt or '.'] += 1

        for key in self._extract_info_keys(info):
            self.info_key_counts[key] += 1

    def _warn(self, code: str, line_number: int, message: str, details: Optional[Dict[str, Any]] = None) -> None:
        warning = f"Line {line_number}: {message}"
        self.warnings.append(warning)
        self.warning_details.append(
            ParserIssue(
                severity='warning',
                code=code,
                message=message,
                line_number=line_number,
                details=details or {},
            )
        )
        logger.warning(warning)

    def _is_gzip_file(self, filepath: str) -> bool:
        if str(filepath).lower().endswith('.gz'):
            return True

        try:
            with open(filepath, 'rb') as handle:
                return handle.read(2) == b'\x1f\x8b'
        except OSError:
            return False

    def _normalize_chromosome(self, chromosome: str) -> Optional[str]:
        value = chromosome.strip()
        if value.lower().startswith('chr'):
            value = value[3:]

        value = value.upper()
        if value == 'M':
            value = 'MT'

        if value in VALID_CHROMOSOMES:
            return value

        return None

    def _extract_rsid(self, identifier: str) -> Optional[str]:
        if not identifier or identifier == '.':
            return None

        for token in identifier.split(';'):
            token = token.strip()
            if RSID_PATTERN.match(token):
                return token

        return None

    def _build_variant_key(
        self,
        chromosome: str,
        position: str,
        ref: str,
        alt_field: str,
        rsid: Optional[str],
    ) -> str:
        if rsid:
            return rsid
        return f'{chromosome}:{position}:{ref}>{alt_field}'

    def _classify_variant(self, ref: str, alt_field: str) -> str:
        alts = [alt for alt in alt_field.split(',') if alt]
        if len(alts) != 1:
            return 'multiallelic'

        alt = alts[0]
        if len(ref) == 1 and len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT':
            return 'snp'

        return 'indel'

    def _is_transition(self, ref: str, alt_field: str) -> bool:
        alt = alt_field.split(',')[0]
        pair = {ref.upper(), alt.upper()}
        return pair in ({'A', 'G'}, {'C', 'T'})

    def _normalize_vcf_genotype(self, raw_gt: str, ref: str, alt_field: str) -> Optional[str]:
        if not raw_gt or raw_gt in {'.', './.', '.|.'}:
            self.missing_genotypes += 1
            return None

        separator = '|' if '|' in raw_gt else '/'
        if separator == '|':
            self.phased_genotypes += 1
        else:
            self.unphased_genotypes += 1

        if separator not in raw_gt:
            if len(raw_gt) == 1 and raw_gt.upper() in {'A', 'C', 'G', 'T'}:
                self.haploid_genotypes += 1
                return raw_gt.upper() * 2
            return None

        alleles = raw_gt.split(separator)
        if len(alleles) != 2 or any(allele == '.' for allele in alleles):
            self.complex_genotypes += 1
            return None

        alt = alt_field.split(',')[0]
        allele_map = {'0': ref.upper(), '1': alt.upper()}

        def translate(token: str) -> Optional[str]:
            token = token.strip()
            if token in allele_map:
                return allele_map[token]
            if token.upper() in {'A', 'C', 'G', 'T'}:
                return token.upper()
            return None

        translated = [translate(token) for token in alleles]
        if any(allele is None for allele in translated):
            return None

        return ''.join(translated)  # type: ignore[arg-type]

    def _record_genotype(self, canonical_genotype: str, raw_gt: str, ref: str) -> None:
        self.genotype_counts[canonical_genotype] += 1

        if len(canonical_genotype) != 2:
            self.complex_genotypes += 1
            return

        separator = '|' if '|' in raw_gt else '/'
        raw_tokens = raw_gt.split(separator)
        if len(raw_tokens) == 2 and raw_tokens[0] == raw_tokens[1]:
            token = raw_tokens[0].strip()
            if token == '0':
                self.homozygous_reference += 1
            elif token == '1':
                self.homozygous_alternate += 1
            elif canonical_genotype == ref.upper() * 2:
                self.homozygous_reference += 1
            else:
                self.homozygous_alternate += 1
        else:
            self.heterozygous += 1

    def _extract_info_keys(self, info: str) -> List[str]:
        if not info or info == '.':
            return []

        keys = []
        for item in info.split(';'):
            item = item.strip()
            if not item:
                continue
            keys.append(item.split('=', 1)[0])
        return keys


def _format_file_size(size_bytes: Optional[int]) -> str:
    if size_bytes is None:
        return 'Unknown'

    size = float(size_bytes)
    units = ['B', 'KB', 'MB', 'GB', 'TB']
    unit_index = 0
    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024.0
        unit_index += 1

    if unit_index == 0:
        return f'{int(size)} {units[unit_index]}'
    return f'{size:.2f} {units[unit_index]}'
