"""
Tests for the genetic file inspection flow.
"""

import gzip
import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.file_inspector import inspect_genetic_file
from backend.parser_common import ParseError


class TestGeneticFileInspection:
    """Tests for inspecting supported genetic file formats."""

    def _write_temp_file(self, content: str, suffix: str) -> str:
        fd, path = tempfile.mkstemp(suffix=suffix)
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(content)
        return path

    def _write_temp_gz_file(self, content: str, suffix: str = '.vcf.gz') -> str:
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        with gzip.open(path, 'wt', encoding='utf-8') as handle:
            handle.write(content)
        return path

    def test_inspect_23andme_file(self):
        content = """# comment
rs3131972\t1\t694713\tGG
rs12124819\t1\t713790\tAG
"""
        path = self._write_temp_file(content, '.txt')
        try:
            result = inspect_genetic_file(path)

            assert result.analysis_supported is True
            assert result.detected_format == '23andMe raw'
            assert result.records is not None
            assert len(result.records) == 2
            assert result.stats['valid_snps'] == 2
            assert result.stats['genetic_records'] == 2
            assert result.stats['comment_lines'] == 1
            assert result.stats['chromosome_counts'] == {'1': 2}
            assert result.stats['genotype_counts']['GG'] == 1
            assert result.stats['homozygous_snps'] == 1
            assert result.stats['heterozygous_snps'] == 1
        finally:
            os.unlink(path)

    def test_inspect_vcf_file(self):
        content = """##fileformat=VCFv4.2
##source=TellmeGen
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT\t0/1
1\t200\t.\tC\tT\t99\tPASS\t.\tGT\t1/1
2\t300\trs3\tG\tGA\t10\tq10\tDP=20\tGT\t0/0
2\t400\trs4\tT\tC,G\t10\tPASS\tDP=10\tGT\t1/2
"""
        path = self._write_temp_file(content, '.vcf')
        try:
            result = inspect_genetic_file(path)
            stats = result.stats

            assert result.analysis_supported is True
            assert result.detected_format == 'VCF'
            assert len(result.records) == 2
            assert stats['genetic_records'] == 2
            assert stats['vcf_version'] == 'VCFv4.2'
            assert stats['metadata_lines'] == 2
            assert stats['header_lines'] == 1
            assert stats['data_lines'] == 4
            assert stats['sample_count'] == 1
            assert stats['samples'] == ['SAMPLE1']
            assert stats['variant_types']['snp'] == 2
            assert stats['variant_types']['indel'] == 1
            assert stats['variant_types']['multiallelic'] == 1
            assert stats['chromosome_counts'] == {'1': 2}
            assert stats['called_genotypes'] == 2
            assert stats['missing_genotypes'] == 0
            assert stats['homozygous_reference_genotypes'] == 0
            assert stats['homozygous_alternate_genotypes'] == 1
            assert stats['heterozygous_genotypes'] == 1
            assert stats['complex_genotypes'] == 0
            assert stats['skipped_non_snp'] == 2
            assert stats['transition_snps'] == 2
            assert stats['transversion_snps'] == 0
        finally:
            os.unlink(path)

    def test_inspect_vcf_gz_file(self):
        content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT\t0/1
"""
        path = self._write_temp_gz_file(content)
        try:
            result = inspect_genetic_file(path)

            assert result.analysis_supported is True
            assert result.is_compressed is True
            assert result.detected_format == 'VCF'
            assert result.stats['genetic_records'] == 1
            assert result.stats['vcf_version'] == 'VCFv4.2'
            assert result.stats['data_lines'] == 1
            assert result.stats['sample_count'] == 1
        finally:
            os.unlink(path)

    def test_missing_file_raises(self):
        with pytest.raises(ParseError):
            inspect_genetic_file('/missing/file.vcf')
