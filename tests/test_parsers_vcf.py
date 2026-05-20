"""Tests for the VCF parser."""

import gzip
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.parsers_vcf import VCFStatsParser


class TestVCFStatsParser:
    def _write_temp_file(self, content: str, suffix: str = '.vcf') -> str:
        fd, path = tempfile.mkstemp(suffix=suffix)
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(content)
        return path

    def _write_temp_gz_file(self, content: str) -> str:
        fd, path = tempfile.mkstemp(suffix='.vcf.gz')
        os.close(fd)
        with gzip.open(path, 'wt', encoding='utf-8') as handle:
            handle.write(content)
        return path

    def test_parse_vcf_records(self):
        content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT\t0/1
1\t200\t.\tC\tT\t99\tPASS\t.\tGT\t1/1
2\t300\trs3\tG\tGA\t10\tq10\tDP=20\tGT\t0/0
"""
        path = self._write_temp_file(content)
        try:
            parser = VCFStatsParser()
            records = parser.parse_file(path)

            assert len(records) == 2
            assert records[0].rsid == 'rs1'
            assert records[0].genotype == 'AG'
            assert records[1].rsid is None
            assert records[1].variant_key == '1:200:C>T'
            assert records[1].genotype == 'TT'

            stats = parser.inspect_file(path).stats
            assert stats['analysis_records'] == 2
            assert stats['skipped_non_snp'] == 1
            assert stats['variant_types']['snp'] == 2
            assert stats['variant_types']['indel'] == 1
        finally:
            os.unlink(path)

    def test_parse_vcf_gz_records(self):
        content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT\t0/1
"""
        path = self._write_temp_gz_file(content)
        try:
            parser = VCFStatsParser()
            records = parser.parse_file(path)

            assert len(records) == 1
            assert records[0].rsid == 'rs1'
            assert records[0].genotype == 'AG'
        finally:
            os.unlink(path)

    def test_vcf_progress_is_reported_as_percentage(self):
        header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""
        record = "1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT\t0/1\n"
        content = header + record * 10000
        path = self._write_temp_file(content)
        progress_values = []

        try:
            parser = VCFStatsParser()
            parser.inspect_file(path, progress_callback=lambda percent, message: progress_values.append(percent))

            assert progress_values
            assert max(progress_values) <= 100
            assert progress_values[-1] == 100
        finally:
            os.unlink(path)
