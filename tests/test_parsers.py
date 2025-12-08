"""
Unit tests for the 23andMe file parser.
"""

import pytest
import tempfile
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from backend.parsers import Parser23andMe, ParseError, parse_23andme_file
from models.data_models import SNPRecord


class TestParser23andMe:
    """Tests for the Parser23andMe class."""
    
    def _create_temp_file(self, content: str) -> str:
        """Create a temporary file with the given content."""
        fd, path = tempfile.mkstemp(suffix='.txt')
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        return path
    
    def test_parse_valid_23andme_file(self):
        """Test parsing a valid 23andMe file."""
        content = """# Comment line
# Another comment
rs3131972	1	694713	GG
rs12124819	1	713790	AG
rs11240777	1	856331	GG
rs6681049	1	909917	TT
rs4970383	1	1019440	AG
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            assert len(records) == 5
            assert records[0].rsid == 'rs3131972'
            assert records[0].chromosome == '1'
            assert records[0].position == 694713
            assert records[0].genotype == 'GG'
            
            stats = parser.get_parse_stats()
            assert stats['valid_snps'] == 5
            assert stats['skipped_undetermined'] == 0
        finally:
            os.unlink(path)
    
    def test_parse_invalid_genotype_skipped(self):
        """Test that lines with invalid genotypes are skipped."""
        content = """rs3131972	1	694713	GG
rs7899632	2	1234567	--
rs12124819	1	713790	AG
rs1234567	3	5000000	ZZ
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            assert len(records) == 2
            assert records[0].rsid == 'rs3131972'
            assert records[1].rsid == 'rs12124819'
            
            stats = parser.get_parse_stats()
            assert stats['skipped_undetermined'] == 1
            assert stats['warnings_count'] == 1
        finally:
            os.unlink(path)
    
    def test_parse_missing_fields_skipped(self):
        """Test that lines with missing fields are skipped."""
        content = """rs3131972	1	694713	GG
rs12124819	1
rs11240777	1	856331	GG
invalid_line
rs6681049	1	909917	TT
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            assert len(records) == 3
            
            stats = parser.get_parse_stats()
            assert stats['warnings_count'] >= 2
        finally:
            os.unlink(path)
    
    def test_parse_invalid_rsid_skipped(self):
        """Test that lines with invalid RSIDs are skipped."""
        content = """rs3131972	1	694713	GG
invalid_id	1	713790	AG
12345	1	856331	GG
rs6681049	1	909917	TT
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            assert len(records) == 2
            rsids = [r.rsid for r in records]
            assert 'rs3131972' in rsids
            assert 'rs6681049' in rsids
        finally:
            os.unlink(path)
    
    def test_parse_counts_correct(self):
        """Test that parse statistics are accurate."""
        content = """# Header comment
rs3131972	1	694713	GG
rs12124819	1	713790	--
rs11240777	1	856331	AG
bad_line
rs6681049	1	909917	CC

# Another comment
rs4970383	1	1019440	TT
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            stats = parser.get_parse_stats()
            
            assert stats['valid_snps'] == 4
            assert stats['skipped_undetermined'] == 1
            assert len(records) == 4
        finally:
            os.unlink(path)
    
    def test_parse_file_not_found(self):
        """Test that ParseError is raised for missing files."""
        parser = Parser23andMe()
        
        with pytest.raises(ParseError):
            parser.parse_file('/nonexistent/path/file.txt')
    
    def test_parse_empty_file_error(self):
        """Test that ParseError is raised for empty files."""
        content = """# Only comments
# No data
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            with pytest.raises(ParseError):
                parser.parse_file(path)
        finally:
            os.unlink(path)
    
    def test_parse_chromosome_validation(self):
        """Test chromosome validation."""
        content = """rs3131972	1	694713	GG
rs12124819	X	713790	AG
rs11240777	Y	856331	GG
rs6681049	MT	909917	TT
rs4970383	99	1019440	AG
"""
        path = self._create_temp_file(content)
        try:
            parser = Parser23andMe()
            records = parser.parse_file(path)
            
            assert len(records) == 4
            chromosomes = [r.chromosome for r in records]
            assert '99' not in chromosomes
            assert 'X' in chromosomes
            assert 'Y' in chromosomes
            assert 'MT' in chromosomes
        finally:
            os.unlink(path)
    
    def test_convenience_function(self):
        """Test the parse_23andme_file convenience function."""
        content = """rs3131972	1	694713	GG
rs12124819	1	713790	AG
"""
        path = self._create_temp_file(content)
        try:
            records, stats = parse_23andme_file(path)
            
            assert len(records) == 2
            assert stats['valid_snps'] == 2
        finally:
            os.unlink(path)


class TestSNPRecord:
    """Tests for the SNPRecord dataclass."""
    
    def test_valid_snp_record(self):
        """Test creating a valid SNP record."""
        record = SNPRecord(
            rsid='rs3131972',
            chromosome='1',
            position=694713,
            genotype='GG'
        )
        
        assert record.rsid == 'rs3131972'
        assert record.chromosome == '1'
        assert record.position == 694713
        assert record.genotype == 'GG'
    
    def test_invalid_rsid_raises_error(self):
        """Test that invalid RSID raises ValueError."""
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='invalid',
                chromosome='1',
                position=694713,
                genotype='GG'
            )
    
    def test_invalid_chromosome_raises_error(self):
        """Test that invalid chromosome raises ValueError."""
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='rs3131972',
                chromosome='99',
                position=694713,
                genotype='GG'
            )
    
    def test_invalid_position_raises_error(self):
        """Test that invalid position raises ValueError."""
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='rs3131972',
                chromosome='1',
                position=0,
                genotype='GG'
            )
        
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='rs3131972',
                chromosome='1',
                position=-100,
                genotype='GG'
            )
    
    def test_invalid_genotype_raises_error(self):
        """Test that invalid genotype raises ValueError."""
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='rs3131972',
                chromosome='1',
                position=694713,
                genotype='ZZ'
            )
        
        with pytest.raises(ValueError):
            SNPRecord(
                rsid='rs3131972',
                chromosome='1',
                position=694713,
                genotype='A'
            )
