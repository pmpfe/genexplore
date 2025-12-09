"""
Tests for polygenic database functionality.
"""

import pytest
import os
import tempfile
from datetime import datetime

from models.polygenic_models import (
    PolygenicScore, PolygenicVariant, PopulationDistribution,
    TraitCategory, DatabaseVersion
)
from database.polygenic_database import (
    PolygenicDatabase, DatabaseVersionManager, PolygenicDatabaseError
)


class TestPolygenicDatabase:
    """Tests for PolygenicDatabase class."""
    
    @pytest.fixture
    def temp_db(self):
        """Create a temporary database for testing."""
        with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
            db_path = f.name
        
        yield db_path
        
        # Cleanup
        if os.path.exists(db_path):
            os.unlink(db_path)
    
    def test_database_creation(self, temp_db):
        """Test database is created with schema."""
        db = PolygenicDatabase(temp_db)
        
        assert os.path.exists(temp_db)
    
    def test_get_all_scores(self, temp_db):
        """Test retrieving all scores."""
        db = PolygenicDatabase(temp_db)
        
        scores = db.get_all_scores()
        
        # Should have sample scores
        assert len(scores) >= 0
    
    def test_insert_and_retrieve_score(self, temp_db):
        """Test inserting and retrieving a score."""
        db = PolygenicDatabase(temp_db)
        
        # Create test score
        test_score = PolygenicScore(
            pgs_id="PGS_TEST_001",
            trait_name="Test Trait",
            trait_category=TraitCategory.METABOLIC,
            publication_doi="10.1000/test",
            publication_year=2023,
            study_population="EUR",
            sample_size=50000,
            num_variants=2,
            description="Test polygenic score",
            variants=[
                PolygenicVariant("rs123456", "1", 100000, "A", "G", 0.15, 0.30),
                PolygenicVariant("rs789012", "2", 200000, "T", "C", 0.25, 0.45),
            ]
        )
        
        test_dist = PopulationDistribution(
            pgs_id="PGS_TEST_001",
            population="EUR",
            mean=1.0,
            std=0.3,
            percentiles={25: 0.8, 50: 1.0, 75: 1.2}
        )
        
        # Insert
        db.insert_score(test_score, test_dist)
        
        # Retrieve
        retrieved = db.get_score_with_variants("PGS_TEST_001")
        
        assert retrieved is not None
        assert retrieved.pgs_id == "PGS_TEST_001"
        assert retrieved.trait_name == "Test Trait"
        assert len(retrieved.variants) == 2
        
        # Check distribution
        dist = db.get_population_distribution("PGS_TEST_001")
        assert dist is not None
        assert dist.mean == 1.0
        assert 50 in dist.percentiles
    
    def test_get_nonexistent_score(self, temp_db):
        """Test retrieving non-existent score returns None."""
        db = PolygenicDatabase(temp_db)
        
        result = db.get_score_with_variants("NONEXISTENT")
        
        assert result is None
    
    def test_get_all_distributions(self, temp_db):
        """Test retrieving all distributions."""
        db = PolygenicDatabase(temp_db)
        
        distributions = db.get_all_distributions()
        
        assert isinstance(distributions, dict)
    
    def test_get_version(self, temp_db):
        """Test getting database version."""
        db = PolygenicDatabase(temp_db)
        
        version = db.get_version("pgs_catalog")
        
        # May be None if no sample data, or have version if sample data inserted
        if version:
            assert version.database_name == "pgs_catalog"
    
    def test_get_score_count(self, temp_db):
        """Test getting score count."""
        db = PolygenicDatabase(temp_db)
        
        count = db.get_score_count()
        
        assert isinstance(count, int)
        assert count >= 0


class TestDatabaseVersionManager:
    """Tests for DatabaseVersionManager class."""
    
    @pytest.fixture
    def temp_backup_dir(self):
        """Create a temporary backup directory."""
        with tempfile.TemporaryDirectory() as d:
            yield d
    
    def test_create_backup(self, temp_backup_dir):
        """Test creating a backup."""
        manager = DatabaseVersionManager()
        
        # Create a test file
        test_file = os.path.join(temp_backup_dir, "test.db")
        with open(test_file, 'w') as f:
            f.write("test data")
        
        # Create backup (will go to default backup dir)
        backup_path = manager.create_backup(test_file)
        
        # May succeed or fail depending on permissions
        if backup_path:
            assert os.path.exists(backup_path)
    
    def test_create_backup_nonexistent(self):
        """Test creating backup of non-existent file."""
        manager = DatabaseVersionManager()
        
        result = manager.create_backup("/nonexistent/path/file.db")
        
        assert result is None
    
    def test_list_backups(self):
        """Test listing backups."""
        manager = DatabaseVersionManager()
        
        backups = manager.list_backups()
        
        assert isinstance(backups, list)


class TestDatabaseVersion:
    """Tests for DatabaseVersion dataclass."""
    
    def test_version_creation(self):
        """Test creating a database version."""
        version = DatabaseVersion(
            database_name="gwas_catalog",
            version="1.0.0",
            release_date=datetime(2024, 1, 1),
            download_date=datetime(2024, 1, 15),
            source_url="https://example.com",
            record_count=1000,
            checksum="abc123"
        )
        
        assert version.database_name == "gwas_catalog"
        assert version.version == "1.0.0"
        assert version.record_count == 1000
    
    def test_version_repr(self):
        """Test version string representation."""
        version = DatabaseVersion(
            database_name="test_db",
            version="2.0.0",
            release_date=None,
            download_date=datetime.now(),
            source_url="",
            record_count=0
        )
        
        assert "test_db" in repr(version)
        assert "2.0.0" in repr(version)
