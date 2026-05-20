"""End-to-end report export tests based on the bundled sample 23andMe file."""

from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest

from backend.file_inspector import inspect_genetic_file
from backend.search_engine import SearchEngine
from config import DATABASE_PATH
from database import polygenic_database as polygenic_database_module
from models.polygenic_models import (
    DatabaseVersion,
    PopulationDistribution,
    PolygenicResult,
    PolygenicScore,
    RiskCategory,
    TraitCategory,
)


SAMPLE_FILE = Path(__file__).resolve().parents[1] / "sample_23andme.txt"


@pytest.fixture
def sample_inspection():
    return inspect_genetic_file(SAMPLE_FILE)


@pytest.fixture
def sample_matches(sample_inspection):
    engine = SearchEngine(DATABASE_PATH)
    return engine.match_user_snps(sample_inspection.records)


def test_sample_23andme_html_export_contains_report_sections(tmp_path, monkeypatch, sample_inspection, sample_matches):
    gwas_version = DatabaseVersion(
        database_name="gwas_catalog",
        version="9.9.9",
        release_date=datetime(2026, 1, 1, 12, 0, 0),
        download_date=datetime(2026, 1, 2, 12, 0, 0),
        source_url="https://example.invalid/gwas",
        record_count=123456,
        checksum="abc123",
    )
    pgs_version = DatabaseVersion(
        database_name="pgs_catalog",
        version="1.2.3",
        release_date=datetime(2026, 1, 5, 12, 0, 0),
        download_date=datetime(2026, 1, 6, 12, 0, 0),
        source_url="https://example.invalid/pgs",
        record_count=42,
        checksum="def456",
    )

    monkeypatch.setattr(polygenic_database_module.DatabaseVersionManager, "get_gwas_version", lambda self: gwas_version)
    monkeypatch.setattr(polygenic_database_module.DatabaseVersionManager, "get_pgs_version", lambda self: pgs_version)
    monkeypatch.setattr(
        polygenic_database_module,
        "get_gwas_database_stats",
        lambda: {"variants": 123456, "traits": 789, "genes": 456},
    )

    from frontend.main_window import MainWindow
    from PyQt6.QtWidgets import QApplication
    import sys

    # Avoid singleton issue if another test already created an app
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    main_window = MainWindow()

    # Stub the actually exported file path
    output_file = tmp_path / "report.html"

    # We mock the export dialog to return our temp path
    monkeypatch.setattr("PyQt6.QtWidgets.QFileDialog.getSaveFileName", lambda *args, **kwargs: (str(output_file), "HTML Files (*.html)"))

    # Trigger export
    # Usually this would be: main_window.export_report()
    # But that might open dialogs. Let's call the internal reporter directly if possible,
    # or ensure the mock works.
    
    # Fill some mock data into the main window if needed
    main_window.current_inspection = sample_inspection
    main_window.current_matches = sample_matches
    # Mock polygenic results
    main_window.polygenic_results = [
        PolygenicResult(
            score=PolygenicScore(
                pgs_id="PGS000001",
                trait_reported="Externalizing behaviors",
                method="LDPred2",
                genome_build="GRCh37",
                variant_count=100,
            ),
            calculated_score=1.23,
            percentile=85.5,
            population_distribution=PopulationDistribution(mean=0.0, std_dev=1.0, cohort="European"),
            risk_category=RiskCategory.HIGH,
            trait_category=TraitCategory.COGNITIVE,
        )
    ]

    main_window.export_report()

    assert output_file.exists()
    html_content = output_file.read_text()

    # Check for major sections
    assert "Genetic Exploration Report" in html_content
    assert "Sample Information" in html_content
    assert "Summary" in html_content
    assert "GWAS Matches" in html_content
    assert "Polygenic Risk Scores" in html_content
    assert "Methodology & Limitations" in html_content
    
    # Check for specific data points from our mocks
    assert "9.9.9" in html_content  # GWAS version
    assert "1.2.3" in html_content  # PGS version
    assert "123,456" in html_content # Variants count (formatted)
    assert "Externalizing behaviors" in html_content
    assert "85.5%" in html_content
