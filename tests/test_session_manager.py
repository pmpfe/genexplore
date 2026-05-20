"""Tests for session persistence and progress reporting."""

import os
import tempfile

import pytest

from backend.session_manager import SessionManager
from models.data_models import GWASMatch, SNPRecord
from models.polygenic_models import PolygenicResult, RiskCategory, TraitCategory


def _build_sample_session():
    snp_records = [
        SNPRecord(
            rsid="rs123456",
            chromosome="1",
            position=123456,
            genotype="AG",
        ),
        SNPRecord(
            rsid="rs789012",
            chromosome="2",
            position=234567,
            genotype="CC",
        ),
    ]

    gwas_matches = [
        GWASMatch(
            rsid="rs123456",
            chromosome="1",
            position=123456,
            user_genotype="AG",
            gene="GENE1",
            trait="Test trait",
            risk_allele="A",
            p_value=1e-8,
            odds_ratio=1.3,
            sample_size=10_000,
            category="Metabolic",
            allele_frequency=0.12,
            impact_score=7.5,
        )
    ]

    polygenic_results = [
        PolygenicResult(
            pgs_id="PGS000001",
            trait_name="Test polygenic trait",
            trait_category=TraitCategory.NEUROPSYCHIATRIC,
            raw_score=1.23,
            normalized_score=0.45,
            percentile=85.5,
            risk_category=RiskCategory.HIGH,
            variants_found=5,
            variants_total=10,
            coverage_percent=50.0,
            population_reference="EUR",
            variant_contributions=[("rs1", 0.12), ("rs2", -0.08)],
            computation_time_ms=12.3,
        )
    ]

    metadata = {
        "source_filename": "sample_23andme.txt",
        "source_format": "23andMe raw",
        "analysis_started_at": "2026-05-20T10:00:00",
    }

    return snp_records, gwas_matches, polygenic_results, metadata


def _is_monotonic(progress_updates):
    values = [value for value, _ in progress_updates]
    return all(first <= second for first, second in zip(values, values[1:]))


class TestSessionManager:
    @pytest.fixture
    def temp_session_path(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            yield os.path.join(temp_dir, "session.gxs")

    def test_save_session_reports_progress(self, temp_session_path):
        snp_records, gwas_matches, polygenic_results, metadata = _build_sample_session()
        progress_updates = []

        saved = SessionManager.save_session(
            temp_session_path,
            snp_records,
            gwas_matches,
            polygenic_results,
            metadata=metadata,
            progress_callback=lambda percent, message: progress_updates.append((percent, message)),
        )

        assert saved is True
        assert os.path.exists(temp_session_path)
        assert progress_updates
        assert progress_updates[0][0] == 0
        assert progress_updates[-1][0] == 100
        assert _is_monotonic(progress_updates)
        assert any("Writing SNP records" in message for _, message in progress_updates)
        assert any("Writing GWAS matches" in message for _, message in progress_updates)
        assert any("Writing polygenic scores" in message for _, message in progress_updates)

    def test_load_session_reports_progress_and_roundtrip(self, temp_session_path):
        snp_records, gwas_matches, polygenic_results, metadata = _build_sample_session()
        SessionManager.save_session(
            temp_session_path,
            snp_records,
            gwas_matches,
            polygenic_results,
            metadata=metadata,
        )

        progress_updates = []
        loaded = SessionManager.load_session(
            temp_session_path,
            progress_callback=lambda percent, message: progress_updates.append((percent, message)),
        )

        loaded_snps, loaded_matches, loaded_results, loaded_metadata = loaded

        assert len(loaded_snps) == len(snp_records)
        assert len(loaded_matches) == len(gwas_matches)
        assert len(loaded_results) == len(polygenic_results)
        assert loaded_metadata["original_summary"]["snp_count"] == len(snp_records)
        assert progress_updates[0][0] == 0
        assert progress_updates[-1][0] == 100
        assert _is_monotonic(progress_updates)
        assert any("Rebuilding SNP records" in message for _, message in progress_updates)
        assert any("Rebuilding GWAS matches" in message for _, message in progress_updates)
        assert any("Rebuilding polygenic scores" in message for _, message in progress_updates)
