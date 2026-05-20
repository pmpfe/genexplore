"""
Session Manager for saving and loading user analysis results.

Provides functionality to serialize and deserialize complete analysis sessions
including genotype data, monogenic matches, and polygenic scores.
"""

import gzip
import json
import os
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple, Callable
from dataclasses import asdict

from models.data_models import SNPRecord, GWASMatch
from models.polygenic_models import PolygenicResult, RiskCategory, TraitCategory
from utils.logging_config import get_logger

logger = get_logger(__name__)

# File format version for compatibility checking
SESSION_FORMAT_VERSION = "1.0"


def _emit_progress(
    progress_callback: Optional[Callable[[int, str], None]],
    percent: int,
    message: str,
) -> None:
    """Emit a progress update if a callback is available."""
    if progress_callback is None:
        return

    try:
        progress_callback(max(0, min(100, int(percent))), message)
    except Exception:
        # Progress reporting must never break save/load operations.
        pass


def _stream_collection_progress(
    progress_callback: Optional[Callable[[int, str], None]],
    total_items: int,
    start_percent: int,
    end_percent: int,
    label: str,
    index: int,
) -> None:
    """Throttle progress updates while streaming large collections."""
    if total_items <= 0:
        _emit_progress(progress_callback, end_percent, f"{label} complete")
        return

    checkpoint = max(1, total_items // 50)
    if index == 1 or index == total_items or index % checkpoint == 0:
        span = max(1, end_percent - start_percent)
        percent = start_percent + int((index / total_items) * span)
        _emit_progress(progress_callback, percent, f"{label}: {index:,}/{total_items:,}")


class SessionManager:
    """
    Manages saving and loading of analysis sessions.
    
    Sessions include:
    - SNP records from parsed 23andMe file
    - Monogenic (GWAS) matches
    - Polygenic score results
    - Metadata (timestamps, versions, etc.)
    """
    
    @staticmethod
    def save_session(
        filepath: str,
        snp_records: List[SNPRecord],
        gwas_matches: List[GWASMatch],
        polygenic_results: List[PolygenicResult],
        metadata: Optional[Dict[str, Any]] = None,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> bool:
        """
        Save complete analysis session to compressed file.
        
        Args:
            filepath: Path for the output file (.gxs extension recommended)
            snp_records: List of SNP records from parsed file
            gwas_matches: List of GWAS matches (monogenic results)
            polygenic_results: List of polygenic score results
            metadata: Optional additional metadata to store
            
        Returns:
            bool: True if save was successful
        """
        try:
            _emit_progress(progress_callback, 0, "Preparing session data...")

            # Ensure directory exists
            os.makedirs(os.path.dirname(filepath) if os.path.dirname(filepath) else '.', exist_ok=True)

            summary = {
                "snp_count": len(snp_records),
                "gwas_match_count": len(gwas_matches),
                "polygenic_score_count": len(polygenic_results)
            }

            # Stream JSON output directly into a gzip file to avoid building
            # large intermediate strings in memory.
            with gzip.open(filepath, 'wt', encoding='utf-8', compresslevel=6) as f:
                _emit_progress(progress_callback, 5, "Writing session header...")
                f.write('{')
                f.write('"format_version":')
                json.dump(SESSION_FORMAT_VERSION, f, ensure_ascii=False)
                f.write(',"created_at":')
                json.dump(datetime.now().isoformat(), f, ensure_ascii=False)
                f.write(',"metadata":')
                json.dump(metadata or {}, f, ensure_ascii=False)
                f.write(',"summary":')
                json.dump(summary, f, ensure_ascii=False)

                # SNP records array (stream each record)
                f.write(',"snp_records":[')
                for i, snp in enumerate(snp_records):
                    if i:
                        f.write(',')
                    snp_obj = {
                        "rsid": snp.rsid,
                        "chromosome": snp.chromosome,
                        "position": snp.position,
                        "genotype": snp.genotype,
                        "source_format": snp.source_format,
                        "source_metadata": snp.source_metadata,
                        "is_valid": snp.is_valid,
                        "validation_notes": snp.validation_notes,
                        "variant_key": snp.variant_key,
                    }
                    json.dump(snp_obj, f, ensure_ascii=False)
                    _stream_collection_progress(
                        progress_callback,
                        len(snp_records),
                        10,
                        55,
                        "Writing SNP records",
                        i + 1,
                    )
                f.write(']')

                if not snp_records:
                    _emit_progress(progress_callback, 55, "Writing SNP records complete")

                # GWAS matches array
                f.write(',"gwas_matches":[')
                for i, m in enumerate(gwas_matches):
                    if i:
                        f.write(',')
                    match_obj = {
                        "rsid": m.rsid,
                        "chromosome": m.chromosome,
                        "position": m.position,
                        "user_genotype": m.user_genotype,
                        "gene": m.gene,
                        "trait": m.trait,
                        "risk_allele": m.risk_allele,
                        "p_value": m.p_value,
                        "odds_ratio": m.odds_ratio,
                        "sample_size": m.sample_size,
                        "category": m.category,
                        "allele_frequency": m.allele_frequency,
                        "impact_score": m.impact_score
                    }
                    json.dump(match_obj, f, ensure_ascii=False)
                    _stream_collection_progress(
                        progress_callback,
                        len(gwas_matches),
                        55,
                        80,
                        "Writing GWAS matches",
                        i + 1,
                    )
                f.write(']')

                if not gwas_matches:
                    _emit_progress(progress_callback, 80, "Writing GWAS matches complete")

                # Polygenic results array
                f.write(',"polygenic_results":[')
                for i, r in enumerate(polygenic_results):
                    if i:
                        f.write(',')
                    result_obj = {
                        "pgs_id": r.pgs_id,
                        "trait_name": r.trait_name,
                        "trait_category": r.trait_category.value,
                        "raw_score": r.raw_score,
                        "normalized_score": r.normalized_score,
                        "percentile": r.percentile,
                        "risk_category": r.risk_category.value,
                        "variants_found": r.variants_found,
                        "variants_total": r.variants_total,
                        "coverage_percent": r.coverage_percent,
                        "population_reference": r.population_reference,
                        "computation_time_ms": r.computation_time_ms
                    }
                    json.dump(result_obj, f, ensure_ascii=False)
                    _stream_collection_progress(
                        progress_callback,
                        len(polygenic_results),
                        80,
                        95,
                        "Writing polygenic scores",
                        i + 1,
                    )
                f.write(']')

                if not polygenic_results:
                    _emit_progress(progress_callback, 95, "Writing polygenic scores complete")

                f.write('}')

            _emit_progress(progress_callback, 100, "Session file written")

            compressed_size = os.path.getsize(filepath)
            logger.info(
                f"Session saved: {filepath} ({compressed_size:,} bytes)"
            )

            return True

        except Exception as e:
            logger.error(f"Failed to save session: {e}")
            raise
    
    @staticmethod
    def load_session(
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> Tuple[List[SNPRecord], List[GWASMatch], List[PolygenicResult], Dict[str, Any]]:
        """
        Load complete analysis session from compressed file.
        
        Args:
            filepath: Path to the session file
            
        Returns:
            Tuple containing:
            - List of SNP records
            - List of GWAS matches
            - List of polygenic results
            - Metadata dictionary
            
        Raises:
            ValueError: If file format is invalid or incompatible
            FileNotFoundError: If file doesn't exist
        """
        try:
            _emit_progress(progress_callback, 0, "Opening session file...")

            # Stream-load the JSON from the gzip file to avoid building large
            # intermediate byte strings in memory.
            with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                _emit_progress(progress_callback, 10, "Reading compressed session data...")
                session_data = json.load(f)

            _emit_progress(progress_callback, 20, "Session data read, rebuilding objects...")
            
            # Check format version
            file_version = session_data.get("format_version", "unknown")
            if file_version != SESSION_FORMAT_VERSION:
                logger.warning(
                    f"Session file version mismatch: {file_version} vs {SESSION_FORMAT_VERSION}"
                )
            
            # Reconstruct SNP records
            snp_records = []
            snp_source = session_data.get("snp_records", [])
            for index, snp_data in enumerate(snp_source, 1):
                try:
                    snp = SNPRecord(
                        rsid=snp_data.get("rsid"),
                        chromosome=snp_data["chromosome"],
                        position=snp_data["position"],
                        genotype=snp_data["genotype"],
                        source_format=snp_data.get("source_format", "23andMe raw"),
                        source_metadata=snp_data.get("source_metadata", {}),
                        is_valid=snp_data.get("is_valid", True),
                        validation_notes=snp_data.get("validation_notes", []),
                    )
                    snp_records.append(snp)
                except (ValueError, KeyError) as e:
                    logger.warning(f"Skipping invalid SNP record: {e}")
                
                _stream_collection_progress(
                    progress_callback,
                    len(snp_source),
                    20,
                    50,
                    "Rebuilding SNP records",
                    index,
                )

            if not snp_source:
                _emit_progress(progress_callback, 50, "Rebuilding SNP records complete")
            
            # Reconstruct GWAS matches
            gwas_matches = []
            gwas_source = session_data.get("gwas_matches", [])
            for index, match_data in enumerate(gwas_source, 1):
                try:
                    match = GWASMatch(
                        rsid=match_data["rsid"],
                        chromosome=match_data["chromosome"],
                        position=match_data["position"],
                        user_genotype=match_data["user_genotype"],
                        gene=match_data.get("gene"),
                        trait=match_data["trait"],
                        risk_allele=match_data["risk_allele"],
                        p_value=match_data["p_value"],
                        odds_ratio=match_data.get("odds_ratio"),
                        sample_size=match_data.get("sample_size"),
                        category=match_data["category"],
                        allele_frequency=match_data["allele_frequency"],
                        impact_score=match_data["impact_score"]
                    )
                    gwas_matches.append(match)
                except (ValueError, KeyError) as e:
                    logger.warning(f"Skipping invalid GWAS match: {e}")
                _stream_collection_progress(
                    progress_callback,
                    len(gwas_source),
                    50,
                    75,
                    "Rebuilding GWAS matches",
                    index,
                )

            if not gwas_source:
                _emit_progress(progress_callback, 75, "Rebuilding GWAS matches complete")
            
            # Reconstruct polygenic results
            polygenic_results = []
            polygenic_source = session_data.get("polygenic_results", [])
            for index, result_data in enumerate(polygenic_source, 1):
                try:
                    # Map trait category string to enum
                    trait_cat = TraitCategory.OTHER
                    for tc in TraitCategory:
                        if tc.value == result_data.get("trait_category"):
                            trait_cat = tc
                            break
                    
                    # Map risk category string to enum
                    risk_cat = RiskCategory.INTERMEDIATE
                    for rc in RiskCategory:
                        if rc.value == result_data.get("risk_category"):
                            risk_cat = rc
                            break
                    
                    result = PolygenicResult(
                        pgs_id=result_data["pgs_id"],
                        trait_name=result_data["trait_name"],
                        trait_category=trait_cat,
                        raw_score=result_data["raw_score"],
                        normalized_score=result_data["normalized_score"],
                        percentile=result_data["percentile"],
                        risk_category=risk_cat,
                        variants_found=result_data["variants_found"],
                        variants_total=result_data["variants_total"],
                        coverage_percent=result_data["coverage_percent"],
                        population_reference=result_data["population_reference"],
                        computation_time_ms=result_data.get("computation_time_ms")
                    )
                    polygenic_results.append(result)
                except (ValueError, KeyError) as e:
                    logger.warning(f"Skipping invalid polygenic result: {e}")
                _stream_collection_progress(
                    progress_callback,
                    len(polygenic_source),
                    75,
                    95,
                    "Rebuilding polygenic scores",
                    index,
                )

            if not polygenic_source:
                _emit_progress(progress_callback, 95, "Rebuilding polygenic scores complete")
            
            # Build metadata
            metadata = session_data.get("metadata", {})
            metadata["file_version"] = file_version
            metadata["created_at"] = session_data.get("created_at")
            metadata["original_summary"] = session_data.get("summary", {})

            _emit_progress(progress_callback, 100, "Session loaded")
            
            logger.info(
                f"Session loaded: {filepath} "
                f"({len(snp_records):,} SNPs, {len(gwas_matches):,} matches, "
                f"{len(polygenic_results)} polygenic scores)"
            )
            
            return snp_records, gwas_matches, polygenic_results, metadata
            
        except gzip.BadGzipFile:
            raise ValueError(f"Invalid session file format (not gzip): {filepath}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid session file format (invalid JSON): {e}")
        except Exception as e:
            logger.error(f"Failed to load session: {e}")
            raise
    
    @staticmethod
    def get_session_info(filepath: str) -> Dict[str, Any]:
        """
        Get summary information about a session file without fully loading it.
        
        Args:
            filepath: Path to the session file
            
        Returns:
            Dictionary with summary information
        """
        try:
            with gzip.open(filepath, 'rt', encoding='utf-8') as f:
                session_data = json.load(f)

            return {
                "filepath": filepath,
                "file_size_bytes": os.path.getsize(filepath),
                "format_version": session_data.get("format_version", "unknown"),
                "created_at": session_data.get("created_at"),
                "summary": session_data.get("summary", {}),
                "metadata": session_data.get("metadata", {})
            }

        except Exception as e:
            return {
                "filepath": filepath,
                "error": str(e)
            }
