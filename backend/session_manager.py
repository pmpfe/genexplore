"""
Session Manager for saving and loading user analysis results.

Provides functionality to serialize and deserialize complete analysis sessions
including genotype data, monogenic matches, and polygenic scores.
"""

import gzip
import json
import os
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import asdict

from models.data_models import SNPRecord, GWASMatch
from models.polygenic_models import PolygenicResult, RiskCategory, TraitCategory
from utils.logging_config import get_logger

logger = get_logger(__name__)

# File format version for compatibility checking
SESSION_FORMAT_VERSION = "1.0"


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
        metadata: Optional[Dict[str, Any]] = None
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
            session_data = {
                "format_version": SESSION_FORMAT_VERSION,
                "created_at": datetime.now().isoformat(),
                "metadata": metadata or {},
                "summary": {
                    "snp_count": len(snp_records),
                    "gwas_match_count": len(gwas_matches),
                    "polygenic_score_count": len(polygenic_results)
                },
                "snp_records": [
                    {
                        "rsid": snp.rsid,
                        "chromosome": snp.chromosome,
                        "position": snp.position,
                        "genotype": snp.genotype
                    }
                    for snp in snp_records
                ],
                "gwas_matches": [
                    {
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
                    for m in gwas_matches
                ],
                "polygenic_results": [
                    {
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
                    for r in polygenic_results
                ]
            }
            
            # Serialize to JSON and compress
            json_data = json.dumps(session_data, ensure_ascii=False)
            compressed_data = gzip.compress(json_data.encode('utf-8'), compresslevel=9)
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(filepath) if os.path.dirname(filepath) else '.', exist_ok=True)
            
            with open(filepath, 'wb') as f:
                f.write(compressed_data)
            
            # Log success
            original_size = len(json_data)
            compressed_size = len(compressed_data)
            compression_ratio = (1 - compressed_size / original_size) * 100 if original_size > 0 else 0
            
            logger.info(
                f"Session saved: {filepath} "
                f"({compressed_size:,} bytes, {compression_ratio:.1f}% compression)"
            )
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to save session: {e}")
            raise
    
    @staticmethod
    def load_session(
        filepath: str
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
            with open(filepath, 'rb') as f:
                compressed_data = f.read()
            
            # Decompress and parse JSON
            json_data = gzip.decompress(compressed_data).decode('utf-8')
            session_data = json.loads(json_data)
            
            # Check format version
            file_version = session_data.get("format_version", "unknown")
            if file_version != SESSION_FORMAT_VERSION:
                logger.warning(
                    f"Session file version mismatch: {file_version} vs {SESSION_FORMAT_VERSION}"
                )
            
            # Reconstruct SNP records
            snp_records = []
            for snp_data in session_data.get("snp_records", []):
                try:
                    snp = SNPRecord(
                        rsid=snp_data["rsid"],
                        chromosome=snp_data["chromosome"],
                        position=snp_data["position"],
                        genotype=snp_data["genotype"]
                    )
                    snp_records.append(snp)
                except (ValueError, KeyError) as e:
                    logger.warning(f"Skipping invalid SNP record: {e}")
            
            # Reconstruct GWAS matches
            gwas_matches = []
            for match_data in session_data.get("gwas_matches", []):
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
            
            # Reconstruct polygenic results
            polygenic_results = []
            for result_data in session_data.get("polygenic_results", []):
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
            
            # Build metadata
            metadata = session_data.get("metadata", {})
            metadata["file_version"] = file_version
            metadata["created_at"] = session_data.get("created_at")
            metadata["original_summary"] = session_data.get("summary", {})
            
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
            with open(filepath, 'rb') as f:
                compressed_data = f.read()
            
            # Decompress and parse JSON
            json_data = gzip.decompress(compressed_data).decode('utf-8')
            session_data = json.loads(json_data)
            
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
