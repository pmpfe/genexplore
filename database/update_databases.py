#!/usr/bin/env python3
"""
Database Update Script for GenExplore
======================================

Downloads and updates GWAS and PGS Catalog databases.

Features:
- Resumable downloads (tracks progress in state file)
- Progress indicators with time estimates
- Separate GWAS and PGS updates
- Smart PGS filtering (excludes scores with >100K variants for performance)

Usage:
    python update_databases.py gwas      # Update GWAS only
    python update_databases.py pgs       # Update PGS only
    python update_databases.py all       # Update both
    python update_databases.py status    # Show current status
"""

import os
import sys
import json
import gzip
import sqlite3
import time
import argparse
import hashlib
from datetime import datetime, timedelta
from typing import Optional, Dict, List, Tuple, Set
from pathlib import Path
from io import BytesIO

import requests
from tqdm import tqdm

# Paths
SCRIPT_DIR = Path(__file__).parent
DB_DIR = SCRIPT_DIR
STATE_FILE = DB_DIR / "update_state.json"
GWAS_DB = DB_DIR / "gwas.db"
PGS_DB = DB_DIR / "pgs.db"

# API endpoints
GWAS_BULK_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv"
GWAS_ZIP_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
PGS_API_BASE = "https://www.pgscatalog.org/rest"

# Limits for PGS (scores with >100K variants are too large)
MAX_VARIANTS_PER_SCORE = 100_000
REQUEST_TIMEOUT = 60
REQUEST_DELAY = 0.5  # Be nice to the API


def load_state() -> Dict:
    """Load update state from file."""
    if STATE_FILE.exists():
        try:
            with open(STATE_FILE, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {
        "gwas": {"last_update": None, "status": "never"},
        "pgs": {
            "last_update": None, 
            "status": "never",
            "processed_scores": [],
            "failed_scores": [],
            "total_scores": 0,
            "total_variants": 0
        }
    }


def save_state(state: Dict) -> None:
    """Save update state to file."""
    with open(STATE_FILE, 'w') as f:
        json.dump(state, f, indent=2, default=str)


def format_time(seconds: float) -> str:
    """Format seconds as human-readable time."""
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        return f"{int(seconds // 60)}m {int(seconds % 60)}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"


def format_size(bytes_size: int) -> str:
    """Format bytes as human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_size < 1024:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024
    return f"{bytes_size:.1f} TB"


# =============================================================================
# GWAS DATABASE
# =============================================================================

def init_gwas_db() -> None:
    """Initialize GWAS database schema."""
    conn = sqlite3.connect(GWAS_DB)
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gwas_associations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            rsid TEXT,
            chromosome TEXT,
            position INTEGER,
            effect_allele TEXT,
            other_allele TEXT,
            p_value REAL,
            odds_ratio REAL,
            beta REAL,
            ci_lower REAL,
            ci_upper REAL,
            trait TEXT,
            trait_uri TEXT,
            study_accession TEXT,
            pubmed_id TEXT,
            first_author TEXT,
            publication_date TEXT,
            journal TEXT,
            sample_size INTEGER,
            ancestry TEXT,
            risk_frequency REAL
        )
    """)
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)
    
    # Create indexes for common queries
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON gwas_associations(rsid)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gwas_trait ON gwas_associations(trait)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gwas_chr_pos ON gwas_associations(chromosome, position)")
    
    conn.commit()
    conn.close()


def update_gwas(state: Dict) -> bool:
    """
    Download and update GWAS database.
    
    Uses the bulk download ZIP file from GWAS Catalog FTP.
    """
    import zipfile
    
    print("\n" + "="*60)
    print("📥 GWAS CATALOG UPDATE")
    print("="*60)
    
    init_gwas_db()
    
    print("\n📊 Downloading GWAS bulk data...")
    print(f"   URL: {GWAS_ZIP_URL}")
    
    try:
        # Stream download with progress
        response = requests.get(GWAS_ZIP_URL, stream=True, timeout=300)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        print(f"   Size: {format_size(total_size)}")
        
        # Download with progress bar
        data = BytesIO()
        with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                data.write(chunk)
                pbar.update(len(chunk))
        
        data.seek(0)
        
        # Extract ZIP
        print("\n📦 Extracting ZIP file...")
        with zipfile.ZipFile(data, 'r') as zf:
            # Find the TSV file
            tsv_files = [f for f in zf.namelist() if f.endswith('.tsv')]
            if not tsv_files:
                print("❌ No TSV file found in ZIP")
                return False
            
            tsv_file = tsv_files[0]
            print(f"   Found: {tsv_file}")
            
            with zf.open(tsv_file) as f:
                content = f.read().decode('utf-8')
        
        lines = content.strip().split('\n')
        print(f"\n✅ Extracted {len(lines):,} lines")
        
    except requests.RequestException as e:
        print(f"\n❌ Download failed: {e}")
        return False
    except zipfile.BadZipFile as e:
        print(f"\n❌ Invalid ZIP file: {e}")
        return False
    
    # Parse and insert data
    print("\n📝 Parsing and inserting into database...")
    
    conn = sqlite3.connect(GWAS_DB)
    cursor = conn.cursor()
    
    # Clear existing data
    cursor.execute("DELETE FROM gwas_associations")
    
    # Parse header
    if not lines:
        print("❌ No data received")
        return False
    
    header = lines[0].split('\t')
    header_map = {col.strip().upper(): i for i, col in enumerate(header)}
    
    # Column mappings (GWAS Catalog format)
    col_mappings = {
        'rsid': ['SNPS', 'SNP_ID_CURRENT', 'STRONGEST SNP-RISK ALLELE'],
        'chromosome': ['CHR_ID', 'CHROMOSOME'],
        'position': ['CHR_POS', 'POSITION'],
        'p_value': ['P-VALUE', 'PVALUE'],
        'odds_ratio': ['OR OR BETA', 'ODDS RATIO', 'OR'],
        'trait': ['DISEASE/TRAIT', 'MAPPED_TRAIT', 'TRAIT'],
        'trait_uri': ['MAPPED_TRAIT_URI'],
        'study_accession': ['STUDY ACCESSION', 'STUDY_ACCESSION'],
        'pubmed_id': ['PUBMEDID', 'PUBMED_ID'],
        'first_author': ['FIRST AUTHOR', 'FIRST_AUTHOR'],
        'sample_size': ['INITIAL SAMPLE SIZE', 'SAMPLE_SIZE'],
        'risk_frequency': ['RISK ALLELE FREQUENCY', 'RAF']
    }
    
    def get_col_idx(names: List[str]) -> Optional[int]:
        for name in names:
            if name.upper() in header_map:
                return header_map[name.upper()]
        return None
    
    # Find column indices
    col_indices = {field: get_col_idx(names) for field, names in col_mappings.items()}
    
    # Insert rows
    inserted = 0
    skipped = 0
    batch = []
    batch_size = 1000
    
    for line in tqdm(lines[1:], desc="Processing", unit="rows"):
        cols = line.split('\t')
        
        try:
            rsid = cols[col_indices['rsid']] if col_indices['rsid'] is not None and col_indices['rsid'] < len(cols) else None
            
            # Extract rsID from "rs123-A" format if needed
            if rsid and '-' in rsid:
                rsid = rsid.split('-')[0]
            
            # Skip if no rsID
            if not rsid or not rsid.startswith('rs'):
                skipped += 1
                continue
            
            # Parse numeric fields safely
            def safe_float(idx):
                if idx is None or idx >= len(cols):
                    return None
                try:
                    val = cols[idx].strip()
                    if not val or val == 'NR' or val == 'NA':
                        return None
                    return float(val)
                except (ValueError, IndexError):
                    return None
            
            def safe_int(idx):
                val = safe_float(idx)
                return int(val) if val is not None else None
            
            def safe_str(idx):
                if idx is None or idx >= len(cols):
                    return None
                val = cols[idx].strip()
                return val if val and val != 'NR' else None
            
            batch.append((
                rsid,
                safe_str(col_indices['chromosome']),
                safe_int(col_indices['position']),
                None,  # effect_allele
                None,  # other_allele
                safe_float(col_indices['p_value']),
                safe_float(col_indices['odds_ratio']),
                None,  # beta
                None,  # ci_lower
                None,  # ci_upper
                safe_str(col_indices['trait']),
                safe_str(col_indices['trait_uri']),
                safe_str(col_indices['study_accession']),
                safe_str(col_indices['pubmed_id']),
                safe_str(col_indices['first_author']),
                None,  # publication_date
                None,  # journal
                safe_int(col_indices['sample_size']),
                None,  # ancestry
                safe_float(col_indices['risk_frequency'])
            ))
            
            if len(batch) >= batch_size:
                cursor.executemany("""
                    INSERT INTO gwas_associations (
                        rsid, chromosome, position, effect_allele, other_allele,
                        p_value, odds_ratio, beta, ci_lower, ci_upper,
                        trait, trait_uri, study_accession, pubmed_id, first_author,
                        publication_date, journal, sample_size, ancestry, risk_frequency
                    ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                """, batch)
                inserted += len(batch)
                batch = []
                
        except Exception as e:
            skipped += 1
            continue
    
    # Insert remaining batch
    if batch:
        cursor.executemany("""
            INSERT INTO gwas_associations (
                rsid, chromosome, position, effect_allele, other_allele,
                p_value, odds_ratio, beta, ci_lower, ci_upper,
                trait, trait_uri, study_accession, pubmed_id, first_author,
                publication_date, journal, sample_size, ancestry, risk_frequency
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, batch)
        inserted += len(batch)
    
    # Update metadata
    cursor.execute(
        "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
        ("last_update", datetime.now().isoformat())
    )
    
    conn.commit()
    conn.close()
    
    # Update state
    state["gwas"]["last_update"] = datetime.now().isoformat()
    state["gwas"]["status"] = "complete"
    state["gwas"]["associations"] = inserted
    save_state(state)
    
    print(f"\n✅ GWAS update complete!")
    print(f"   Inserted: {inserted:,} associations")
    print(f"   Skipped: {skipped:,} rows")
    
    return True


# =============================================================================
# PGS CATALOG DATABASE
# =============================================================================

def init_pgs_db() -> None:
    """Initialize PGS database schema."""
    conn = sqlite3.connect(PGS_DB)
    cursor = conn.cursor()
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS polygenic_scores (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pgs_id TEXT UNIQUE NOT NULL,
            trait_name TEXT NOT NULL,
            trait_category TEXT,
            trait_efo_id TEXT,
            publication_doi TEXT,
            publication_year INTEGER,
            publication_title TEXT,
            sample_size INTEGER,
            num_variants INTEGER,
            num_variants_downloaded INTEGER DEFAULT 0,
            ancestry TEXT,
            genome_build TEXT,
            download_status TEXT DEFAULT 'pending',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pgs_variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pgs_id TEXT NOT NULL,
            rsid TEXT,
            chromosome TEXT,
            position INTEGER,
            effect_allele TEXT,
            other_allele TEXT,
            effect_weight REAL NOT NULL,
            allele_frequency REAL,
            FOREIGN KEY (pgs_id) REFERENCES polygenic_scores(pgs_id)
        )
    """)
    
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)
    
    # Indexes
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_pgs_variants_pgs_id ON pgs_variants(pgs_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_pgs_variants_rsid ON pgs_variants(rsid)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_pgs_trait ON polygenic_scores(trait_name)")
    
    conn.commit()
    conn.close()


def get_all_pgs_metadata() -> List[Dict]:
    """Fetch metadata for all PGS scores from API."""
    print("\n📊 Fetching PGS Catalog index...")
    
    all_scores = []
    next_url = f"{PGS_API_BASE}/score/all?limit=100"
    page = 0
    
    with tqdm(desc="Fetching score list", unit="scores") as pbar:
        while next_url:
            try:
                response = requests.get(next_url, timeout=REQUEST_TIMEOUT)
                response.raise_for_status()
                data = response.json()
                
                results = data.get('results', [])
                all_scores.extend(results)
                pbar.update(len(results))
                
                next_url = data.get('next')
                page += 1
                
                # Rate limiting
                time.sleep(REQUEST_DELAY)
                
            except requests.RequestException as e:
                print(f"\n⚠️ Error fetching page {page}: {e}")
                break
    
    return all_scores


def download_pgs_scoring_file(pgs_id: str) -> Optional[List[Dict]]:
    """
    Download scoring file for a PGS score.
    
    Returns list of variant dictionaries or None if failed.
    """
    # Scoring file URL format
    url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/{pgs_id}.txt.gz"
    
    try:
        response = requests.get(url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        
        # Decompress gzip
        content = gzip.decompress(response.content).decode('utf-8')
        lines = content.strip().split('\n')
        
        # Skip header comments
        data_lines = [l for l in lines if not l.startswith('#')]
        if not data_lines:
            return None
        
        # Parse header
        header = data_lines[0].split('\t')
        header_lower = [h.lower().strip() for h in header]
        
        # Map columns
        col_map = {}
        for i, col in enumerate(header_lower):
            if col in ['rsid', 'snp', 'snp_id']:
                col_map['rsid'] = i
            elif col in ['chr_name', 'chromosome', 'chr']:
                col_map['chromosome'] = i
            elif col in ['chr_position', 'position', 'pos', 'bp']:
                col_map['position'] = i
            elif col in ['effect_allele', 'a1', 'allele1', 'ea']:
                col_map['effect_allele'] = i
            elif col in ['other_allele', 'a2', 'allele2', 'oa', 'reference_allele']:
                col_map['other_allele'] = i
            elif col in ['effect_weight', 'weight', 'beta']:
                col_map['weight'] = i
            elif col in ['allelefrequency_effect', 'eaf', 'effect_allele_frequency']:
                col_map['frequency'] = i
        
        # Must have weight column
        if 'weight' not in col_map:
            return None
        
        variants = []
        for line in data_lines[1:]:
            cols = line.split('\t')
            
            try:
                weight = float(cols[col_map['weight']])
            except (ValueError, IndexError):
                continue
            
            variant = {
                'rsid': cols[col_map.get('rsid', -1)] if col_map.get('rsid', -1) < len(cols) else None,
                'chromosome': cols[col_map.get('chromosome', -1)] if col_map.get('chromosome', -1) < len(cols) else None,
                'position': None,
                'effect_allele': cols[col_map.get('effect_allele', -1)] if col_map.get('effect_allele', -1) < len(cols) else None,
                'other_allele': cols[col_map.get('other_allele', -1)] if col_map.get('other_allele', -1) < len(cols) else None,
                'weight': weight,
                'frequency': None
            }
            
            # Parse position
            if col_map.get('position') is not None and col_map['position'] < len(cols):
                try:
                    variant['position'] = int(cols[col_map['position']])
                except ValueError:
                    pass
            
            # Parse frequency
            if col_map.get('frequency') is not None and col_map['frequency'] < len(cols):
                try:
                    variant['frequency'] = float(cols[col_map['frequency']])
                except ValueError:
                    pass
            
            variants.append(variant)
        
        return variants
        
    except Exception as e:
        return None


def update_pgs(state: Dict, resume: bool = True) -> bool:
    """
    Download and update PGS database.
    
    Features:
    - Skips scores with >100K variants (too large)
    - Resumable (tracks processed scores)
    - Progress with time estimates
    """
    print("\n" + "="*60)
    print("📥 PGS CATALOG UPDATE")
    print("="*60)
    
    init_pgs_db()
    
    # Get all score metadata
    all_scores = get_all_pgs_metadata()
    if not all_scores:
        print("❌ Failed to fetch PGS catalog")
        return False
    
    print(f"\n📊 Found {len(all_scores):,} scores in PGS Catalog")
    
    # Filter out huge scores
    eligible_scores = [
        s for s in all_scores 
        if (s.get('variants_number') or 0) <= MAX_VARIANTS_PER_SCORE
    ]
    
    huge_scores = len(all_scores) - len(eligible_scores)
    total_variants = sum(s.get('variants_number', 0) for s in eligible_scores)
    
    print(f"   Eligible scores (≤{MAX_VARIANTS_PER_SCORE:,} variants): {len(eligible_scores):,}")
    print(f"   Skipped (too large): {huge_scores:,}")
    print(f"   Total variants to download: {total_variants:,}")
    
    # Load already processed scores
    processed: Set[str] = set()
    if resume:
        processed = set(state["pgs"].get("processed_scores", []))
        failed = set(state["pgs"].get("failed_scores", []))
        print(f"   Already processed: {len(processed):,}")
        print(f"   Previously failed: {len(failed):,}")
    
    # Filter to only unprocessed
    to_process = [s for s in eligible_scores if s.get('id') not in processed]
    
    if not to_process:
        print("\n✅ All eligible scores already downloaded!")
        return True
    
    # Calculate remaining work
    remaining_variants = sum(s.get('variants_number', 0) for s in to_process)
    print(f"\n📥 Scores to download: {len(to_process):,}")
    print(f"   Variants remaining: {remaining_variants:,}")
    
    # Estimate time (based on ~500 variants/second average)
    est_seconds = remaining_variants / 500
    print(f"   Estimated time: {format_time(est_seconds)}")
    
    # Connect to database
    conn = sqlite3.connect(PGS_DB)
    cursor = conn.cursor()
    
    # Process scores with weighted progress
    scores_done = 0
    scores_failed = 0
    variants_done = 0
    start_time = time.time()
    
    # Progress bar based on variants (weighted)
    pbar = tqdm(total=remaining_variants, desc="Downloading", unit="var", unit_scale=True)
    
    for score in to_process:
        pgs_id = score.get('id')
        num_variants = score.get('variants_number', 0)
        trait = score.get('trait_reported', 'Unknown')
        
        if not pgs_id:
            continue
        
        # Download scoring file
        variants = download_pgs_scoring_file(pgs_id)
        
        if variants is None:
            scores_failed += 1
            state["pgs"]["failed_scores"] = list(set(state["pgs"].get("failed_scores", [])) | {pgs_id})
            pbar.update(num_variants)  # Still advance progress
            continue
        
        # Insert score metadata
        try:
            # Safe extraction of nested fields
            trait_efo = score.get('trait_efo') or []
            trait_efo_id = trait_efo[0].get('id') if isinstance(trait_efo, list) and trait_efo else None
            
            publication = score.get('publication') or {}
            pub_doi = publication.get('doi')
            pub_date = publication.get('date_publication') or ''
            pub_year = int(pub_date[:4]) if pub_date and len(pub_date) >= 4 else None
            pub_title = publication.get('title')
            
            samples = score.get('samples_variants') or []
            sample_size = samples[0].get('sample_number') if isinstance(samples, list) and samples else None
            
            # Ancestry: extract from distribution dict
            ancestry_dist = score.get('ancestry_distribution') or {}
            gwas_data = ancestry_dist.get('gwas') or {}
            dist = gwas_data.get('dist') if isinstance(gwas_data, dict) else {}
            # Get primary ancestry (highest percentage)
            ancestry = max(dist.keys(), key=lambda k: dist.get(k, 0)) if dist else None
            
            cursor.execute("""
                INSERT OR REPLACE INTO polygenic_scores (
                    pgs_id, trait_name, trait_category, trait_efo_id,
                    publication_doi, publication_year, publication_title,
                    sample_size, num_variants, num_variants_downloaded,
                    ancestry, genome_build, download_status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                pgs_id,
                trait[:500] if trait else 'Unknown',
                score.get('trait_category', 'Other'),
                trait_efo_id,
                pub_doi,
                pub_year,
                pub_title,
                sample_size,
                num_variants,
                len(variants),
                ancestry,
                score.get('genome_build'),
                'complete'
            ))
            
            # Delete existing variants for this score (for re-downloads)
            cursor.execute("DELETE FROM pgs_variants WHERE pgs_id = ?", (pgs_id,))
            
            # Insert variants in batches
            batch_size = 1000
            for i in range(0, len(variants), batch_size):
                batch = variants[i:i+batch_size]
                cursor.executemany("""
                    INSERT INTO pgs_variants (
                        pgs_id, rsid, chromosome, position,
                        effect_allele, other_allele, effect_weight, allele_frequency
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, [(
                    pgs_id,
                    v['rsid'],
                    v['chromosome'],
                    v['position'],
                    v['effect_allele'],
                    v['other_allele'],
                    v['weight'],
                    v['frequency']
                ) for v in batch])
            
            conn.commit()
            
            scores_done += 1
            variants_done += len(variants)
            
            # Update state (save periodically)
            state["pgs"]["processed_scores"] = list(set(state["pgs"].get("processed_scores", [])) | {pgs_id})
            if scores_done % 10 == 0:
                save_state(state)
            
        except sqlite3.Error as e:
            scores_failed += 1
            conn.rollback()
        
        # Update progress bar
        pbar.update(num_variants)
        
        # Update ETA in description
        elapsed = time.time() - start_time
        if variants_done > 0:
            rate = variants_done / elapsed
            remaining = remaining_variants - pbar.n
            eta = remaining / rate if rate > 0 else 0
            pbar.set_postfix({
                'scores': scores_done,
                'failed': scores_failed,
                'ETA': format_time(eta)
            })
        
        # Rate limiting
        time.sleep(REQUEST_DELAY)
    
    pbar.close()
    
    # Final state update
    state["pgs"]["last_update"] = datetime.now().isoformat()
    state["pgs"]["status"] = "complete"
    state["pgs"]["total_scores"] = scores_done
    state["pgs"]["total_variants"] = variants_done
    save_state(state)
    
    # Update metadata
    cursor.execute(
        "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
        ("last_update", datetime.now().isoformat())
    )
    conn.commit()
    conn.close()
    
    # Summary
    elapsed = time.time() - start_time
    print(f"\n✅ PGS update complete!")
    print(f"   Scores downloaded: {scores_done:,}")
    print(f"   Scores failed: {scores_failed:,}")
    print(f"   Variants stored: {variants_done:,}")
    print(f"   Time taken: {format_time(elapsed)}")
    
    return True


# =============================================================================
# STATUS
# =============================================================================

def show_status() -> None:
    """Show current database status."""
    state = load_state()
    
    print("\n" + "="*60)
    print("📊 DATABASE STATUS")
    print("="*60)
    
    # GWAS status
    print("\n🧬 GWAS Catalog:")
    gwas_state = state.get("gwas", {})
    if GWAS_DB.exists():
        size = GWAS_DB.stat().st_size
        conn = sqlite3.connect(GWAS_DB)
        cursor = conn.cursor()
        try:
            cursor.execute("SELECT COUNT(*) FROM gwas_associations")
            count = cursor.fetchone()[0]
            cursor.execute("SELECT value FROM metadata WHERE key='last_update'")
            row = cursor.fetchone()
            last_update = row[0] if row else "Unknown"
            
            print(f"   Status: ✅ Available")
            print(f"   Associations: {count:,}")
            print(f"   Database size: {format_size(size)}")
            print(f"   Last update: {last_update}")
        except sqlite3.OperationalError:
            print(f"   Status: ⚠️ Database exists but not initialized")
            print(f"   Run 'python update_databases.py gwas' to initialize")
        conn.close()
    else:
        print(f"   Status: ❌ Not downloaded")
    
    # PGS status
    print("\n📈 PGS Catalog:")
    pgs_state = state.get("pgs", {})
    if PGS_DB.exists():
        size = PGS_DB.stat().st_size
        conn = sqlite3.connect(PGS_DB)
        cursor = conn.cursor()
        try:
            cursor.execute("SELECT COUNT(*) FROM polygenic_scores")
            scores = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM pgs_variants")
            variants = cursor.fetchone()[0]
            cursor.execute("SELECT value FROM metadata WHERE key='last_update'")
            row = cursor.fetchone()
            last_update = row[0] if row else "Unknown"
            
            processed = len(pgs_state.get("processed_scores", []))
            failed = len(pgs_state.get("failed_scores", []))
            
            print(f"   Status: ✅ Available")
            print(f"   Scores downloaded: {scores:,}")
            print(f"   Total variants: {variants:,}")
            print(f"   Database size: {format_size(size)}")
            print(f"   Last update: {last_update}")
            if failed > 0:
                print(f"   ⚠️ Failed downloads: {failed}")
        except sqlite3.OperationalError:
            print(f"   Status: ⚠️ Database exists but not initialized")
            print(f"   Run 'python update_databases.py pgs' to initialize")
        conn.close()
    else:
        print(f"   Status: ❌ Not downloaded")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Update GenExplore genetic databases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python update_databases.py gwas      # Update GWAS only
    python update_databases.py pgs       # Update PGS only  
    python update_databases.py all       # Update both
    python update_databases.py status    # Show current status

Notes:
    - Downloads are resumable (safe to interrupt)
    - PGS scores with >100K variants are skipped (too large)
    - GWAS uses bulk download (~100MB)
    - Full PGS download may take several hours
        """
    )
    
    parser.add_argument(
        'command',
        choices=['gwas', 'pgs', 'all', 'status'],
        help='What to update'
    )
    
    parser.add_argument(
        '--fresh',
        action='store_true',
        help='Start fresh (ignore previous progress)'
    )
    
    args = parser.parse_args()
    
    # Load or initialize state
    state = load_state()
    
    if args.fresh:
        if args.command in ['pgs', 'all']:
            state["pgs"] = {
                "last_update": None,
                "status": "never",
                "processed_scores": [],
                "failed_scores": [],
                "total_scores": 0,
                "total_variants": 0
            }
        if args.command in ['gwas', 'all']:
            state["gwas"] = {"last_update": None, "status": "never"}
        save_state(state)
    
    if args.command == 'status':
        show_status()
        return
    
    print("\n🧬 GenExplore Database Updater")
    print("="*60)
    
    success = True
    
    if args.command in ['gwas', 'all']:
        success = update_gwas(state) and success
    
    if args.command in ['pgs', 'all']:
        success = update_pgs(state, resume=not args.fresh) and success
    
    if success:
        print("\n" + "="*60)
        print("✅ Update complete!")
        show_status()
    else:
        print("\n⚠️ Some updates failed. Check errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
