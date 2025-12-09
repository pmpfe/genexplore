# GenExplore - Genetic Analysis Desktop Application

A PyQt6-based desktop application for analyzing 23andMe genetic data, featuring both monogenic (single-variant) analysis against the GWAS Catalog and polygenic risk score (PRS) analysis using the PGS Catalog.

![App Screenshot](utils/screenshot.png)

## Quick Start

```bash
# Setup
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python database/setup_database.py

# Download full databases (optional but recommended)
python database/update_databases.py --all

# Run
python main.py
```

## Features

### 🧬 Monogenic Analysis
- Upload 23andMe raw data files
- Match SNPs against 888K+ GWAS variants
- Impact score calculation (0-10)
- Filtering by score, p-value, category, carrier status
- Detailed explanations with external links

### 📊 Polygenic Risk Scores
- 2,968 polygenic scores covering 660+ traits
- 25M+ variant weights from PGS Catalog
- Background computation with progress indicators
- Population distribution visualization
- Risk categories (Low/Intermediate/High)
- Coverage quality warnings

### ⚡ Performance
- Three-stage parallel loading: File → Monogenic → Polygenic
- Non-blocking UI during calculations
- Efficient SQLite queries with proper indexing

---

## Installation

### Requirements
- Python 3.9+
- ~4GB disk space (for full databases)

### Dependencies
```
PyQt6>=6.7.0
pandas>=2.2.0
numpy>=1.24.4
requests>=2.31.0
tqdm>=4.66.0
```

### Setup Commands
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Setup sample databases
python database/setup_database.py

# Download full databases (recommended)
python database/update_databases.py --all

# Run application
python main.py
```

---

## Database Update Script

The `database/update_databases.py` script downloads and updates the scientific databases.

### Usage
```bash
# Update both GWAS and PGS databases
python database/update_databases.py --all

# Update only GWAS
python database/update_databases.py --gwas

# Update only PGS
python database/update_databases.py --pgs

# Limit PGS scores (for testing)
python database/update_databases.py --pgs --limit 100
```

### What Gets Downloaded

| Database | Source | Records | Size | Time |
|----------|--------|---------|------|------|
| GWAS | GWAS Catalog | 888K variants | ~170MB | ~5 min |
| PGS | PGS Catalog | 2,968 scores, 25M variants | ~3.5GB | ~4-6 hours |

### Resume Capability
- Both downloads support interruption and resume
- Progress saved incrementally to database
- Incomplete scores detected and re-downloaded
- Time estimates shown during download

---

## Usage

### Monogenic Analysis
1. Click "Upload 23andMe File"
2. View matches in results table
3. Use filters (score, p-value, category, search)
4. Click "Explain" for detailed information

### Polygenic Analysis
1. Load 23andMe file (same as above)
2. Switch to "📊 Polygenic Scores" tab
3. Scores compute automatically in background
4. Click "📊 View" for detailed analysis
5. Use filters to find specific traits

### Progress Indicators
When loading a file, three progress bars show:
- **File**: Reading and parsing genetic data
- **Mono**: Matching against GWAS database
- **Poly**: Computing polygenic scores (background)

---

## Scientific Background

### Monogenic Impact Score (0-10)
```
score = p_value_component + allele_frequency_component
      = min(-log10(p_value)/10, 1) × 7 + (1 - AF) × 3
```

### Polygenic Risk Score
```
PRS = Σ (effect_weight × effect_allele_count)
```

Results normalized to population distribution and converted to percentiles.

### Risk Categories
| Percentile | Category |
|------------|----------|
| < 20% | Low Risk |
| 20-79% | Intermediate |
| ≥ 80% | High Risk |

---

## Data Sources

| Source | Description | Link |
|--------|-------------|------|
| GWAS Catalog | Curated GWAS associations | https://www.ebi.ac.uk/gwas/ |
| PGS Catalog | Polygenic score repository | https://www.pgscatalog.org/ |
| dbSNP | SNP reference | https://www.ncbi.nlm.nih.gov/snp/ |

---

## Limitations

1. **Coverage**: Not all PGS variants present in 23andMe data (~40-80% typical)
2. **Population bias**: Most scores derived from European populations
3. **Not diagnostic**: Educational/research purposes only
4. **Static data**: Requires manual database updates

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Database not found | `python database/setup_database.py` |
| No matches found | Check file format, ensure database populated |
| UI freezes | Wait for background tasks, check logs |
| Low coverage warning | Normal - not all variants in genotype file |

---

## Disclaimer

⚠️ **For educational and research purposes only.**

Results should NOT be used for medical diagnosis or treatment. Consult healthcare professionals for interpretation of genetic data.

---

# Technical Documentation

*The following sections are for developers and AI assistants working on this codebase.*

## Project Structure

```
genexplore/
├── main.py                      # Entry point
├── config.py                    # Configuration constants
├── requirements.txt             # Dependencies
├── sample_23andme.txt           # Test data
│
├── database/
│   ├── setup_database.py        # Initial DB setup with sample data
│   ├── polygenic_database.py    # PGS database operations (877 lines)
│   ├── update_databases.py      # Download script (920 lines)
│   ├── gwas.db                  # GWAS SQLite (~170MB)
│   └── pgs.db                   # PGS SQLite (~3.5GB)
│
├── backend/
│   ├── parsers.py               # 23andMe file parser (159 lines)
│   ├── search_engine.py         # GWAS matching engine (310 lines)
│   ├── scoring.py               # Monogenic scoring (138 lines)
│   ├── polygenic_scoring.py     # PRS calculation (324 lines)
│   └── validators.py            # Input validation (170 lines)
│
├── frontend/
│   ├── main_window.py           # Main UI, tabs, monogenic (1400+ lines)
│   └── polygenic_widgets.py     # Polygenic UI components (1200+ lines)
│
├── models/
│   ├── data_models.py           # Monogenic dataclasses
│   └── polygenic_models.py      # Polygenic dataclasses
│
├── utils/
│   ├── logging_config.py        # Logging setup
│   └── file_utils.py            # File utilities
│
├── tests/                       # pytest tests
└── logs/                        # Runtime logs
```

## Database Schema

### gwas.db
```sql
gwas_associations (
    rsid TEXT PRIMARY KEY,
    gene TEXT,
    trait TEXT,
    risk_allele TEXT,
    p_value REAL,
    odds_ratio REAL,
    category TEXT,
    af_overall REAL, af_eur REAL, af_afr REAL, af_eas REAL, af_amr REAL
)
```

### pgs.db
```sql
polygenic_scores (
    pgs_id TEXT PRIMARY KEY,
    trait TEXT,
    publication TEXT,
    num_variants INTEGER,
    category TEXT,
    ancestry TEXT
)

pgs_variants (
    id INTEGER PRIMARY KEY,
    pgs_id TEXT,
    rsid TEXT,
    effect_allele TEXT,
    effect_weight REAL,
    FOREIGN KEY (pgs_id) REFERENCES polygenic_scores(pgs_id)
)
-- Indexes on pgs_id and rsid for performance

population_distributions (
    pgs_id TEXT PRIMARY KEY,
    mean REAL,
    std REAL,
    percentiles TEXT  -- JSON
)
```

## Key Classes

### Backend
- `GeneticDataParser`: Parses 23andMe files → dict[rsid, genotype]
- `SearchEngine`: Matches user SNPs against GWAS DB
- `PolygenicScoringEngine`: Calculates PRS scores
- `PolygenicDatabase`: PGS database operations

### Frontend
- `MainWindow`: Main application window with tabs
- `PolygenicTab`: Polygenic scores browser
- `PolygenicDetailDialog`: Score detail view with distribution plot

### Threading
- `FileLoadWorker`: Background file loading
- `MonogenicComputeWorker`: GWAS matching
- `PolygenicComputeWorker`: PRS calculation (non-blocking)

## Performance Considerations

1. **Database indexing**: rsid indexed in both databases
2. **Batch queries**: Variants fetched in chunks
3. **Background computation**: UI remains responsive
4. **Progress signals**: Qt signals for UI updates

## Update Script Architecture

The `update_databases.py` script:

1. **GWAS Update**:
   - Downloads TSV from GWAS Catalog FTP
   - Streams and parses incrementally
   - Inserts with batched transactions

2. **PGS Update**:
   - Fetches score metadata from PGS Catalog API
   - Downloads scoring files (.txt.gz) individually
   - Parses variant weights
   - Detects and cleans incomplete scores on resume
   - Progress based on variant count (weighted)

## Testing

```bash
python -m pytest tests/ -v
python -m pytest tests/test_polygenic_scoring.py -v
```

## Logging

- Console: INFO level
- File (`logs/app.log`): DEBUG, 10MB rotation, 5 backups
- Errors: Also to `error.log`

## Future Improvements

1. **Ancestry-specific scoring**: Use population-matched distributions
2. **Score quality metrics**: Incorporate PGS Catalog quality indicators
3. **Automatic updates**: Scheduled database refresh
4. **Export functionality**: PDF reports, CSV export
5. **Additional file formats**: Ancestry, MyHeritage support
