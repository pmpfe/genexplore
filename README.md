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

### 💾 Save & Load Sessions
- Save complete analysis to compressed .gxs files
- Load previous sessions instantly
- No re-computation needed when loading
- Portable session files (~2-5 MB)

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

### Save & Load Sessions
Analysis results can be saved and loaded to avoid re-computation:

1. **Save**: After analysis completes, click "💾 Save" to export results
2. **Load**: Click "📂 Load" to restore a previously saved session

**Session files (.gxs)**:
- Compressed JSON format (~2-5 MB per session)
- Contains: SNP records, GWAS matches, polygenic scores
- Loads instantly without re-computation
- Portable between installations

### Progress Indicators
When loading a file, three progress bars show:
- **File**: Reading and parsing genetic data
- **Mono**: Matching against GWAS database
- **Poly**: Computing polygenic scores (background)

---

## Scientific Background & Calculation Methods

### Monogenic Impact Score (0-10)

The monogenic impact score combines two components to rank variant significance:

```
impact_score = p_value_component + allele_frequency_component

Where:
  p_value_component = min(-log10(p_value) / 10, 1.0) × 7.0
  allele_frequency_component = (1 - allele_frequency) × 3.0

Final score clamped to [0, 10]
```

**Rationale**:
- **P-value component (0-7 points)**: More significant associations (lower p-values) score higher. A p-value of 10^-10 gives maximum 7 points.
- **Allele frequency component (0-3 points)**: Rarer variants score higher, as rare risk alleles often have larger effects.

**Score Interpretation**:
| Score | Interpretation |
|-------|----------------|
| ≥ 8.0 | Very High Impact |
| 6.0-7.9 | High Impact |
| 4.0-5.9 | Moderate Impact |
| 2.0-3.9 | Low Impact |
| < 2.0 | Minimal Impact |

### Polygenic Risk Score (PRS) Calculation

PRS aggregates the effects of many variants, each with small individual effect:

```
PRS = Σ (effect_weight × effect_allele_count)

Where:
  effect_weight = β coefficient from GWAS/PGS study
  effect_allele_count = 0, 1, or 2 (copies of effect allele in genotype)
```

**Allele Counting Process**:
1. For each PGS variant, get user's genotype (e.g., "AG")
2. Count how many copies match the effect allele
3. Handle strand flips using complement mapping (A↔T, C↔G)
4. If genotype doesn't match expected alleles, variant is skipped

### Population Distribution Estimation

When pre-computed population distributions are not available (most cases), we estimate them using Hardy-Weinberg equilibrium theory:

```
For each variant with effect weight β and effect allele frequency p:

Expected contribution:  E[X] = 2 × p × β
Variance contribution:  Var[X] = 2 × p × (1-p) × β²

Population mean:  μ = Σ E[X]  (sum over all variants)
Population std:   σ = √(Σ Var[X])
```

**Approximations Used**:
1. **Missing allele frequencies**: When effect allele frequency is not available in PGS data, we assume p = 0.5 (maximum uncertainty)
2. **Independence assumption**: Variants are assumed to be independent (no linkage disequilibrium correction)
3. **Hardy-Weinberg equilibrium**: Assumes random mating population

### Percentile Calculation

Raw scores are converted to percentiles using z-score normalization:

```
z_score = (raw_score - population_mean) / population_std
percentile = Φ(z_score) × 100

Where Φ is the standard normal CDF
```

**Implementation**: Uses `percentiles` dict if available, otherwise linear interpolation or normal approximation.

### Risk Categories
| Percentile | Category | Interpretation |
|------------|----------|----------------|
| < 20% | Low Risk | Lower genetic predisposition than 80% of population |
| 20-79% | Intermediate | Within average range |
| ≥ 80% | High Risk | Higher genetic predisposition than 80% of population |

### Coverage Quality

Coverage = (variants_found / variants_total) × 100%

| Coverage | Quality | Notes |
|----------|---------|-------|
| ≥ 70% | Good | Results are reliable |
| 50-69% | Moderate | Results should be interpreted with caution |
| < 50% | Low | ⚠️ Warning displayed, results may be unreliable |

**Why coverage varies**: 23andMe genotyping chips don't include all variants used in PGS studies. Typical coverage is 40-80% depending on the score.

---

## Data Sources

| Source | Description | Link |
|--------|-------------|------|
| GWAS Catalog | Curated GWAS associations | https://www.ebi.ac.uk/gwas/ |
| PGS Catalog | Polygenic score repository | https://www.pgscatalog.org/ |
| dbSNP | SNP reference | https://www.ncbi.nlm.nih.gov/snp/ |

---

## Limitations & Approximations

### Data Coverage Limitations
1. **Variant coverage**: 23andMe chips include ~600K-700K SNPs. PGS scores may require variants not on the chip (40-80% coverage typical)
2. **Missing allele frequencies**: When not provided in PGS data, p=0.5 is assumed, which may over/underestimate variance
3. **No imputation**: Missing variants are simply skipped, not imputed from nearby variants

### Statistical Approximations
1. **Independence assumption**: Variants are treated as independent; linkage disequilibrium not corrected
2. **Hardy-Weinberg equilibrium**: Population distribution estimates assume HWE
3. **Normal distribution**: Percentiles assume scores are normally distributed in the population
4. **Single ancestry**: No adjustment for ancestry-specific allele frequencies

### Scientific Limitations
1. **Population bias**: Most GWAS/PGS studies are from European populations; accuracy may be lower for other ancestries
2. **Environmental factors**: Polygenic scores don't account for lifestyle, diet, or environmental exposures
3. **Gene-gene interactions**: Epistatic effects are not modeled
4. **Rare variants**: Focus on common variants; rare high-impact variants may be missed

### Practical Limitations
1. **Not diagnostic**: Results are for educational/research purposes only
2. **Static data**: Databases require manual updates via `update_databases.py`
3. **No clinical validation**: Scores not validated for clinical use

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
│   ├── session_manager.py       # Save/Load sessions (280 lines)
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

## Algorithm Details

### Allele Counting with Strand Flip Handling

When matching user genotype to PGS variant:

```python
# 1. Check if genotype alleles match expected
if allele1 in {effect_allele, other_allele}:
    count directly
else:
    # 2. Try complement (strand flip)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if complement[allele1] in {effect_allele, other_allele}:
        use complemented alleles
    else:
        skip variant (ambiguous)
```

### Population Distribution Estimation Algorithm

```python
for each variant in score:
    p = effect_allele_frequency or 0.5  # Default if missing
    β = effect_weight
    
    mean += 2 * p * β           # Expected diploid contribution
    variance += 2 * p * (1-p) * β²  # Binomial variance

std = sqrt(variance)
```

### Percentile from Z-Score

```python
z_score = (raw_score - mean) / std
percentile = norm.cdf(z_score) * 100  # Standard normal CDF
```

### Three-Stage Loading Process

```
Stage 1: File Loading (blocking)
├── Parse 23andMe file format
├── Validate SNP records
└── Build genotype lookup dict {rsid: genotype}

Stage 2: Monogenic Analysis (blocking)
├── Query GWAS database for matching rsids
├── Calculate impact scores
├── Sort by impact score
└── Display in results table

Stage 3: Polygenic Analysis (background, non-blocking)
├── Load PGS score definitions from database
├── For each score:
│   ├── Fetch variants from database
│   ├── Match to user genotypes
│   ├── Compute weighted sum
│   ├── Estimate population distribution
│   ├── Calculate percentile
│   └── Assign risk category
├── Emit progress signals to UI
└── Display results when complete
```

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

## Session File Format (.gxs)

Session files are gzip-compressed JSON with the following structure:

```json
{
  "format_version": "1.0",
  "created_at": "2024-01-15T10:30:00",
  "metadata": { "app_version": "1.0.0" },
  "summary": {
    "snp_count": 700000,
    "gwas_match_count": 1500,
    "polygenic_score_count": 660
  },
  "snp_records": [
    {"rsid": "rs123", "chromosome": "1", "position": 12345, "genotype": "AG"}
  ],
  "gwas_matches": [
    {"rsid": "...", "trait": "...", "impact_score": 7.5, ...}
  ],
  "polygenic_results": [
    {"pgs_id": "PGS000001", "trait_name": "...", "percentile": 65.2, ...}
  ]
}
```

Typical compressed sizes: 2-5 MB per session.
