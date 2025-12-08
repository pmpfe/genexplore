# Genetic Analysis Desktop Application

A PyQt6-based desktop application for analyzing 23andMe genetic data against the GWAS Catalog.

## Features

### Level 1 - MVP
- Upload 23andMe raw data files (.txt format)
- Parse and validate SNP data
- Match user SNPs against GWAS Catalog database
- Display results in sortable, paginated table
- Show: SNP ID, Gene, Trait, User Genotype, Risk Allele, P-value

### Level 2 - Advanced Features
- Impact Score calculation (0-10 scale)
- Real-time filtering:
  - Minimum impact score slider
  - Maximum p-value slider (log scale)
  - Trait category dropdown
  - Free text search (traits, genes, SNPs)
- Sort by any column
- Visual indicators for risk alleles

## Installation

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Setup database (creates GWAS sample data)
python database/setup_database.py

# Run the application
python main.py
```

## Requirements

- Python 3.9+
- PyQt6 6.7.0
- pandas 2.2.0
- numpy 1.24.4

## Project Structure

```
genetic_app/
├── main.py                    # Application entry point
├── config.py                  # Configuration constants
├── requirements.txt           # Python dependencies
├── sample_23andme.txt         # Sample data file for testing
│
├── database/
│   ├── setup_database.py      # Database setup script
│   └── gwas.db               # SQLite database (generated)
│
├── backend/
│   ├── parsers.py            # 23andMe file parser
│   ├── search_engine.py      # GWAS matching and search
│   ├── scoring.py            # Impact score calculation
│   └── validators.py         # Input validation
│
├── frontend/
│   └── main_window.py        # PyQt6 main window
│
├── models/
│   └── data_models.py        # Data classes (SNPRecord, GWASMatch, etc.)
│
├── utils/
│   ├── logging_config.py     # Logging setup
│   └── file_utils.py         # File utilities
│
├── tests/
│   ├── test_parsers.py       # Parser tests
│   ├── test_scoring.py       # Scoring tests
│   └── test_search_engine.py # Search engine tests
│
└── logs/
    └── app.log               # Application log (generated)
```

## Usage

1. **Start the application**: `python main.py`
2. **Upload file**: Click "Upload 23andMe File" and select your raw data file
3. **View results**: Browse matches in the results table
4. **Filter results**:
   - Adjust the "Min Impact Score" slider
   - Adjust the "Max P-value" slider
   - Select a trait category
   - Type in the search box
5. **Sort results**: Click any column header to sort
6. **Navigate pages**: Use Previous/Next buttons for large result sets

## Impact Score Calculation

The impact score (0-10) combines:

1. **P-value significance**: Lower p-values = higher impact
2. **Allele frequency**: Rarer variants = higher impact

```
score_p_value = min(-log10(p_value) / 10, 1.0) × 7.0
score_af = (1 - allele_frequency) × 3.0
impact_score = clamp(score_p_value + score_af, 0, 10)
```

## Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_parsers.py -v

# Run with coverage
python -m pytest tests/ --cov=backend --cov=models
```

## 23andMe File Format

The parser accepts tab-separated files with format:
```
# Comments start with #
rsid	chromosome	position	genotype
rs3131972	1	694713	GG
rs12124819	1	713790	AG
```

- **rsid**: SNP identifier (rs followed by numbers)
- **chromosome**: 1-22, X, Y, or MT
- **position**: Genomic position (positive integer)
- **genotype**: Two alleles (e.g., AA, AT, GG)
- Lines with genotype "--" (undetermined) are skipped

## Database

The GWAS database contains:
- 100+ variant entries from GWAS Catalog
- Categories: Metabolic, Cardiovascular, Neuropsychiatric, Physical Trait, Oncology, Immune, Infectious, Other
- Population allele frequencies (overall, EUR, AFR, EAS, AMR)

To recreate the database:
```bash
python database/setup_database.py --drop
```

## Logging

Logs are written to `logs/app.log`:
- Console: INFO level
- File: DEBUG level with rotation (10MB max, 5 backups)

## Troubleshooting

**Database not found error:**
```bash
python database/setup_database.py
```

**No matches found:**
- Ensure the 23andMe file format is correct
- Check that variants exist in the GWAS database

**UI not responding:**
- The application uses background threading for processing
- Large files may take a few seconds to process

## License

This project is for educational and research purposes.

## Disclaimer

This application is for informational purposes only and should not be used for medical diagnosis or treatment decisions. Always consult with healthcare professionals for genetic counseling.
