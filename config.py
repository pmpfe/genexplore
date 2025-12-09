"""
Configuration constants for the Genetic Analysis Application.
"""

import logging
import os
import re

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATABASE_PATH = os.path.join(BASE_DIR, 'database', 'gwas.db')
PGS_DATABASE_PATH = os.path.join(BASE_DIR, 'database', 'pgs.db')
BACKUP_DIR = os.path.join(BASE_DIR, 'database', 'backups')
LOG_DIR = os.path.join(BASE_DIR, 'logs')

# Logging configuration
LOG_LEVEL = logging.DEBUG
LOG_FILE_SIZE = 10 * 1024 * 1024  # 10MB
LOG_BACKUP_COUNT = 5

# Pagination
RESULTS_PER_PAGE = 50

# Default values for scoring
DEFAULT_ALLELE_FREQUENCY = 0.5
DEFAULT_IMPACT_SCORE = 5.0

# Valid chromosome values
VALID_CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

# Trait categories
TRAIT_CATEGORIES = [
    'ALL',
    'Metabolic',
    'Cardiovascular',
    'Neuropsychiatric',
    'Physical Trait',
    'Infectious',
    'Immune',
    'Oncology',
    'Other'
]

# Polygenic score categories (for PGS analysis)
PGS_CATEGORIES = [
    'Metabolic',
    'Cardiovascular',
    'Neuropsychiatric',
    'Oncology',
    'Immune',
    'Physical Trait',
    'Other'
]

# Validation patterns
RSID_PATTERN = re.compile(r'^rs\d+$')
GENOTYPE_PATTERN = re.compile(r'^[ATCG]{2}$')
GENOTYPE_UNDETERMINED = '--'

# Polygenic scoring thresholds
PGS_LOW_COVERAGE_THRESHOLD = 80.0  # Warn if coverage below this %
PGS_HIGH_RISK_PERCENTILE = 80  # ≥80th percentile = high risk
PGS_LOW_RISK_PERCENTILE = 20   # <20th percentile = low risk

# Application settings
APP_NAME = 'Genetic Analysis'
APP_VERSION = '2.0.0'  # Updated for polygenic integration
