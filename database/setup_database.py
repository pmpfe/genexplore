"""
Database setup script for the Genetic Analysis Application.

Creates and populates the GWAS database with sample data.
"""

import sqlite3
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import DATABASE_PATH

# Sample GWAS data - Real variants from GWAS Catalog
SAMPLE_GWAS_DATA = [
    # Oncology
    ("rs6983267", "8", 128413305, "A", "G", "G", "Colorectal cancer", "MYC", 5.2e-11, 1.27, 201518, "Oncology", "17618284"),
    ("rs1042522", "17", 7579472, "C", "G", "C", "Lung cancer", "TP53", 3.1e-8, 1.18, 15000, "Oncology", "19336370"),
    ("rs2981582", "10", 123337335, "A", "G", "A", "Breast cancer", "FGFR2", 2.0e-76, 1.26, 96000, "Oncology", "17529967"),
    ("rs889312", "5", 56031884, "A", "C", "C", "Breast cancer", "MAP3K1", 6.5e-20, 1.13, 96000, "Oncology", "17529967"),
    ("rs13281615", "8", 128355618, "A", "G", "G", "Breast cancer", "MYC", 5.1e-9, 1.08, 96000, "Oncology", "17529967"),
    ("rs4430796", "17", 36098040, "A", "G", "A", "Prostate cancer", "HNF1B", 1.4e-10, 1.22, 50000, "Oncology", "17603485"),
    ("rs1447295", "8", 128554220, "A", "C", "A", "Prostate cancer", "CASC8", 7.7e-14, 1.72, 50000, "Oncology", "17401363"),
    ("rs10993994", "10", 51549496, "C", "T", "T", "Prostate cancer", "MSMB", 7.8e-31, 1.25, 50000, "Oncology", "18264097"),
    ("rs6504950", "17", 53056471, "A", "G", "G", "Breast cancer", "STXBP4", 4.1e-8, 1.05, 96000, "Oncology", "20453838"),
    ("rs614367", "11", 69328764, "C", "T", "T", "Breast cancer", "CCND1", 1.4e-13, 1.15, 96000, "Oncology", "20453838"),
    
    # Cardiovascular
    ("rs1333049", "9", 22125500, "C", "G", "C", "Coronary artery disease", "CDKN2A", 4.8e-14, 1.36, 256000, "Cardiovascular", "17554300"),
    ("rs10757278", "9", 22124478, "A", "G", "G", "Myocardial infarction", "CDKN2A", 1.2e-20, 1.28, 50000, "Cardiovascular", "17478679"),
    ("rs1801133", "1", 11856378, "C", "T", "T", "Homocysteine levels", "MTHFR", 2.3e-15, 1.15, 30000, "Cardiovascular", "18439552"),
    ("rs662799", "11", 116663707, "A", "G", "A", "Triglycerides", "APOA5", 1.0e-50, 1.25, 100000, "Cardiovascular", "18193043"),
    ("rs964184", "11", 116648917, "C", "G", "G", "LDL cholesterol", "ZNF259", 3.4e-40, 1.15, 100000, "Cardiovascular", "20686565"),
    ("rs12740374", "1", 109821511, "G", "T", "T", "LDL cholesterol", "CELSR2", 1.0e-170, 1.20, 100000, "Cardiovascular", "20686565"),
    ("rs6511720", "19", 11202306, "G", "T", "T", "LDL cholesterol", "LDLR", 5.0e-117, 1.25, 100000, "Cardiovascular", "20686565"),
    ("rs1800961", "20", 43042364, "C", "T", "T", "HDL cholesterol", "HNF4A", 1.0e-14, 1.10, 100000, "Cardiovascular", "18193044"),
    ("rs4420638", "19", 45422946, "A", "G", "G", "Total cholesterol", "APOE", 1.0e-300, 1.30, 100000, "Cardiovascular", "18193043"),
    ("rs1799983", "7", 150999023, "G", "T", "T", "Hypertension", "NOS3", 2.5e-8, 1.12, 50000, "Cardiovascular", "19430483"),
    ("rs5186", "3", 148459988, "A", "C", "C", "Hypertension", "AGTR1", 1.8e-9, 1.15, 50000, "Cardiovascular", "20522523"),
    ("rs17465637", "1", 56962821, "A", "C", "C", "Coronary artery disease", "MIA3", 1.2e-9, 1.14, 100000, "Cardiovascular", "17634449"),
    
    # Metabolic
    ("rs7903146", "10", 114758349, "C", "T", "T", "Type 2 diabetes", "TCF7L2", 2.3e-36, 1.37, 120000, "Metabolic", "17293876"),
    ("rs7903146", "10", 114758349, "C", "T", "T", "Fasting glucose", "TCF7L2", 1.1e-20, 1.25, 80000, "Metabolic", "20081858"),
    ("rs7903146", "10", 114758349, "C", "T", "T", "Insulin resistance", "TCF7L2", 5.5e-15, 1.20, 60000, "Metabolic", "22158537"),
    ("rs1801282", "3", 12393125, "C", "G", "C", "Type 2 diabetes", "PPARG", 1.7e-6, 1.14, 80000, "Metabolic", "18372903"),
    ("rs5219", "11", 17409572, "C", "T", "T", "Type 2 diabetes", "KCNJ11", 5.0e-11, 1.14, 80000, "Metabolic", "18372903"),
    ("rs13266634", "8", 118184783, "C", "T", "T", "Type 2 diabetes", "SLC30A8", 5.3e-8, 1.12, 80000, "Metabolic", "17460697"),
    ("rs10811661", "9", 22134095, "C", "T", "T", "Type 2 diabetes", "CDKN2A", 7.8e-15, 1.20, 80000, "Metabolic", "17463246"),
    ("rs7756992", "6", 20679709, "A", "G", "G", "Type 2 diabetes", "CDKAL1", 4.1e-11, 1.12, 80000, "Metabolic", "17463249"),
    ("rs9939609", "16", 53820527, "A", "T", "A", "Obesity", "FTO", 1.0e-42, 1.31, 200000, "Metabolic", "17434869"),
    ("rs17782313", "18", 57851097, "C", "T", "C", "BMI", "MC4R", 2.0e-15, 1.12, 200000, "Metabolic", "18454148"),
    ("rs1558902", "16", 53803574, "A", "T", "A", "BMI", "FTO", 4.8e-120, 1.39, 339000, "Metabolic", "25673413"),
    ("rs2943641", "2", 227093745, "C", "T", "C", "Insulin resistance", "IRS1", 5.4e-20, 1.19, 50000, "Metabolic", "22158537"),
    ("rs780094", "2", 27741237, "C", "T", "T", "Fasting glucose", "GCKR", 1.0e-30, 1.10, 80000, "Metabolic", "20081858"),
    ("rs560887", "2", 169763148, "C", "T", "C", "Fasting glucose", "G6PC2", 1.0e-75, 1.20, 80000, "Metabolic", "20081858"),
    ("rs174547", "11", 61570783, "C", "T", "T", "Triglycerides", "FADS1", 4.5e-24, 1.08, 100000, "Metabolic", "20686565"),
    ("rs328", "8", 19819724, "C", "G", "G", "Triglycerides", "LPL", 1.0e-27, 1.15, 100000, "Metabolic", "18193043"),
    ("rs1260326", "2", 27730940, "C", "T", "T", "Triglycerides", "GCKR", 8.0e-133, 1.15, 100000, "Metabolic", "20686565"),
    ("rs429358", "19", 45411941, "C", "T", "C", "Alzheimer disease", "APOE", 1.0e-200, 3.68, 50000, "Metabolic", "19734902"),
    ("rs7412", "19", 45412079, "C", "T", "T", "LDL cholesterol", "APOE", 1.0e-150, 1.40, 100000, "Metabolic", "20686565"),
    
    # Neuropsychiatric
    ("rs1006737", "3", 53127857, "A", "G", "A", "Bipolar disorder", "CACNA1C", 7.0e-8, 1.18, 50000, "Neuropsychiatric", "18711365"),
    ("rs2251219", "16", 9975495, "C", "T", "T", "Major depression", "GRIN2A", 6.0e-8, 1.12, 50000, "Neuropsychiatric", "21926974"),
    ("rs6265", "11", 27679916, "C", "T", "T", "Depression", "BDNF", 3.0e-6, 1.10, 40000, "Neuropsychiatric", "21112890"),
    ("rs1800497", "11", 113400106, "C", "T", "T", "Alcohol dependence", "DRD2", 2.5e-7, 1.25, 30000, "Neuropsychiatric", "18227835"),
    ("rs4680", "22", 19963748, "A", "G", "G", "Schizophrenia", "COMT", 4.0e-6, 1.08, 40000, "Neuropsychiatric", "21926972"),
    ("rs1344706", "2", 185778428, "A", "C", "C", "Schizophrenia", "ZNF804A", 1.6e-7, 1.10, 40000, "Neuropsychiatric", "18711365"),
    ("rs12807809", "11", 113412746, "C", "T", "T", "Nicotine dependence", "NCAM1", 1.3e-8, 1.15, 50000, "Neuropsychiatric", "20418890"),
    ("rs16969968", "15", 78882925, "A", "G", "A", "Nicotine dependence", "CHRNA5", 5.0e-19, 1.32, 50000, "Neuropsychiatric", "18385739"),
    ("rs6313", "13", 47471478, "C", "T", "T", "Depression response", "HTR2A", 1.0e-5, 1.15, 20000, "Neuropsychiatric", "18073774"),
    ("rs10494561", "1", 78433414, "A", "G", "G", "Anxiety disorders", "PTBP2", 3.2e-6, 1.12, 30000, "Neuropsychiatric", "26754954"),
    ("rs12124819", "1", 713790, "A", "G", "G", "Cognitive function", "LINC01128", 2.1e-8, 1.05, 100000, "Neuropsychiatric", "25869804"),
    
    # Physical Traits
    ("rs12913832", "15", 28365618, "A", "G", "G", "Eye color", "HERC2", 1.0e-300, 25.0, 10000, "Physical Trait", "18252222"),
    ("rs1667394", "15", 28530182, "A", "G", "A", "Eye color", "OCA2", 3.5e-40, 5.0, 10000, "Physical Trait", "18252222"),
    ("rs1426654", "15", 48426484, "A", "G", "A", "Skin pigmentation", "SLC24A5", 1.0e-200, 10.0, 5000, "Physical Trait", "16357253"),
    ("rs16891982", "5", 33951693, "C", "G", "C", "Skin pigmentation", "SLC45A2", 5.0e-90, 8.0, 5000, "Physical Trait", "17999355"),
    ("rs1805007", "16", 89986117, "C", "T", "T", "Red hair", "MC1R", 1.0e-150, 20.0, 10000, "Physical Trait", "11692016"),
    ("rs1805008", "16", 89986144, "C", "T", "T", "Red hair", "MC1R", 2.0e-80, 12.0, 10000, "Physical Trait", "11692016"),
    ("rs143384", "20", 34025756, "A", "G", "A", "Height", "GDF5", 5.0e-30, 1.08, 250000, "Physical Trait", "20881960"),
    ("rs1042725", "12", 66339827, "C", "T", "T", "Height", "HMGA2", 4.2e-16, 1.05, 250000, "Physical Trait", "18391950"),
    ("rs6060369", "20", 6604619, "C", "T", "T", "Height", "BMP2", 7.5e-11, 1.04, 250000, "Physical Trait", "20881960"),
    ("rs2284746", "2", 55308273, "C", "G", "G", "Height", "EFEMP1", 1.0e-21, 1.06, 250000, "Physical Trait", "20881960"),
    ("rs6060373", "20", 6608684, "G", "T", "T", "Height", "BMP2", 2.3e-12, 1.04, 250000, "Physical Trait", "20881960"),
    ("rs1799971", "6", 154039662, "A", "G", "G", "Pain sensitivity", "OPRM1", 1.5e-6, 1.20, 10000, "Physical Trait", "15057820"),
    ("rs4988235", "2", 136608646, "C", "T", "T", "Lactose intolerance", "MCM6", 1.0e-100, 15.0, 50000, "Physical Trait", "12594458"),
    ("rs713598", "7", 141972604, "C", "G", "C", "Bitter taste perception", "TAS2R38", 1.0e-50, 5.0, 5000, "Physical Trait", "12595690"),
    
    # Immune/Infectious
    ("rs3087243", "2", 204738919, "A", "G", "G", "Rheumatoid arthritis", "CTLA4", 1.0e-8, 1.12, 30000, "Immune", "17804836"),
    ("rs2476601", "1", 114377568, "A", "G", "A", "Rheumatoid arthritis", "PTPN22", 9.0e-25, 1.75, 30000, "Immune", "15208781"),
    ("rs6897932", "5", 35910332, "C", "T", "C", "Multiple sclerosis", "IL7R", 2.9e-7, 1.18, 20000, "Immune", "17660530"),
    ("rs3135388", "6", 32681631, "A", "G", "A", "Multiple sclerosis", "HLA-DRA", 1.0e-60, 2.50, 20000, "Immune", "17660530"),
    ("rs11209026", "1", 67705958, "A", "G", "A", "Psoriasis", "IL23R", 4.0e-15, 1.40, 25000, "Immune", "17554261"),
    ("rs2201841", "1", 67649663, "C", "T", "T", "Crohn disease", "IL23R", 2.0e-12, 1.30, 25000, "Immune", "17435756"),
    ("rs17234657", "5", 40437946, "G", "T", "T", "Crohn disease", "PTGER4", 1.0e-10, 1.25, 25000, "Immune", "17435756"),
    ("rs11465804", "1", 67703015, "G", "T", "T", "Ulcerative colitis", "IL23R", 3.0e-9, 1.35, 25000, "Immune", "18587394"),
    ("rs2066847", "16", 50745926, "C", "T", "T", "Crohn disease", "NOD2", 1.0e-40, 2.50, 25000, "Immune", "11385576"),
    ("rs10781499", "9", 117552851, "A", "G", "A", "Inflammatory bowel disease", "CARD9", 8.0e-12, 1.20, 25000, "Immune", "17435756"),
    ("rs2395029", "6", 31431780, "G", "T", "G", "HIV progression", "HLA-B", 5.0e-16, 2.10, 10000, "Infectious", "17767157"),
    ("rs334", "11", 5227002, "A", "T", "T", "Malaria resistance", "HBB", 1.0e-100, 10.0, 5000, "Infectious", "12364793"),
    
    # Other
    ("rs5030737", "10", 54531242, "A", "G", "A", "Vitamin D deficiency", "CYP2R1", 1.5e-10, 1.15, 30000, "Other", "20541252"),
    ("rs12785878", "11", 71167449, "G", "T", "T", "Vitamin D levels", "NADSYN1", 3.0e-12, 1.12, 30000, "Other", "20541252"),
    ("rs2282679", "4", 72608383, "A", "C", "C", "Vitamin D levels", "GC", 1.0e-50, 1.30, 30000, "Other", "20541252"),
    ("rs12203592", "6", 396321, "C", "T", "T", "Freckling", "IRF4", 1.0e-80, 3.00, 10000, "Other", "18488028"),
    ("rs1800562", "6", 26093141, "A", "G", "A", "Hemochromatosis", "HFE", 1.0e-100, 8.00, 20000, "Other", "8696333"),
    ("rs855791", "22", 37462936, "A", "G", "A", "Iron levels", "TMPRSS6", 1.0e-40, 1.25, 50000, "Other", "19553259"),
    ("rs4820268", "22", 37470224, "A", "G", "G", "Hemoglobin levels", "TMPRSS6", 5.0e-35, 1.20, 50000, "Other", "19553259"),
    ("rs1805087", "1", 237048500, "A", "G", "G", "Folate levels", "MTR", 2.0e-10, 1.15, 30000, "Other", "18463370"),
    ("rs12272669", "11", 5248232, "A", "G", "A", "Fetal hemoglobin", "HBG2", 1.0e-45, 1.50, 10000, "Other", "17903302"),
    ("rs2228570", "12", 48272895, "C", "T", "T", "Bone mineral density", "VDR", 3.5e-7, 1.10, 30000, "Other", "19079261"),
    ("rs1801131", "1", 11854476, "A", "C", "C", "Folate metabolism", "MTHFR", 1.0e-8, 1.08, 30000, "Other", "18439552"),
    ("rs4654748", "1", 23754255, "C", "T", "T", "Vitamin B6 levels", "NBPF3", 2.0e-12, 1.12, 20000, "Other", "21878437"),
    ("rs602662", "19", 49206674, "A", "G", "G", "Vitamin B12 levels", "FUT2", 1.0e-25, 1.20, 30000, "Other", "21878437"),
    ("rs3760775", "15", 78806023, "A", "G", "A", "Aging", "CHRNA3", 4.0e-8, 1.08, 50000, "Other", "25740864"),
    ("rs2187668", "6", 32610875, "C", "T", "T", "Celiac disease", "HLA-DQA1", 1.0e-200, 6.00, 30000, "Other", "18311140"),
    ("rs6822844", "4", 123377980, "G", "T", "T", "Celiac disease", "IL2", 1.0e-12, 1.30, 30000, "Other", "18311140"),
    
    # Additional variants for diversity
    ("rs4148323", "2", 234668879, "A", "G", "A", "Bilirubin levels", "UGT1A1", 1.0e-300, 2.50, 50000, "Metabolic", "18179887"),
    ("rs6025", "1", 169549811, "C", "T", "T", "Venous thrombosis", "F5", 1.0e-50, 5.00, 30000, "Cardiovascular", "7989264"),
    ("rs1799963", "11", 46739505, "A", "G", "A", "Venous thrombosis", "F2", 1.0e-20, 2.80, 30000, "Cardiovascular", "8619974"),
    ("rs4149056", "12", 21331549, "C", "T", "C", "Statin response", "SLCO1B1", 5.0e-20, 4.50, 10000, "Other", "18650507"),
    ("rs1057910", "10", 96702047, "A", "C", "C", "Warfarin dose", "CYP2C9", 1.0e-30, 1.80, 10000, "Other", "15930419"),
    ("rs9923231", "16", 31107689, "C", "T", "T", "Warfarin dose", "VKORC1", 1.0e-100, 2.50, 10000, "Other", "15930419"),
    ("rs12979860", "19", 39738787, "C", "T", "C", "Hepatitis C response", "IFNL3", 1.0e-30, 2.00, 10000, "Infectious", "20639878"),
    ("rs8099917", "19", 39733758, "G", "T", "T", "Hepatitis C response", "IFNL3", 5.0e-25, 1.80, 10000, "Infectious", "19749758"),
    ("rs3131972", "1", 694713, "A", "G", "G", "Immune function", "SAMD11", 2.0e-6, 1.05, 50000, "Immune", "21833088"),
    ("rs11240777", "1", 856331, "A", "G", "G", "Gene expression", "KLHL17", 1.5e-5, 1.03, 40000, "Other", "22446963"),
    ("rs4970383", "1", 1019440, "A", "G", "G", "Blood pressure", "ISG15", 3.0e-6, 1.06, 80000, "Cardiovascular", "21909115"),
    ("rs6681049", "1", 909917, "C", "T", "T", "Lipid levels", "PLEKHN1", 2.5e-5, 1.04, 60000, "Metabolic", "20686565"),
]

# Sample allele frequency data
SAMPLE_AF_DATA = [
    ("rs6983267", 0.48, 0.50, 0.42, 0.38, 0.45),
    ("rs1042522", 0.72, 0.70, 0.55, 0.60, 0.65),
    ("rs2981582", 0.38, 0.40, 0.35, 0.30, 0.35),
    ("rs889312", 0.28, 0.30, 0.25, 0.20, 0.27),
    ("rs13281615", 0.40, 0.42, 0.38, 0.35, 0.40),
    ("rs4430796", 0.48, 0.50, 0.45, 0.42, 0.48),
    ("rs1447295", 0.11, 0.12, 0.10, 0.08, 0.10),
    ("rs10993994", 0.41, 0.45, 0.30, 0.35, 0.40),
    ("rs6504950", 0.27, 0.28, 0.25, 0.22, 0.26),
    ("rs614367", 0.16, 0.18, 0.12, 0.10, 0.14),
    ("rs1333049", 0.47, 0.50, 0.40, 0.35, 0.45),
    ("rs10757278", 0.48, 0.50, 0.42, 0.38, 0.46),
    ("rs1801133", 0.35, 0.38, 0.12, 0.30, 0.45),
    ("rs662799", 0.15, 0.08, 0.35, 0.30, 0.12),
    ("rs964184", 0.18, 0.15, 0.40, 0.35, 0.20),
    ("rs12740374", 0.22, 0.25, 0.08, 0.15, 0.20),
    ("rs6511720", 0.12, 0.15, 0.02, 0.08, 0.10),
    ("rs1800961", 0.05, 0.04, 0.02, 0.08, 0.03),
    ("rs4420638", 0.17, 0.20, 0.05, 0.08, 0.15),
    ("rs1799983", 0.33, 0.35, 0.25, 0.15, 0.30),
    ("rs5186", 0.28, 0.30, 0.05, 0.08, 0.22),
    ("rs17465637", 0.72, 0.75, 0.55, 0.60, 0.70),
    ("rs7903146", 0.30, 0.35, 0.28, 0.05, 0.25),
    ("rs1801282", 0.12, 0.15, 0.02, 0.04, 0.10),
    ("rs5219", 0.35, 0.38, 0.42, 0.45, 0.40),
    ("rs13266634", 0.70, 0.75, 0.40, 0.60, 0.65),
    ("rs10811661", 0.79, 0.82, 0.90, 0.55, 0.75),
    ("rs7756992", 0.31, 0.35, 0.42, 0.50, 0.38),
    ("rs9939609", 0.42, 0.45, 0.12, 0.15, 0.35),
    ("rs17782313", 0.24, 0.28, 0.18, 0.20, 0.22),
    ("rs1558902", 0.41, 0.45, 0.12, 0.15, 0.35),
    ("rs2943641", 0.63, 0.65, 0.70, 0.85, 0.68),
    ("rs780094", 0.39, 0.42, 0.55, 0.48, 0.45),
    ("rs560887", 0.30, 0.35, 0.08, 0.20, 0.28),
    ("rs174547", 0.33, 0.35, 0.55, 0.65, 0.40),
    ("rs328", 0.09, 0.10, 0.02, 0.05, 0.08),
    ("rs1260326", 0.40, 0.45, 0.48, 0.55, 0.42),
    ("rs429358", 0.14, 0.15, 0.22, 0.08, 0.10),
    ("rs7412", 0.08, 0.08, 0.10, 0.12, 0.05),
    ("rs1006737", 0.33, 0.35, 0.28, 0.30, 0.32),
    ("rs2251219", 0.38, 0.40, 0.35, 0.32, 0.36),
    ("rs6265", 0.20, 0.18, 0.45, 0.50, 0.25),
    ("rs1800497", 0.19, 0.18, 0.50, 0.35, 0.22),
    ("rs4680", 0.48, 0.50, 0.28, 0.35, 0.45),
    ("rs1344706", 0.32, 0.35, 0.45, 0.40, 0.38),
    ("rs12807809", 0.40, 0.42, 0.35, 0.30, 0.38),
    ("rs16969968", 0.35, 0.38, 0.05, 0.03, 0.15),
    ("rs6313", 0.43, 0.45, 0.52, 0.48, 0.40),
    ("rs10494561", 0.28, 0.30, 0.22, 0.18, 0.25),
    ("rs12124819", 0.25, 0.28, 0.15, 0.20, 0.22),
    ("rs12913832", 0.78, 0.80, 0.02, 0.05, 0.45),
    ("rs1667394", 0.70, 0.75, 0.02, 0.03, 0.35),
    ("rs1426654", 0.99, 0.99, 0.03, 0.05, 0.55),
    ("rs16891982", 0.95, 0.98, 0.02, 0.05, 0.45),
    ("rs1805007", 0.08, 0.10, 0.00, 0.00, 0.03),
    ("rs1805008", 0.05, 0.06, 0.00, 0.00, 0.02),
    ("rs143384", 0.58, 0.60, 0.48, 0.45, 0.55),
    ("rs1042725", 0.49, 0.52, 0.38, 0.42, 0.48),
    ("rs6060369", 0.35, 0.38, 0.28, 0.25, 0.32),
    ("rs2284746", 0.60, 0.62, 0.70, 0.75, 0.58),
    ("rs6060373", 0.33, 0.35, 0.25, 0.22, 0.30),
    ("rs1799971", 0.15, 0.12, 0.02, 0.40, 0.18),
    ("rs4988235", 0.50, 0.75, 0.15, 0.02, 0.35),
    ("rs713598", 0.45, 0.48, 0.35, 0.30, 0.42),
    ("rs3087243", 0.42, 0.45, 0.68, 0.55, 0.48),
    ("rs2476601", 0.09, 0.10, 0.01, 0.00, 0.03),
    ("rs6897932", 0.23, 0.25, 0.15, 0.12, 0.20),
    ("rs3135388", 0.13, 0.15, 0.05, 0.02, 0.08),
    ("rs11209026", 0.07, 0.08, 0.01, 0.00, 0.03),
    ("rs2201841", 0.08, 0.09, 0.02, 0.01, 0.04),
    ("rs17234657", 0.15, 0.18, 0.05, 0.03, 0.10),
    ("rs11465804", 0.07, 0.08, 0.01, 0.00, 0.03),
    ("rs2066847", 0.04, 0.05, 0.00, 0.00, 0.01),
    ("rs10781499", 0.10, 0.12, 0.25, 0.30, 0.15),
    ("rs2395029", 0.05, 0.06, 0.02, 0.01, 0.03),
    ("rs334", 0.02, 0.00, 0.10, 0.00, 0.02),
    ("rs5030737", 0.08, 0.10, 0.02, 0.03, 0.05),
    ("rs12785878", 0.35, 0.38, 0.28, 0.25, 0.32),
    ("rs2282679", 0.27, 0.30, 0.08, 0.15, 0.22),
    ("rs12203592", 0.16, 0.18, 0.01, 0.00, 0.08),
    ("rs1800562", 0.06, 0.08, 0.00, 0.00, 0.02),
    ("rs855791", 0.45, 0.48, 0.52, 0.68, 0.50),
    ("rs4820268", 0.42, 0.45, 0.48, 0.65, 0.48),
    ("rs1805087", 0.20, 0.22, 0.12, 0.10, 0.18),
    ("rs12272669", 0.15, 0.08, 0.30, 0.25, 0.12),
    ("rs2228570", 0.35, 0.38, 0.52, 0.45, 0.40),
    ("rs1801131", 0.32, 0.35, 0.12, 0.25, 0.40),
    ("rs4654748", 0.45, 0.48, 0.38, 0.35, 0.42),
    ("rs602662", 0.45, 0.48, 0.42, 0.38, 0.40),
    ("rs3760775", 0.38, 0.40, 0.08, 0.05, 0.15),
    ("rs2187668", 0.15, 0.18, 0.02, 0.01, 0.08),
    ("rs6822844", 0.12, 0.15, 0.02, 0.03, 0.08),
    ("rs4148323", 0.10, 0.02, 0.30, 0.08, 0.05),
    ("rs6025", 0.03, 0.04, 0.01, 0.00, 0.02),
    ("rs1799963", 0.01, 0.02, 0.00, 0.00, 0.01),
    ("rs4149056", 0.15, 0.18, 0.12, 0.02, 0.10),
    ("rs1057910", 0.07, 0.08, 0.02, 0.04, 0.05),
    ("rs9923231", 0.38, 0.42, 0.08, 0.92, 0.35),
    ("rs12979860", 0.67, 0.75, 0.38, 0.95, 0.60),
    ("rs8099917", 0.20, 0.25, 0.10, 0.05, 0.15),
    ("rs3131972", 0.17, 0.20, 0.05, 0.08, 0.12),
    ("rs11240777", 0.30, 0.35, 0.22, 0.18, 0.28),
    ("rs4970383", 0.35, 0.38, 0.28, 0.25, 0.32),
    ("rs6681049", 0.40, 0.42, 0.48, 0.52, 0.45),
]


def create_database(db_path: str, drop_existing: bool = False) -> bool:
    """
    Create the GWAS database with tables and sample data.
    
    Args:
        db_path: Path to the database file.
        drop_existing: If True, drop existing tables before creating.
        
    Returns:
        bool: True if successful.
    """
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        if drop_existing:
            cursor.execute("DROP TABLE IF EXISTS gwas_fts")
            cursor.execute("DROP TABLE IF EXISTS allele_frequencies")
            cursor.execute("DROP TABLE IF EXISTS gwas_variants")
        
        # Create gwas_variants table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS gwas_variants (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                position INTEGER NOT NULL,
                ref_allele TEXT,
                alt_allele TEXT,
                risk_allele TEXT,
                reported_trait TEXT NOT NULL,
                mapped_gene TEXT,
                p_value REAL NOT NULL,
                odds_ratio REAL,
                sample_size INTEGER,
                category TEXT,
                pubmed_id TEXT,
                UNIQUE(variant_id, reported_trait)
            )
        """)
        
        # Create allele_frequencies table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS allele_frequencies (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id TEXT NOT NULL UNIQUE,
                af_overall REAL,
                af_eur REAL,
                af_afr REAL,
                af_eas REAL,
                af_amr REAL,
                FOREIGN KEY (variant_id) REFERENCES gwas_variants(variant_id)
            )
        """)
        
        # Create indexes
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_variant_id ON gwas_variants(variant_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_p_value ON gwas_variants(p_value)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_category ON gwas_variants(category)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_trait ON gwas_variants(reported_trait)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_af_variant ON allele_frequencies(variant_id)")
        
        # Create FTS5 virtual table for full-text search
        cursor.execute("DROP TABLE IF EXISTS gwas_fts")
        cursor.execute("""
            CREATE VIRTUAL TABLE gwas_fts USING fts5(
                variant_id,
                reported_trait,
                mapped_gene,
                content=gwas_variants,
                content_rowid=id
            )
        """)
        
        # Insert sample GWAS data
        cursor.executemany("""
            INSERT OR IGNORE INTO gwas_variants 
            (variant_id, chromosome, position, ref_allele, alt_allele, risk_allele,
             reported_trait, mapped_gene, p_value, odds_ratio, sample_size, category, pubmed_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, SAMPLE_GWAS_DATA)
        
        # Insert sample allele frequency data
        cursor.executemany("""
            INSERT OR IGNORE INTO allele_frequencies 
            (variant_id, af_overall, af_eur, af_afr, af_eas, af_amr)
            VALUES (?, ?, ?, ?, ?, ?)
        """, SAMPLE_AF_DATA)
        
        # Populate FTS table
        cursor.execute("""
            INSERT INTO gwas_fts(rowid, variant_id, reported_trait, mapped_gene)
            SELECT id, variant_id, reported_trait, mapped_gene FROM gwas_variants
        """)
        
        conn.commit()
        
        # Verify insertion
        cursor.execute("SELECT COUNT(*) FROM gwas_variants")
        variant_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM allele_frequencies")
        af_count = cursor.fetchone()[0]
        
        print(f"Database created successfully at {db_path}")
        print(f"  - {variant_count} variant entries")
        print(f"  - {af_count} allele frequency entries")
        
        return True
        
    except sqlite3.Error as e:
        print(f"Error creating database: {e}")
        conn.rollback()
        return False
    finally:
        conn.close()


def verify_database(db_path: str) -> bool:
    """
    Verify that the database exists and has data.
    
    Args:
        db_path: Path to the database file.
        
    Returns:
        bool: True if database is valid.
    """
    if not os.path.exists(db_path):
        return False
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM gwas_variants")
        count = cursor.fetchone()[0]
        
        conn.close()
        return count > 0
    except sqlite3.Error:
        return False


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Setup GWAS database')
    parser.add_argument('--drop', action='store_true', 
                        help='Drop existing tables before creating')
    parser.add_argument('--path', type=str, default=DATABASE_PATH,
                        help='Path to database file')
    
    args = parser.parse_args()
    
    if verify_database(args.path) and not args.drop:
        print(f"Database already exists at {args.path}")
        print("Use --drop to recreate")
    else:
        create_database(args.path, drop_existing=args.drop)
