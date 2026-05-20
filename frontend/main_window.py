"""
Main window UI for the Genetic Analysis Application.

Implements PyQt6 interface with tabbed layout for monogenic and polygenic analysis.
"""

import os
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import List, Optional
from enum import Enum
import re
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QTableWidget, QTableWidgetItem, QFileDialog, QLabel, QLineEdit,
    QComboBox, QSlider, QProgressBar, QStatusBar, QMessageBox,
    QHeaderView, QGroupBox, QSpinBox, QFrame, QSplitter, QApplication,
    QDialog, QTextEdit, QTextBrowser, QScrollArea, QDialogButtonBox, QTabWidget,
    QProgressDialog
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QUrl
from PyQt6.QtGui import QFont, QColor, QDesktopServices

from models.data_models import SNPRecord, GWASMatch, FilterCriteria
from models.polygenic_models import PolygenicResult
from backend.parsers import ParseError
from backend.parser_factory import ParserFactory
from backend.search_engine import SearchEngine, DatabaseError
from backend.scoring import get_score_interpretation
from backend.session_manager import SessionManager
from frontend.polygenic_widgets import (
    PolygenicBrowserWidget, DatabaseSettingsWidget
)
from config import (
    DATABASE_PATH, RESULTS_PER_PAGE, TRAIT_CATEGORIES,
    APP_NAME, APP_VERSION
)
from utils.logging_config import get_logger

logger = get_logger(__name__)


class ProcessingWorker(QThread):
    """
    Background worker thread for processing genetic input files.

    Signals:
        finished: Emitted when processing is complete with results (matches, stats, snp_records).
        error: Emitted when an error occurs.
        file_progress: Emitted with file loading progress (0-100, message).
        mono_progress: Emitted with monogenic analysis progress (0-100, message).
    """
    finished = pyqtSignal(list, dict, list)
    error = pyqtSignal(str)
    file_progress = pyqtSignal(int, str)
    mono_progress = pyqtSignal(int, str)
    
    def __init__(self, filepath: str, search_engine: SearchEngine) -> None:
        super().__init__()
        self.filepath = filepath
        self.search_engine = search_engine
        self._is_cancelled = False
    
    def run(self) -> None:
        """Execute the processing pipeline."""
        try:
            self.file_progress.emit(5, "Opening file...")
            self.mono_progress.emit(0, "Waiting...")
            
            parser = ParserFactory.get_parser(self.filepath)

            # Progress adapters for different parser callback signatures
            def parsing_progress_lines(lines_processed: int, total_lines: int) -> None:
                if total_lines > 0:
                    progress_percent = int((lines_processed / total_lines) * 100)
                    self.file_progress.emit(progress_percent, f"{lines_processed:,}/{total_lines:,} lines")

            def parsing_progress_percent(percent: int, message: str) -> None:
                self.file_progress.emit(percent, message)

            # Use inspect_file() when possible to obtain both records and stats.
            if hasattr(parser, 'set_progress_callback'):
                # Parser23andMe-style: accepts set_progress_callback(lines_processed, total_lines)
                parser.set_progress_callback(parsing_progress_lines)
                inspect_result = parser.inspect_file(self.filepath)
                # Parser23andMe.inspect_file returns (records, stats)
                if isinstance(inspect_result, tuple) and len(inspect_result) == 2:
                    snp_records, stats = inspect_result
                else:
                    # Fallback: try parse_file and get_parse_stats
                    snp_records = parser.parse_file(self.filepath)
                    stats = parser.get_parse_stats() if hasattr(parser, 'get_parse_stats') else {}
            else:
                # VCFStatsParser-style: inspect_file returns an object with .records and .stats
                try:
                    vcf_result = parser.inspect_file(self.filepath, progress_callback=parsing_progress_percent)
                    snp_records = vcf_result.records
                    stats = vcf_result.stats
                except TypeError:
                    # Fallback: no progress callback accepted
                    vcf_result = parser.inspect_file(self.filepath)
                    snp_records = getattr(vcf_result, 'records', [])
                    stats = getattr(vcf_result, 'stats', {})
            
            self.file_progress.emit(100, f"✓ {len(snp_records):,} SNPs loaded")
            
            if self._is_cancelled:
                return
            
            # Monogenic analysis
            self.mono_progress.emit(10, f"Matching {len(snp_records):,} SNPs...")
            
            matches = self.search_engine.match_user_snps(snp_records)
            
            if self._is_cancelled:
                return
            
            self.mono_progress.emit(100, f"✓ {len(matches):,} matches found")
            
            stats['matches_found'] = len(matches)
            self.finished.emit(matches, stats, snp_records)
            
        except ParseError as e:
            logger.error(f"Parse error: {e}")
            self.error.emit(f"Error reading file: {str(e)}")
        except DatabaseError as e:
            logger.error(f"Database error: {e}")
            self.error.emit(f"Database error: {str(e)}")
        except Exception as e:
            logger.exception(f"Unexpected error: {e}")
            self.error.emit(f"Unexpected error: {str(e)}")
    
    def cancel(self) -> None:
        """Cancel the processing operation."""
        self._is_cancelled = True


class SessionOperationWorker(QThread):
    """Background worker for saving and loading analysis sessions."""

    progress = pyqtSignal(int, str)
    finished = pyqtSignal(object)
    error = pyqtSignal(object)

    def __init__(self, operation: str, filepath: str, payload: Optional[dict] = None) -> None:
        super().__init__()
        self.operation = operation
        self.filepath = filepath
        self.payload = payload or {}

    def run(self) -> None:
        try:
            if self.operation == "save":
                SessionManager.save_session(
                    filepath=self.filepath,
                    snp_records=self.payload.get("snp_records", []),
                    gwas_matches=self.payload.get("gwas_matches", []),
                    polygenic_results=self.payload.get("polygenic_results", []),
                    metadata=self.payload.get("metadata", {}),
                    progress_callback=self._on_progress,
                )
                self.finished.emit(self.filepath)
                return

            if self.operation == "load":
                result = SessionManager.load_session(
                    self.filepath,
                    progress_callback=self._on_progress,
                )
                self.finished.emit(result)
                return

            raise ValueError(f"Unsupported session operation: {self.operation}")
        except Exception as exc:
            logger.exception("Session %s failed: %s", self.operation, exc)
            self.error.emit(exc)

    def _on_progress(self, percent: int, message: str) -> None:
        self.progress.emit(percent, message)


class HelpDialog(QDialog):
    """
    Help dialog with scrollable explanation of the application.
    """
    
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Help - Genetic Analysis Application")
        self.setMinimumSize(700, 600)
        self._init_ui()
    
    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        
        help_text = QTextBrowser()
        help_text.setReadOnly(True)
        help_text.setOpenExternalLinks(True)
        help_text.setHtml(self._get_help_content())
        help_text.setMinimumHeight(500)
        
        layout.addWidget(help_text)
        
        # Close button
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)
    
    def _get_help_content(self) -> str:
        return """
        <h1>Genetic Analysis Application</h1>
        
        <p style="text-align:center"><a href="https://github.com/pmpfe/genexplore">github.com/pmpfe/genexplore</a><br/>
        Author: <a href="mailto:pferreira@gmail.com">pferreira@gmail.com</a></p>

        <h2>What is this program?</h2>
        <p>This application analyzes your raw genetic data from 23andMe and provides two types of analysis:</p>
        <ul>
            <li><b>Monogenic Analysis:</b> Identifies individual genetic variants (SNPs) associated with traits 
            and diseases using the GWAS Catalog.</li>
            <li><b>Polygenic Analysis:</b> Calculates Polygenic Risk Scores (PRS) that combine the effects of 
            many genetic variants to predict your genetic predisposition for complex traits.</li>
        </ul>
        
        <h2>Types of Analysis</h2>
        
        <h3>🧬 Monogenic Analysis (GWAS)</h3>
        <p>This analysis looks at individual SNPs and their associations with traits. Each variant is 
        analyzed independently using data from the <b>GWAS Catalog</b> - a curated database of genome-wide 
        association studies.</p>
        <p>Results show individual variants with their statistical significance (p-value) and an 
        Impact Score (0-10) that combines significance with allele rarity.</p>
        
        <h3>📊 Polygenic Analysis (PRS)</h3>
        <p>Polygenic Risk Scores aggregate the effects of many genetic variants to estimate your 
        genetic predisposition for complex traits like height, BMI, or disease risk.</p>
        <p>Unlike monogenic analysis, PRS considers that most traits are influenced by hundreds or 
        thousands of variants, each with a small effect. The combined score provides a more 
        comprehensive picture of genetic risk.</p>
        <p><b>How PRS is calculated:</b></p>
        <ol>
            <li>For each variant in the score, your genotype is matched against the effect allele</li>
            <li>Each match contributes a weight (beta value) derived from scientific studies</li>
            <li>The weighted sum produces your raw polygenic score</li>
            <li>Scores are converted to percentiles based on population distributions</li>
        </ol>
        
        <h2>How to Use</h2>
        <ol>
            <li><b>Upload your data:</b> Click "upload gene data file" and select your raw data file</li>
            <li><b>Wait for analysis:</b> Three progress bars show:
                <ul>
                    <li>📁 File loading progress</li>
                    <li>🧬 Monogenic (GWAS) analysis progress</li>
                    <li>📊 Polygenic (PRS) analysis progress</li>
                </ul>
            </li>
            <li><b>Explore results:</b> Use the tabs to switch between Monogenic and Polygenic views</li>
            <li><b>Filter and search:</b> Use the filter controls to find specific traits or categories</li>
        </ol>
        
        <h2>Understanding Monogenic Results</h2>
        <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
            <tr style="background-color: #e0e0e0;">
                <th>Column</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><b>SNP ID</b></td>
                <td>The unique identifier (e.g., rs6983267). "rs" = Reference SNP.</td>
            </tr>
            <tr>
                <td><b>Gene</b></td>
                <td>The gene where this variant is located or nearest to.</td>
            </tr>
            <tr>
                <td><b>Trait</b></td>
                <td>The disease or characteristic associated with this variant.</td>
            </tr>
            <tr>
                <td><b>User Genotype</b></td>
                <td>Your genotype - two alleles, one from each parent (e.g., "AG").</td>
            </tr>
            <tr>
                <td><b>Risk Allele</b></td>
                <td>The allele associated with increased risk/effect. Highlighted in red if present.</td>
            </tr>
            <tr>
                <td><b>P-value</b></td>
                <td>Statistical significance. Lower = more significant. p < 5e-8 is genome-wide significant.</td>
            </tr>
            <tr>
                <td><b>Category</b></td>
                <td>Trait classification: Metabolic, Cardiovascular, Neuropsychiatric, Physical Trait, 
                Oncology, Immune, Infectious, or Other.</td>
            </tr>
            <tr>
                <td><b>Impact Score</b></td>
                <td>Score 0-10 combining p-value significance (70%) and allele rarity (30%).</td>
            </tr>
            <tr>
                <td><b>Interpretation</b></td>
                <td>Very High (≥8), High (6-8), Moderate (4-6), Low (2-4), Minimal (<2).</td>
            </tr>
        </table>
        
        <h2>Understanding Polygenic Results</h2>
        <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
            <tr style="background-color: #e0e0e0;">
                <th>Column</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><b>Trait</b></td>
                <td>The trait or condition for which the polygenic score was calculated.</td>
            </tr>
            <tr>
                <td><b>Category</b></td>
                <td>Classification of the trait type.</td>
            </tr>
            <tr>
                <td><b>Raw Score</b></td>
                <td>Your calculated polygenic score (sum of weighted allele effects).</td>
            </tr>
            <tr>
                <td><b>Percentile</b></td>
                <td>Where you fall in the population distribution (0-100%). 
                Higher percentile = higher genetic predisposition.</td>
            </tr>
            <tr>
                <td><b>Variants Used</b></td>
                <td>Number of variants from the score that matched your data.</td>
            </tr>
            <tr>
                <td><b>Risk Level</b></td>
                <td>Visual indicator: 🟢 Low (<25%), 🟡 Average (25-75%), 🔴 High (>75%).</td>
            </tr>
        </table>
        
        <h2>Using Filters</h2>
        <ul>
            <li><b>Min Impact Score:</b> Show only results above this threshold (Monogenic)</li>
            <li><b>Max P-value:</b> Filter by statistical significance (Monogenic)</li>
            <li><b>Category:</b> Filter by trait category</li>
            <li><b>Search:</b> Free text search in traits, genes, and SNP IDs</li>
        </ul>
        
        <h2>⚠️ Important Disclaimer</h2>
        <p style="color: #b00000;"><b>This application is for educational and informational purposes only.</b></p>
        <p>The results should NOT be used for medical diagnosis or treatment decisions. 
        Genetic associations are complex and influenced by many factors:</p>
        <ul>
            <li><b>Environment and lifestyle</b> often have greater effects than genetics</li>
            <li><b>Gene-gene interactions</b> are not fully captured</li>
            <li><b>Population-specific effects</b> - scores may be less accurate for non-European ancestries</li>
            <li><b>Scientific uncertainty</b> - our understanding of genetics continues to evolve</li>
        </ul>
        <p><b>Always consult with qualified healthcare professionals and genetic counselors 
        for interpretation of genetic data.</b></p>
        
        <h2>Data Sources</h2>
        <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
            <tr style="background-color: #e0e0e0;">
                <th>Source</th>
                <th>Description</th>
                <th>Link</th>
            </tr>
            <tr>
                <td><b>GWAS Catalog</b></td>
                <td>Curated database of genome-wide association studies, maintained by EMBL-EBI and NHGRI. 
                Contains thousands of variant-trait associations from published research.</td>
                <td><a href="https://www.ebi.ac.uk/gwas/">www.ebi.ac.uk/gwas/</a></td>
            </tr>
            <tr>
                <td><b>PGS Catalog</b></td>
                <td>Open database of polygenic scores and their metadata, including variant-level scoring files. 
                Contains ~4,000 scores for ~660 traits with ~150 million variant weights.</td>
                <td><a href="https://www.pgscatalog.org/">www.pgscatalog.org/</a></td>
            </tr>
            <tr>
                <td><b>dbSNP</b></td>
                <td>NCBI's database of genetic variation. Provides reference information for each SNP 
                including genomic location and population frequencies.</td>
                <td><a href="https://www.ncbi.nlm.nih.gov/snp/">www.ncbi.nlm.nih.gov/snp/</a></td>
            </tr>
            <tr>
                <td><b>gnomAD</b></td>
                <td>Genome Aggregation Database with allele frequencies from 140,000+ individuals. 
                Used to assess how rare variants are in the general population.</td>
                <td><a href="https://gnomad.broadinstitute.org/">gnomad.broadinstitute.org/</a></td>
            </tr>
            <tr>
                <td><b>Ensembl</b></td>
                <td>Genome database providing gene annotations, variant information, and cross-references 
                to other databases.</td>
                <td><a href="https://www.ensembl.org/">www.ensembl.org/</a></td>
            </tr>
            <tr>
                <td><b>ClinVar</b></td>
                <td>NCBI database of clinically relevant genetic variants and their relationships 
                to human health.</td>
                <td><a href="https://www.ncbi.nlm.nih.gov/clinvar/">www.ncbi.nlm.nih.gov/clinvar/</a></td>
            </tr>
        </table>
        
        <h2>Further Reading</h2>
        <ul>
            <li><b>Understanding GWAS:</b> <a href="https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies">genome.gov - GWAS Glossary</a></li>
            <li><b>Polygenic Scores Explained:</b> <a href="https://www.pgscatalog.org/about/">PGS Catalog - About</a></li>
            <li><b>Genetics Education:</b> <a href="https://www.genome.gov/For-Patients-and-Families">NIH Genetics Education Resources</a></li>
            <li><b>Scientific Papers:</b>
                <ul>
                    <li><a href="https://doi.org/10.1093/nar/gkaa1061">GWAS Catalog 2023 Update (Nucleic Acids Research)</a></li>
                    <li><a href="https://doi.org/10.1038/s41588-021-00783-5">PGS Catalog Paper (Nature Genetics)</a></li>
                </ul>
            </li>
        </ul>
        
        <h2>Database Updates</h2>
        <p>To update the local databases with the latest data from GWAS Catalog and PGS Catalog, 
        run the update script from the command line:</p>
        <pre style="background-color: #f0f0f0; padding: 10px; border-radius: 5px;">
python update_databases.py --all     # Update both databases
python update_databases.py --gwas    # Update only GWAS
python update_databases.py --pgs     # Update only PGS Catalog
        </pre>
        <p>This process downloads the latest data and may take several hours for the full PGS Catalog 
        (~150 million variants).</p>
        
        <h2>Technical Information</h2>
        <ul>
            <li><b>Application Version:</b> 1.0.0</li>
            <li><b>Supported Input Formats:</b> 23andMe raw data (.txt)</li>
            <li><b>Database Format:</b> SQLite with optimized indexes</li>
            <li><b>Source Code:</b> <a href="https://github.com/">Available on GitHub</a></li>
        </ul>
        """


class ExplainDialog(QDialog):
    """
    Dialog showing detailed explanation for a specific GWAS match.
    """
    
    def __init__(self, match: GWASMatch, parent=None) -> None:
        super().__init__(parent)
        self.match = match
        self.setWindowTitle(f"Details: {match.rsid} - {match.trait[:50]}")
        self.setMinimumSize(650, 550)
        self._init_ui()
    
    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        
        # Scrollable content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        content = QWidget()
        content_layout = QVBoxLayout(content)
        
        explain_text = QTextEdit()
        explain_text.setReadOnly(True)
        explain_text.setHtml(self._get_explanation())
        explain_text.setMinimumHeight(400)
        
        content_layout.addWidget(explain_text)
        scroll.setWidget(content)
        layout.addWidget(scroll)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        # External links
        dbsnp_btn = QPushButton("🔗 Open in dbSNP")
        dbsnp_btn.clicked.connect(self._open_dbsnp)
        button_layout.addWidget(dbsnp_btn)
        
        gwas_btn = QPushButton("🔗 Open in GWAS Catalog")
        gwas_btn.clicked.connect(self._open_gwas_catalog)
        button_layout.addWidget(gwas_btn)
        
        if self.match.gene:
            gene_btn = QPushButton(f"🔗 Gene: {self.match.gene}")
            gene_btn.clicked.connect(self._open_gene_info)
            button_layout.addWidget(gene_btn)
        
        button_layout.addStretch()
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        button_layout.addWidget(close_btn)
        
        layout.addLayout(button_layout)
    
    def _get_explanation(self) -> str:
        m = self.match
        
        # Calculate score components for explanation
        import math
        neg_log_p = -math.log10(m.p_value) if m.p_value > 0 else 0
        p_score = min(neg_log_p / 10, 1.0) * 7.0
        af_score = (1 - m.allele_frequency) * 3.0
        
        # Risk allele analysis
        has_risk = m.risk_allele in m.user_genotype
        risk_count = m.user_genotype.count(m.risk_allele)
        if risk_count == 2:
            risk_status = f"<span style='color: #b00000;'><b>Homozygous for risk allele</b> (2 copies)</span>"
        elif risk_count == 1:
            risk_status = f"<span style='color: #b06000;'><b>Heterozygous</b> (1 copy of risk allele)</span>"
        else:
            risk_status = "<span style='color: #006000;'><b>No risk allele</b> (0 copies)</span>"
        
        # Format odds ratio
        odds_ratio_str = f"{m.odds_ratio:.2f}" if m.odds_ratio else "Not available"
        
        # Format sample size
        sample_size_str = f"{m.sample_size:,} individuals" if m.sample_size else "Not available"
        
        return f"""
        <h2>Variant Information</h2>
        <table border="0" cellpadding="5" style="width: 100%;">
            <tr><td><b>SNP ID:</b></td><td>{m.rsid}</td></tr>
            <tr><td><b>Chromosome:</b></td><td>{m.chromosome}</td></tr>
            <tr><td><b>Position:</b></td><td>{m.position:,}</td></tr>
            <tr><td><b>Gene:</b></td><td>{m.gene or 'Intergenic (not in a gene)'}</td></tr>
        </table>
        
        <h2>Association Details</h2>
        <table border="0" cellpadding="5" style="width: 100%;">
            <tr><td><b>Associated Trait:</b></td><td>{m.trait}</td></tr>
            <tr><td><b>Category:</b></td><td>{m.category}</td></tr>
            <tr><td><b>Risk Allele:</b></td><td>{m.risk_allele}</td></tr>
            <tr><td><b>Odds Ratio:</b></td><td>{odds_ratio_str}</td></tr>
            <tr><td><b>Study Sample Size:</b></td><td>{sample_size_str}</td></tr>
        </table>
        
        <h2>Your Genotype Analysis</h2>
        <table border="0" cellpadding="5" style="width: 100%;">
            <tr><td><b>Your Genotype:</b></td><td><b style="font-size: 14pt;">{m.user_genotype}</b></td></tr>
            <tr><td><b>Risk Allele Status:</b></td><td>{risk_status}</td></tr>
        </table>
        
        <h2>Statistical Significance</h2>
        <table border="0" cellpadding="5" style="width: 100%;">
            <tr><td><b>P-value:</b></td><td>{m.p_value:.2e}</td></tr>
            <tr><td><b>-log₁₀(p-value):</b></td><td>{neg_log_p:.2f}</td></tr>
            <tr><td><b>Interpretation:</b></td><td>
                {'Highly significant (p < 5×10⁻⁸)' if m.p_value < 5e-8 else 
                 'Significant (p < 0.001)' if m.p_value < 0.001 else
                 'Suggestive (p < 0.05)' if m.p_value < 0.05 else 'Not significant'}</td></tr>
        </table>
        
        <h2>Population Frequency</h2>
        <table border="0" cellpadding="5" style="width: 100%;">
            <tr><td><b>Allele Frequency:</b></td><td>{m.allele_frequency:.1%}</td></tr>
            <tr><td><b>Interpretation:</b></td><td>
                {'Very rare variant (<1%)' if m.allele_frequency < 0.01 else
                 'Rare variant (1-5%)' if m.allele_frequency < 0.05 else
                 'Low frequency (5-10%)' if m.allele_frequency < 0.10 else
                 'Common variant (≥10%)'}</td></tr>
        </table>
        
        <h2>Impact Score Calculation</h2>
        <p>The Impact Score combines statistical significance with variant rarity:</p>
        <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
            <tr style="background-color: #e0e0e0;">
                <th>Component</th>
                <th>Formula</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>P-value Score</td>
                <td>min(-log₁₀(p) / 10, 1) × 7</td>
                <td>{p_score:.2f}</td>
            </tr>
            <tr>
                <td>Rarity Score</td>
                <td>(1 - allele_frequency) × 3</td>
                <td>{af_score:.2f}</td>
            </tr>
            <tr style="background-color: #ffffcc;">
                <td><b>Total Impact Score</b></td>
                <td>P-value + Rarity (clamped 0-10)</td>
                <td><b>{m.impact_score:.2f}</b></td>
            </tr>
        </table>
        
        <h2>What This Means</h2>
        <p>This variant ({m.rsid}) has been associated with <b>{m.trait}</b> in genome-wide 
        association studies{f' involving {m.sample_size:,} participants' if m.sample_size else ''}.</p>
        
        {'<p style="color: #b00000;">⚠️ You carry ' + str(risk_count) + ' cop' + ('y' if risk_count == 1 else 'ies') + 
         ' of the risk allele. This may indicate a slightly increased statistical risk, but genetics is only one factor among many.</p>' 
         if has_risk else 
         '<p style="color: #006000;">✓ You do not carry the risk allele for this association.</p>'}
        
        <p><b>Remember:</b> Statistical association does not mean causation. Many genetic 
        and environmental factors influence traits and disease risk. Consult healthcare 
        professionals for personalized interpretation.</p>
        """
    
    def _open_dbsnp(self) -> None:
        url = f"https://www.ncbi.nlm.nih.gov/snp/{self.match.rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _open_gwas_catalog(self) -> None:
        url = f"https://www.ebi.ac.uk/gwas/variants/{self.match.rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _open_gene_info(self) -> None:
        if self.match.gene:
            url = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={self.match.gene}"
            QDesktopServices.openUrl(QUrl(url))


class MainWindow(QMainWindow):
    """
    Main application window for the Genetic Analysis Application.
    
    Provides:
    - Tabbed interface for monogenic and polygenic analysis
    - File upload for 23andMe data (shared between tabs)
    - Results table with sorting and pagination
    - Real-time filtering by score, p-value, category, and text search
    """
    
    def __init__(self) -> None:
        super().__init__()
        
        self.all_matches: List[GWASMatch] = []
        self.current_matches: List[GWASMatch] = []
        self.filtered_matches: List[GWASMatch] = []
        self.current_page = 0
        self.worker: Optional[ProcessingWorker] = None
        self.snp_records: List[SNPRecord] = []
        self.polygenic_results: List[PolygenicResult] = []
        self.current_session_metadata: dict = {}
        self.current_inspection = None
        self.analysis_started_at = None
        self.analysis_completed_at = None
        self._last_loaded_file = None
        self.session_worker: Optional[SessionOperationWorker] = None
        self.session_progress_dialog = None
        self._current_session_operation = None
        self._current_session_filepath = None
        
        self.search_engine = SearchEngine(DATABASE_PATH)
        
        self._init_ui()
        self._setup_connections()
        self._verify_database()
    
    def _init_ui(self) -> None:
        """Initialize the user interface with tabbed layout."""
        self.setWindowTitle(f"{APP_NAME} v{APP_VERSION}")
        self.setMinimumSize(1200, 800)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(10)
        main_layout.setContentsMargins(15, 15, 15, 15)
        
        # Header section with upload controls (shared across tabs)
        header_layout = QHBoxLayout()
        
        title_label = QLabel(APP_NAME)
        title_label.setFont(QFont('Arial', 18, QFont.Weight.Bold))
        header_layout.addWidget(title_label)
        
        header_layout.addStretch()
        
        # Genotype status indicator
        self.genotype_status = QLabel("No genetic data loaded")
        self.genotype_status.setStyleSheet("color: #666; font-style: italic;")
        header_layout.addWidget(self.genotype_status)
        
        header_layout.addWidget(QLabel(" | "))
        
        self.help_btn = QPushButton("❓ Help")
        self.help_btn.setMinimumWidth(80)
        self.help_btn.setMinimumHeight(40)
        self.help_btn.setFont(QFont('Arial', 11))
        header_layout.addWidget(self.help_btn)
        
        # Save/Load session buttons
        self.save_btn = QPushButton("💾 Save")
        self.save_btn.setMinimumWidth(80)
        self.save_btn.setMinimumHeight(40)
        self.save_btn.setFont(QFont('Arial', 11))
        self.save_btn.setEnabled(False)
        self.save_btn.setToolTip("Save analysis results to file")
        header_layout.addWidget(self.save_btn)
        
        self.export_btn = QPushButton("🌐 Export HTML")
        self.export_btn.setMinimumWidth(120)
        self.export_btn.setMinimumHeight(40)
        self.export_btn.setFont(QFont('Arial', 11))
        self.export_btn.setEnabled(False)
        self.export_btn.setToolTip("Export analysis report as a self-contained HTML file")
        header_layout.addWidget(self.export_btn)
        
        self.load_btn = QPushButton("📂 Load")
        self.load_btn.setMinimumWidth(80)
        self.load_btn.setMinimumHeight(40)
        self.load_btn.setFont(QFont('Arial', 11))
        self.load_btn.setToolTip("Load previously saved analysis")
        header_layout.addWidget(self.load_btn)
        
        self.upload_btn = QPushButton("📁 upload gene data file")
        self.upload_btn.setMinimumWidth(200)
        self.upload_btn.setMinimumHeight(40)
        self.upload_btn.setFont(QFont('Arial', 11))
        header_layout.addWidget(self.upload_btn)

        # Connect export button
        self.export_btn.clicked.connect(self.export_report)
        
        main_layout.addLayout(header_layout)
        
        # Progress section with 3 separate progress bars
        self.progress_frame = QFrame()
        self.progress_frame.setVisible(False)
        progress_main_layout = QVBoxLayout(self.progress_frame)
        progress_main_layout.setContentsMargins(0, 5, 0, 5)
        progress_main_layout.setSpacing(5)
        
        # Status label
        self.progress_label = QLabel("Processing...")
        progress_main_layout.addWidget(self.progress_label)
        
        # Three progress bars in a horizontal layout
        progress_bars_layout = QHBoxLayout()
        progress_bars_layout.setSpacing(15)
        
        # 1. File parsing progress
        file_progress_group = QVBoxLayout()
        file_progress_group.setSpacing(2)
        self.file_progress_label = QLabel("📁 File Loading")
        self.file_progress_label.setFont(QFont('Arial', 9))
        file_progress_group.addWidget(self.file_progress_label)
        self.file_progress_bar = QProgressBar()
        self.file_progress_bar.setRange(0, 100)
        self.file_progress_bar.setMinimumWidth(200)
        self.file_progress_bar.setMaximumHeight(20)
        file_progress_group.addWidget(self.file_progress_bar)
        progress_bars_layout.addLayout(file_progress_group)
        
        # 2. Monogenic analysis progress
        mono_progress_group = QVBoxLayout()
        mono_progress_group.setSpacing(2)
        self.mono_progress_label = QLabel("🧬 Monogenic Analysis")
        self.mono_progress_label.setFont(QFont('Arial', 9))
        mono_progress_group.addWidget(self.mono_progress_label)
        self.mono_progress_bar = QProgressBar()
        self.mono_progress_bar.setRange(0, 100)
        self.mono_progress_bar.setMinimumWidth(200)
        self.mono_progress_bar.setMaximumHeight(20)
        mono_progress_group.addWidget(self.mono_progress_bar)
        progress_bars_layout.addLayout(mono_progress_group)
        
        # 3. Polygenic analysis progress
        poly_progress_group = QVBoxLayout()
        poly_progress_group.setSpacing(2)
        self.poly_progress_label = QLabel("📊 Polygenic Analysis")
        self.poly_progress_label.setFont(QFont('Arial', 9))
        poly_progress_group.addWidget(self.poly_progress_label)
        self.poly_progress_bar = QProgressBar()
        self.poly_progress_bar.setRange(0, 100)
        self.poly_progress_bar.setMinimumWidth(200)
        self.poly_progress_bar.setMaximumHeight(20)
        poly_progress_group.addWidget(self.poly_progress_bar)
        progress_bars_layout.addLayout(poly_progress_group)
        
        progress_bars_layout.addStretch()
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setMaximumWidth(80)
        self.cancel_btn.setMaximumHeight(35)
        progress_bars_layout.addWidget(self.cancel_btn)
        
        progress_main_layout.addLayout(progress_bars_layout)
        
        main_layout.addWidget(self.progress_frame)
        
        # Tab widget
        self.tab_widget = QTabWidget()
        self.tab_widget.setFont(QFont('Arial', 11))
        
        # Tab 1: Monogenic Analysis (existing functionality)
        self.monogenic_tab = QWidget()
        self._init_monogenic_tab()
        self.tab_widget.addTab(self.monogenic_tab, "🧬 Monogenic Analysis")
        
        # Tab 2: Polygenic Analysis (new)
        self.polygenic_widget = PolygenicBrowserWidget()
        self.tab_widget.addTab(self.polygenic_widget, "📊 Polygenic Scores")
        
        # Tab 3: Database Settings
        self.settings_widget = DatabaseSettingsWidget()
        self.tab_widget.addTab(self.settings_widget, "⚙️ Database Settings")
        
        main_layout.addWidget(self.tab_widget, stretch=1)
        
        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready - upload gene data file to begin")
        
        self._apply_styles()
    
    def _init_monogenic_tab(self) -> None:
        """Initialize the monogenic analysis tab content."""
        layout = QVBoxLayout(self.monogenic_tab)
        layout.setContentsMargins(0, 10, 0, 0)
        
        # Filters section
        self.filters_group = QGroupBox("Filters & Search")
        filters_layout = QVBoxLayout(self.filters_group)
        
        # First row of filters
        filter_row1 = QHBoxLayout()
        
        # Score filter
        score_group = QVBoxLayout()
        score_label = QLabel("Min Impact Score:")
        self.score_slider = QSlider(Qt.Orientation.Horizontal)
        self.score_slider.setRange(0, 100)
        self.score_slider.setValue(0)
        self.score_value_label = QLabel("0.0")
        score_group.addWidget(score_label)
        score_group.addWidget(self.score_slider)
        score_group.addWidget(self.score_value_label)
        filter_row1.addLayout(score_group)
        
        # P-value filter
        pvalue_group = QVBoxLayout()
        pvalue_label = QLabel("Max P-value (log scale):")
        self.pvalue_slider = QSlider(Qt.Orientation.Horizontal)
        self.pvalue_slider.setRange(0, 100)
        self.pvalue_slider.setValue(73) # Corresponds to p-value of 5e-8
        self.pvalue_value_label = QLabel("1.0")
        pvalue_group.addWidget(pvalue_label)
        pvalue_group.addWidget(self.pvalue_slider)
        pvalue_group.addWidget(self.pvalue_value_label)
        filter_row1.addLayout(pvalue_group)
        
        # Category filter
        category_group = QVBoxLayout()
        category_label = QLabel("Category:")
        self.category_combo = QComboBox()
        self.category_combo.addItems(TRAIT_CATEGORIES)
        category_group.addWidget(category_label)
        category_group.addWidget(self.category_combo)
        filter_row1.addLayout(category_group)
        
        # Carrier status filter
        carrier_group = QVBoxLayout()
        carrier_label = QLabel("Carrier Status:")
        self.carrier_combo = QComboBox()
        self.carrier_combo.addItems(['ALL', 'non-carrier', 'heterozygous', 'homozygous', 'carrier'])
        carrier_group.addWidget(carrier_label)
        carrier_group.addWidget(self.carrier_combo)
        filter_row1.addLayout(carrier_group)
        
        filters_layout.addLayout(filter_row1)
        
        # Second row - search
        filter_row2 = QHBoxLayout()
        search_label = QLabel("Search (traits, genes, SNPs):")
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Type to search... e.g., diabetes, MTHFR, rs7903146")
        self.search_input.setClearButtonEnabled(True)
        filter_row2.addWidget(search_label)
        filter_row2.addWidget(self.search_input, stretch=1)
        
        self.reset_filters_btn = QPushButton("Reset Filters")
        filter_row2.addWidget(self.reset_filters_btn)
        
        filters_layout.addLayout(filter_row2)
        
        layout.addWidget(self.filters_group)
        self.filters_group.setEnabled(False)
        
        # Results count label
        self.results_label = QLabel("No results to display")
        self.results_label.setFont(QFont('Arial', 11))
        layout.addWidget(self.results_label)
        
        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(10)
        self.results_table.setHorizontalHeaderLabels([
            "SNP ID", "Gene", "Trait", "User Genotype", 
            "Risk Allele", "P-value", "Category", "Impact Score", "Interpretation", ""
        ])
        
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(6, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(7, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(8, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(9, QHeaderView.ResizeMode.ResizeToContents)
        
        self.results_table.setSortingEnabled(True)
        self.results_table.setAlternatingRowColors(True)
        
        layout.addWidget(self.results_table, stretch=1)
        
        # Pagination section
        pagination_layout = QHBoxLayout()
        
        self.prev_btn = QPushButton("◀ Previous")
        self.prev_btn.setEnabled(False)
        pagination_layout.addWidget(self.prev_btn)
        
        pagination_layout.addStretch()
        
        self.page_label = QLabel("Page 0 of 0")
        pagination_layout.addWidget(self.page_label)
        
        pagination_layout.addStretch()
        
        self.next_btn = QPushButton("Next ▶")
        self.next_btn.setEnabled(False)
        pagination_layout.addWidget(self.next_btn)
        
        layout.addLayout(pagination_layout)
    
    def _apply_styles(self) -> None:
        """Apply stylesheet to the application."""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f5;
            }
            QTabWidget::pane {
                border: 1px solid #ccc;
                border-radius: 4px;
                background-color: white;
            }
            QTabBar::tab {
                background-color: #e0e0e0;
                border: 1px solid #ccc;
                border-bottom: none;
                padding: 10px 20px;
                margin-right: 2px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
            }
            QTabBar::tab:selected {
                background-color: white;
                border-bottom: 1px solid white;
            }
            QTabBar::tab:hover:!selected {
                background-color: #d0d0d0;
            }
            QGroupBox {
                font-weight: bold;
                border: 1px solid #ccc;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
            QPushButton {
                background-color: #4a90d9;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #357abd;
            }
            QPushButton:disabled {
                background-color: #cccccc;
            }
            QTableWidget {
                gridline-color: #ddd;
                background-color: white;
                alternate-background-color: #f9f9f9;
            }
            QTableWidget::item {
                padding: 5px;
            }
            QHeaderView::section {
                background-color: #e0e0e0;
                padding: 8px;
                border: 1px solid #ccc;
                font-weight: bold;
            }
            QLineEdit {
                padding: 8px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }
            QSlider::groove:horizontal {
                height: 8px;
                background: #ddd;
                border-radius: 4px;
            }
            QSlider::handle:horizontal {
                background: #4a90d9;
                width: 16px;
                margin: -4px 0;
                border-radius: 8px;
            }
            QProgressBar {
                border: 1px solid #ccc;
                border-radius: 4px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #4a90d9;
            }
        """)
    
    def _setup_connections(self) -> None:
        """Set up signal-slot connections."""
        self.upload_btn.clicked.connect(self._on_upload_clicked)
        self.help_btn.clicked.connect(self._on_help_clicked)
        self.save_btn.clicked.connect(self._on_save_clicked)
        self.load_btn.clicked.connect(self._on_load_clicked)
        self.cancel_btn.clicked.connect(self._on_cancel_clicked)
        self.prev_btn.clicked.connect(self._on_prev_page)
        self.next_btn.clicked.connect(self._on_next_page)
        self.reset_filters_btn.clicked.connect(self._reset_filters)
        
        # Real-time filter connections
        self.score_slider.valueChanged.connect(self._on_score_changed)
        self.pvalue_slider.valueChanged.connect(self._on_pvalue_changed)
        self.category_combo.currentTextChanged.connect(self._apply_filters)
        self.carrier_combo.currentTextChanged.connect(self._apply_filters)
        self.search_input.textChanged.connect(self._on_search_changed)
        
        # Debounce timer for search
        self._search_timer = QTimer()
        self._search_timer.setSingleShot(True)
        self._search_timer.timeout.connect(self._apply_filters)
    
    def _verify_database(self) -> None:
        """Verify the database exists and is valid."""
        if not self.search_engine.verify_database():
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setWindowTitle("Database Not Found")
            msg.setText("GWAS database not found.")
            msg.setInformativeText(
                "Please run 'python database/setup_database.py' to create the database."
            )
            msg.exec()
            self.status_bar.showMessage("Warning: Database not available")
            self.upload_btn.setEnabled(False)
    
    def _on_upload_clicked(self) -> None:
        """Handle file upload button click."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Select genetic data file (23andMe or VCF)",
            "",
            "All Files (*);;23andMe raw text (*.txt);;VCF files (*.vcf *.vcf.gz)",
            "All Files (*)"
        )
        
        if not filepath:
            return
        
        self._start_processing(filepath)
    
    def _start_processing(self, filepath: str) -> None:
        """Start processing a genetic data file in background."""
        self._last_loaded_file = filepath
        self.analysis_started_at = datetime.now()
        try:
            parser = ParserFactory.get_parser(filepath)
            detected_format = parser.format_name
        except Exception:
            detected_format = 'Unknown'
        self.current_inspection = SimpleNamespace(
            filepath=filepath,
            filename=Path(filepath).name,
            detected_format=detected_format,
            file_size_bytes=os.path.getsize(filepath) if os.path.exists(filepath) else None,
        )
        self.current_session_metadata = {
            'snp_file': filepath,
            'source_filename': Path(filepath).name,
            'source_format': detected_format,
            'analysis_started_at': self.analysis_started_at.isoformat(),
        }
        self.upload_btn.setEnabled(False)
        self.filters_group.setEnabled(False)
        self.progress_frame.setVisible(True)
        
        # Initialize all progress bars
        self.file_progress_bar.setValue(0)
        self.mono_progress_bar.setValue(0)
        self.poly_progress_bar.setValue(0)
        self.file_progress_label.setText("📁 File Loading")
        self.mono_progress_label.setText("🧬 Monogenic Analysis")
        self.poly_progress_label.setText("📊 Polygenic Analysis")
        self.progress_label.setText("Loading genetic data...")
        
        self.worker = ProcessingWorker(filepath, self.search_engine)
        self.worker.file_progress.connect(self._on_file_progress)
        self.worker.mono_progress.connect(self._on_mono_progress)
        self.worker.finished.connect(self._on_processing_finished)
        self.worker.error.connect(self._on_processing_error)
        self.worker.start()
    
    def _on_file_progress(self, value: int, message: str) -> None:
        """Handle file loading progress updates."""
        self.file_progress_bar.setValue(value)
        self.file_progress_label.setText(f"📁 {message}")
        if value < 100:
            self.progress_label.setText("Loading genetic data...")
            self.status_bar.showMessage(f"File loading: {message}")
    
    def _on_mono_progress(self, value: int, message: str) -> None:
        """Handle monogenic analysis progress updates."""
        self.mono_progress_bar.setValue(value)
        self.mono_progress_label.setText(f"🧬 {message}")
        if value > 0 and value < 100:
            self.progress_label.setText("Analyzing monogenic variants...")
            self.status_bar.showMessage(f"Monogenic analysis: {message}")
        elif value == 100:
            self.progress_label.setText("Monogenic analysis complete!")
    
    def _on_poly_progress(self, current: int, total: int, message: str) -> None:
        """Handle polygenic analysis progress updates."""
        if total > 0:
            percent = int(current / total * 100)
            self.poly_progress_bar.setValue(percent)
            self.poly_progress_label.setText(f"📊 {message}")
            self.status_bar.showMessage(f"Polygenic: {current}/{total} scores")
    
    def _on_processing_finished(self, matches: List[GWASMatch], stats: dict, snp_records: List[SNPRecord]) -> None:
        """Handle successful processing completion."""
        # Keep progress frame visible but enable interface
        self.upload_btn.setEnabled(True)
        self.filters_group.setEnabled(True)
        
        self.all_matches = matches
        self.current_matches = matches
        self.snp_records = snp_records
        self.current_page = 0
        self.analysis_completed_at = datetime.now()
        self.current_session_metadata.update({
            'analysis_completed_at': self.analysis_completed_at.isoformat(),
            'source_record_count': len(snp_records),
            'gwas_match_count': len(matches),
            'polygenic_score_count': len(self.polygenic_results),
        })
        
        logger.info(f"Processing complete: {stats}")
        
        self._reset_filters()
        # Allow exporting monogenic results immediately
        try:
            self.export_btn.setEnabled(True)
        except Exception:
            pass
        
        # Update genotype status indicator
        self.genotype_status.setText(f"✓ {len(snp_records):,} SNPs loaded")
        self.genotype_status.setStyleSheet("color: #006000; font-weight: bold;")
        
        # Share genotype data with polygenic analysis tab and start background computation
        self.polygenic_widget.set_genotype_data(snp_records)
        
        # Start polygenic computation in background (non-blocking)
        self._start_polygenic_computation()
        
        self.status_bar.showMessage(
            f"Loaded {stats.get('valid_snps', 0):,} SNPs | "
            f"Found {stats.get('matches_found', 0):,} GWAS matches | "
            f"Polygenic computing in background..."
        )
    
    def _start_polygenic_computation(self) -> None:
        """Start polygenic score computation in background."""
        if not self.snp_records:
            return
        
        # Update progress bar
        self.poly_progress_bar.setValue(0)
        self.poly_progress_label.setText("📊 Starting computation...")
        self.progress_label.setText("Computing polygenic scores in background...")
        
        # Connect to polygenic widget's worker signals
        self.polygenic_widget.compute_btn.setEnabled(False)
        
        # Start computation via the polygenic widget
        # We hook into its signals for progress
        if hasattr(self.polygenic_widget, 'worker') and self.polygenic_widget.worker:
            self.polygenic_widget.worker.progress.connect(self._on_poly_progress)
            self.polygenic_widget.worker.finished.connect(self._on_polygenic_finished)
        
        # Trigger computation
        self.polygenic_widget._start_computation()
        
        # Re-connect signals after worker is created
        if self.polygenic_widget.worker:
            try:
                self.polygenic_widget.worker.progress.disconnect(self._on_poly_progress)
            except:
                pass
            self.polygenic_widget.worker.progress.connect(self._on_poly_progress)
            try:
                self.polygenic_widget.worker.finished.disconnect(self._on_polygenic_finished)
            except:
                pass
            self.polygenic_widget.worker.finished.connect(self._on_polygenic_finished)
    
    def _on_polygenic_finished(self, results) -> None:
        """Handle polygenic computation completion."""
        self.polygenic_results = results
        self.poly_progress_bar.setValue(100)
        self.poly_progress_label.setText(f"📊 ✓ {len(results)} scores computed")
        self.progress_label.setText("All analyses complete!")
        
        # Enable save button now that we have complete results
        self.save_btn.setEnabled(True)
        # Enable export button as well
        try:
            self.export_btn.setEnabled(True)
        except Exception:
            pass
        self.analysis_completed_at = datetime.now()
        self.current_session_metadata.update({
            'analysis_completed_at': self.analysis_completed_at.isoformat(),
            'polygenic_score_count': len(results),
        })
        
        # Update status bar with final summary
        self.status_bar.showMessage(
            f"✓ {len(self.snp_records):,} SNPs | "
            f"{len(self.all_matches):,} GWAS matches | "
            f"{len(results)} polygenic scores"
        )
        
        # Hide progress frame after a short delay
        QTimer.singleShot(2000, self._hide_progress_if_complete)
    
    def _hide_progress_if_complete(self) -> None:
        """Hide progress frame if all computations are complete."""
        # Check if all progress bars are at 100%
        if (self.file_progress_bar.value() == 100 and 
            self.mono_progress_bar.value() == 100 and
            self.poly_progress_bar.value() == 100):
            self.progress_frame.setVisible(False)
    
    def _on_processing_error(self, error_message: str) -> None:
        """Handle processing error."""
        self.progress_frame.setVisible(False)
        self.upload_btn.setEnabled(True)
        
        QMessageBox.critical(
            self,
            "Processing Error",
            error_message
        )
        
        self.status_bar.showMessage(f"Error: {error_message}")
    
    def _on_cancel_clicked(self) -> None:
        """Handle cancel button click."""
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.worker.wait()
        
        # Also cancel polygenic worker if running
        if (hasattr(self.polygenic_widget, 'worker') and 
            self.polygenic_widget.worker and 
            self.polygenic_widget.worker.isRunning()):
            self.polygenic_widget._cancel_computation()
        
        self.progress_frame.setVisible(False)
        self.upload_btn.setEnabled(True)
        self.status_bar.showMessage("Processing cancelled")
    
    def _on_score_changed(self, value: int) -> None:
        """Handle score slider change."""
        score = value / 10.0
        self.score_value_label.setText(f"{score:.1f}")
        self._apply_filters()
    
    def _on_pvalue_changed(self, value: int) -> None:
        """Handle p-value slider change."""
        if value == 100:
            pvalue = 1.0
        else:
            exponent = -value / 10.0
            pvalue = 10 ** exponent
        self.pvalue_value_label.setText(f"{pvalue:.2e}")
        self._apply_filters()
    
    def _on_search_changed(self) -> None:
        """Handle search input change with debounce."""
        self._search_timer.start(300)
    
    def _get_current_filter_criteria(self) -> FilterCriteria:
        """Get current filter settings as FilterCriteria object."""
        min_score = self.score_slider.value() / 10.0
        
        pvalue_slider = self.pvalue_slider.value()
        if pvalue_slider == 100:
            max_pvalue = 1.0
        else:
            max_pvalue = 10 ** (-pvalue_slider / 10.0)
        
        return FilterCriteria(
            min_score=min_score,
            max_pvalue=max_pvalue,
            category=self.category_combo.currentText(),
            carrier_status=self.carrier_combo.currentText(),
            search_text=self.search_input.text().strip(),
            sort_by='score',
            sort_ascending=False
        )
    
    def _apply_filters(self) -> None:
        """Apply current filter settings and refresh table."""
        if not self.all_matches:
            return
        
        criteria = self._get_current_filter_criteria()
        self.filtered_matches = criteria.apply_to_matches(self.all_matches)
        self.current_page = 0
        self._update_table()
    
    def _reset_filters(self) -> None:
        """Reset all filters to default values."""
        self.score_slider.setValue(0)
        self.pvalue_slider.setValue(100)
        self.category_combo.setCurrentIndex(0)
        self.carrier_combo.setCurrentIndex(0)
        self.search_input.clear()
        
        self.filtered_matches = list(self.all_matches)
        self.filtered_matches.sort(key=lambda m: m.impact_score, reverse=True)
        self.current_page = 0
        self._update_table()
    
    def _update_table(self) -> None:
        """Update the results table with current filtered data."""
        total_results = len(self.filtered_matches)
        total_pages = max(1, (total_results + RESULTS_PER_PAGE - 1) // RESULTS_PER_PAGE)
        
        start_idx = self.current_page * RESULTS_PER_PAGE
        end_idx = min(start_idx + RESULTS_PER_PAGE, total_results)
        
        page_matches = self.filtered_matches[start_idx:end_idx]
        
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(page_matches))
        
        for row, match in enumerate(page_matches):
            self._set_table_row(row, match)
        
        self.results_table.setSortingEnabled(True)
        
        # Update labels and pagination
        if total_results > 0:
            self.results_label.setText(
                f"Showing {start_idx + 1}-{end_idx} of {total_results} results "
                f"(from {len(self.all_matches)} total matches)"
            )
        else:
            self.results_label.setText("No matching results found")
        
        self.page_label.setText(f"Page {self.current_page + 1} of {total_pages}")
        self.prev_btn.setEnabled(self.current_page > 0)
        self.next_btn.setEnabled(self.current_page < total_pages - 1)
    
    def _set_table_row(self, row: int, match: GWASMatch) -> None:
        """Set data for a single table row."""
        items = [
            match.rsid,
            match.gene or "-",
            match.trait,
            match.user_genotype,
            match.risk_allele,
            f"{match.p_value:.2e}",
            match.category,
            f"{match.impact_score:.2f}",
            get_score_interpretation(match.impact_score)
        ]
        
        for col, value in enumerate(items):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            
            # Color coding for impact score
            if col == 7:
                if match.impact_score >= 8:
                    item.setBackground(QColor(255, 200, 200))
                elif match.impact_score >= 6:
                    item.setBackground(QColor(255, 230, 200))
                elif match.impact_score >= 4:
                    item.setBackground(QColor(255, 255, 200))
            
            # Highlight if user has risk allele
            if col == 3 and match.has_risk_allele():
                item.setForeground(QColor(180, 0, 0))
                item.setFont(QFont('Arial', 10, QFont.Weight.Bold))
            
            self.results_table.setItem(row, col, item)
        
        # Add Explain button in last column
        explain_btn = QPushButton("🔍 Explain")
        explain_btn.setStyleSheet("""
            QPushButton {
                background-color: #5cb85c;
                color: white;
                border: none;
                padding: 4px 8px;
                border-radius: 3px;
                font-size: 11px;
            }
            QPushButton:hover {
                background-color: #449d44;
            }
        """)
        explain_btn.clicked.connect(lambda checked, m=match: self._show_explain_dialog(m))
        self.results_table.setCellWidget(row, 9, explain_btn)
    
    def _show_explain_dialog(self, match: GWASMatch) -> None:
        """Show the explain dialog for a specific match."""
        dialog = ExplainDialog(match, self)
        dialog.exec()
    
    def _on_help_clicked(self) -> None:
        """Show the help dialog."""
        dialog = HelpDialog(self)
        dialog.exec()
    
    def _on_save_clicked(self) -> None:
        """Save current analysis session to file."""
        if not self.snp_records:
            QMessageBox.warning(
                self, "No Data", 
                "No analysis data to save. Please upload a genetic data file first."
            )
            return
        
        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Save Analysis Session",
            "",
            "GenExplore Session (*.gxs);;All Files (*)"
        )
        
        if not filepath:
            return
        
        # Ensure .gxs extension
        if not filepath.endswith('.gxs'):
            filepath += '.gxs'
        metadata = self._build_session_metadata()
        self._start_session_operation(
            "save",
            filepath,
            {
                "snp_records": self.snp_records,
                "gwas_matches": self.all_matches,
                "polygenic_results": self.polygenic_results,
                "metadata": metadata,
            },
        )
    
    def _on_load_clicked(self) -> None:
        """Load a previously saved analysis session."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Load Analysis Session",
            "",
            "GenExplore Session (*.gxs);;All Files (*)"
        )
        
        if not filepath:
            return
        self._start_session_operation("load", filepath)

    def _build_session_metadata(self) -> dict:
        """Assemble metadata for session persistence."""
        metadata = dict(self.current_session_metadata)
        metadata.update({
            "app_version": APP_VERSION,
            "snp_file": self._last_loaded_file,
            "source_filename": Path(self._last_loaded_file).name if self._last_loaded_file else metadata.get('source_filename'),
            "source_format": getattr(self.current_inspection, 'detected_format', metadata.get('source_format')),
            "source_file_size_bytes": getattr(self.current_inspection, 'file_size_bytes', metadata.get('source_file_size_bytes')),
            "analysis_started_at": self.analysis_started_at.isoformat() if self.analysis_started_at else metadata.get('analysis_started_at'),
            "analysis_completed_at": self.analysis_completed_at.isoformat() if self.analysis_completed_at else metadata.get('analysis_completed_at'),
        })
        return metadata

    def _start_session_operation(self, operation: str, filepath: str, payload: Optional[dict] = None) -> None:
        """Run a session save/load operation in the background."""
        if self.session_worker and self.session_worker.isRunning():
            QMessageBox.information(
                self,
                "Operation in progress",
                "Another session operation is already running. Please wait for it to finish.",
            )
            return

        self._current_session_operation = operation
        self._current_session_filepath = filepath
        self._set_session_controls_busy(True)

        title = "Saving session" if operation == "save" else "Loading session"
        self.session_progress_dialog = QProgressDialog(f"{title}...", None, 0, 100, self)
        self.session_progress_dialog.setWindowTitle(title)
        self.session_progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
        self.session_progress_dialog.setMinimumDuration(0)
        self.session_progress_dialog.setAutoClose(False)
        self.session_progress_dialog.setAutoReset(False)
        self.session_progress_dialog.setCancelButton(None)
        self.session_progress_dialog.setValue(0)
        self.session_progress_dialog.show()

        self.status_bar.showMessage(f"{title}...")
        QApplication.processEvents()

        self.session_worker = SessionOperationWorker(operation, filepath, payload)
        self.session_worker.progress.connect(self._on_session_progress)
        self.session_worker.finished.connect(self._on_session_operation_finished)
        self.session_worker.error.connect(self._on_session_operation_error)
        self.session_worker.start()

    def _set_session_controls_busy(self, busy: bool) -> None:
        """Enable or disable session-related actions while a worker runs."""
        self.save_btn.setEnabled(not busy and bool(self.snp_records))
        self.load_btn.setEnabled(not busy)
        self.upload_btn.setEnabled(not busy)
        try:
            self.export_btn.setEnabled(not busy and (bool(self.polygenic_results) or bool(self.all_matches) or bool(self.snp_records)))
        except Exception:
            pass
        self.filters_group.setEnabled(not busy)

    def _close_session_progress_dialog(self) -> None:
        """Close any session progress dialog currently shown."""
        if self.session_progress_dialog is not None:
            self.session_progress_dialog.setValue(100)
            self.session_progress_dialog.close()
            self.session_progress_dialog.deleteLater()
            self.session_progress_dialog = None

    def _on_session_progress(self, value: int, message: str) -> None:
        """Update the session progress dialog and status bar."""
        if self.session_progress_dialog is not None:
            self.session_progress_dialog.setLabelText(message)
            self.session_progress_dialog.setValue(value)
        self.status_bar.showMessage(message)

    def _apply_loaded_session(self, filepath: str, snp_records: List[SNPRecord], gwas_matches: List[GWASMatch], polygenic_results: List[PolygenicResult], metadata: dict) -> None:
        """Apply a loaded session to the current UI state."""
        self.current_session_metadata = dict(metadata or {})
        self._last_loaded_file = metadata.get('snp_file') or metadata.get('source_file') or metadata.get('source_filename') or filepath
        self.analysis_started_at = metadata.get('analysis_started_at')
        self.analysis_completed_at = metadata.get('analysis_completed_at')
        source_format = metadata.get('source_format') or metadata.get('detected_format') or 'Session'
        source_filename = metadata.get('source_filename') or (Path(self._last_loaded_file).name if self._last_loaded_file else Path(filepath).name)
        self.current_inspection = SimpleNamespace(
            filepath=self._last_loaded_file or filepath,
            filename=source_filename,
            detected_format=source_format,
            file_size_bytes=metadata.get('source_file_size_bytes'),
        )

        self.snp_records = snp_records
        self.all_matches = gwas_matches
        self.current_matches = gwas_matches
        self.polygenic_results = polygenic_results
        self.current_page = 0

        self.genotype_status.setText(f"✓ {len(snp_records):,} SNPs loaded")
        self.genotype_status.setStyleSheet("color: #006000; font-weight: bold;")

        self.filters_group.setEnabled(True)
        self._set_session_controls_busy(False)

        self._reset_filters()
        self.polygenic_widget.set_genotype_data(snp_records)
        self.polygenic_widget.display_loaded_results(polygenic_results)

    def _on_session_operation_finished(self, result: object) -> None:
        """Handle completion of a save or load session worker."""
        operation = self._current_session_operation
        filepath = self._current_session_filepath or ""
        self._close_session_progress_dialog()
        self.session_worker = None
        self._set_session_controls_busy(False)

        if operation == "save":
            self.status_bar.showMessage(f"Session saved to {filepath}")
            QMessageBox.information(
                self,
                "Session Saved",
                f"Analysis session saved successfully.\n\n"
                f"File: {filepath}\n"
                f"SNPs: {len(self.snp_records):,}\n"
                f"GWAS matches: {len(self.all_matches):,}\n"
                f"Polygenic scores: {len(self.polygenic_results)}",
            )
        elif operation == "load":
            try:
                snp_records, gwas_matches, polygenic_results, metadata = result  # type: ignore[misc]
                self._apply_loaded_session(filepath, snp_records, gwas_matches, polygenic_results, metadata)

                created_at = metadata.get('created_at', 'Unknown')
                self.status_bar.showMessage(
                    f"✓ Session loaded | {len(snp_records):,} SNPs | "
                    f"{len(gwas_matches):,} matches | {len(polygenic_results)} scores"
                )
                QMessageBox.information(
                    self,
                    "Session Loaded",
                    f"Analysis session loaded successfully.\n\n"
                    f"Created: {created_at}\n"
                    f"SNPs: {len(snp_records):,}\n"
                    f"GWAS matches: {len(gwas_matches):,}\n"
                    f"Polygenic scores: {len(polygenic_results)}",
                )
            except Exception as exc:
                logger.exception("Loaded session result could not be applied: %s", exc)
                QMessageBox.critical(
                    self,
                    "Load Error",
                    f"Failed to apply loaded session:\n{exc}",
                )
                self.status_bar.showMessage("Load failed")
        else:
            self.status_bar.showMessage("Session operation finished")

    def _on_session_operation_error(self, error: object) -> None:
        """Handle save/load session errors."""
        operation = self._current_session_operation
        self._close_session_progress_dialog()
        self.session_worker = None
        self._set_session_controls_busy(False)

        if operation == "load" and isinstance(error, ValueError):
            QMessageBox.critical(
                self,
                "Invalid File",
                f"The file is not a valid GenExplore session:\n{error}",
            )
            self.status_bar.showMessage("Load failed - invalid file")
            return

        logger.error("Failed to %s session: %s", operation or "process", error)
        title = "Save Error" if operation == "save" else "Load Error"
        action = "save" if operation == "save" else "load"
        QMessageBox.critical(
            self,
            title,
            f"Failed to {action} session:\n{error}",
        )
        self.status_bar.showMessage(f"{title.replace(' Error', '')} failed")

    def export_report(self) -> None:
        """Export current analysis to a self-contained HTML report with filtering, sorting and detailed anchors."""
        try:
            filepath, _ = QFileDialog.getSaveFileName(
                self,
                "Export Report",
                "",
                "HTML Files (*.html);;All Files (*)"
            )

            if not filepath:
                return

            if not filepath.lower().endswith('.html'):
                filepath += '.html'

            import html as html_lib

            from database.polygenic_database import DatabaseVersionManager, get_gwas_database_stats

            def _fmt_text(value, default='N/A'):
                if value is None or value == '':
                    return default
                if isinstance(value, Enum):
                    return value.value
                return str(value)

            def _fmt_datetime(value):
                if value is None or value == '':
                    return 'N/A'
                if isinstance(value, datetime):
                    return value.strftime('%Y-%m-%d %H:%M:%S')
                if isinstance(value, str):
                    try:
                        return datetime.fromisoformat(value).strftime('%Y-%m-%d %H:%M:%S')
                    except Exception:
                        return value
                return str(value)

            def _slugify(value):
                slug = re.sub(r'[^a-zA-Z0-9_-]+', '-', str(value)).strip('-').lower()
                return slug or 'item'

            def _safe_float_text(value):
                if value is None or value == '':
                    return 'N/A'
                try:
                    return f'{float(value):.1f}%'
                except Exception:
                    return html_lib.escape(str(value))

            def _human_file_size(size_bytes):
                if size_bytes in (None, '', 0):
                    return 'N/A' if size_bytes in (None, '') else '0 B'
                try:
                    size = float(size_bytes)
                except Exception:
                    return 'N/A'
                units = ['B', 'KB', 'MB', 'GB', 'TB']
                unit_index = 0
                while size >= 1024 and unit_index < len(units) - 1:
                    size /= 1024.0
                    unit_index += 1
                return f'{int(size)} {units[unit_index]}' if unit_index == 0 else f'{size:.2f} {units[unit_index]}'

            def _resolve_current_context():
                metadata = dict(getattr(self, 'current_session_metadata', {}) or {})
                inspection = getattr(self, 'current_inspection', None)
                source_path = (
                    getattr(inspection, 'filepath', None)
                    or metadata.get('snp_file')
                    or self._last_loaded_file
                )
                source_filename = (
                    getattr(inspection, 'filename', None)
                    or metadata.get('source_filename')
                    or (Path(source_path).name if source_path else 'Unknown')
                )
                detected_format = (
                    getattr(inspection, 'detected_format', None)
                    or metadata.get('source_format')
                    or metadata.get('detected_format')
                    or 'Unknown'
                )
                file_size_bytes = (
                    getattr(inspection, 'file_size_bytes', None)
                    or metadata.get('source_file_size_bytes')
                )
                return metadata, source_path, source_filename, detected_format, file_size_bytes

            vm = DatabaseVersionManager()
            gwas_ver = vm.get_gwas_version()
            pgs_ver = vm.get_pgs_version()
            gwas_stats = get_gwas_database_stats() or {}

            metadata, source_path, source_filename, detected_format, file_size_bytes = _resolve_current_context()
            matches = list(getattr(self, 'current_matches', None) or self.all_matches or [])
            poly_results = list(self.polygenic_results or [])
            poly_score_lookup = {
                getattr(score, 'pgs_id', None): score
                for score in getattr(self.polygenic_widget, 'scores', []) or []
                if getattr(score, 'pgs_id', None)
            }

            analysis_started = getattr(self, 'analysis_started_at', None) or metadata.get('analysis_started_at')
            analysis_completed = getattr(self, 'analysis_completed_at', None) or metadata.get('analysis_completed_at')

            css = """
            :root { color-scheme: light; }
            html { scroll-behavior: smooth; }
            body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,'Helvetica Neue',Arial,sans-serif;margin:24px;color:#1f2328;background:#fbfcfe}
            .shell{max-width:1280px;margin:0 auto}
            header{display:flex;align-items:flex-start;justify-content:space-between;gap:20px;padding:18px 20px;background:linear-gradient(135deg,#ffffff 0%,#f2f7ff 100%);border:1px solid #dde6f3;border-radius:16px;box-shadow:0 10px 30px rgba(20,40,80,.05)}
            h1{margin:0;font-size:2rem;letter-spacing:-0.02em}
            h2{margin:0 0 10px 0;font-size:1.35rem}
            h3{margin:18px 0 8px 0}
            .muted{color:#5b6573}
            .badge{display:inline-block;padding:5px 10px;border-radius:999px;background:#e9f2ff;color:#124ea0;font-weight:700}
            .index-links{display:flex;flex-wrap:wrap;gap:8px}
            .index-links a{display:inline-block;padding:8px 12px;border-radius:999px;background:#eef5ff;color:#114a9c;text-decoration:none;font-weight:700}
            .index-links a:hover{background:#dbeafe}
            .grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px;margin-top:16px}
            .card{background:#fff;border:1px solid #e3e8ef;border-radius:14px;padding:16px 18px;box-shadow:0 8px 24px rgba(20,40,80,.04)}
            .card.wide{grid-column:1/-1}
            .section{margin-top:18px}
            .toolbar{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin:18px 0;padding:14px 16px;background:#fff;border:1px solid #e3e8ef;border-radius:14px}
            .toolbar label{display:flex;align-items:center;gap:8px;flex-wrap:wrap}
            input[type=search], input[type=number]{padding:8px 10px;border:1px solid #cbd5e1;border-radius:8px;min-width:220px}
            select{padding:8px 10px;border:1px solid #cbd5e1;border-radius:8px;min-width:220px;background:#fff}
            table{width:100%;border-collapse:separate;border-spacing:0;background:#fff;border:1px solid #e3e8ef;border-radius:14px;overflow:hidden}
            thead th{position:sticky;top:0;background:#f6f8fb;padding:11px 10px;border-bottom:1px solid #e3e8ef;text-align:left;cursor:pointer;white-space:nowrap}
            tbody td{padding:10px;border-bottom:1px solid #eef2f7;vertical-align:top}
            tbody tr:hover{background:#fbfdff}
            .sortable::after{content:' ↕';color:#8a94a3;font-weight:400}
            .details-link,.action-link{display:inline-block;padding:6px 10px;border-radius:999px;background:#0b5fff;color:#fff;text-decoration:none;font-weight:700;line-height:1}
            .details-link:hover,.action-link:hover{background:#084bc0}
            .details-link.secondary{background:#0f9d58}
            .details-link.secondary:hover{background:#0b8043}
            .tag{display:inline-block;padding:3px 8px;border-radius:999px;background:#f1f5f9;color:#334155;font-size:.85rem}
            .small{font-size:.92rem;color:#5b6573}
            .detail-anchor{scroll-margin-top:18px}
            .detail-block{margin-top:14px;padding:16px 18px;border:1px solid #e3e8ef;border-radius:14px;background:#fff}
            .detail-meta{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin-top:10px}
            .detail-meta div{padding:10px 12px;background:#f8fafc;border:1px solid #edf2f7;border-radius:10px}
            .detail-meta strong{display:block;margin-bottom:3px}
            .detail-links{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0 4px 0}
            .detail-note{margin-top:10px;padding:10px 12px;border-left:4px solid #0b5fff;background:#f4f8ff;color:#1f3a5f;border-radius:8px}
            .contributors{margin-top:10px;width:100%}
            .contributors th{cursor:default}
            footer{margin:26px 0 10px;color:#667085;font-size:.9rem}
            @media print{.toolbar,.details-link,.action-link,header .actions{display:none}}
            """

            js = """
            function sortTable(table){
                const headers = Array.from(table.querySelectorAll('thead th.sortable'));
                headers.forEach((header, index) => {
                    header.addEventListener('click', () => {
                        const tbody = table.querySelector('tbody');
                        const rows = Array.from(tbody.querySelectorAll('tr'));
                        const numeric = header.dataset.type === 'numeric';
                        const nextOrder = header.dataset.order === 'asc' ? 'desc' : 'asc';
                        headers.forEach(h => { delete h.dataset.order; });
                        header.dataset.order = nextOrder;
                        rows.sort((left, right) => {
                            const leftCell = left.children[index];
                            const rightCell = right.children[index];
                            const leftValue = leftCell.dataset.value ?? leftCell.innerText;
                            const rightValue = rightCell.dataset.value ?? rightCell.innerText;
                            if (numeric) {
                                const leftNumber = parseFloat(leftValue) || 0;
                                const rightNumber = parseFloat(rightValue) || 0;
                                return nextOrder === 'asc' ? leftNumber - rightNumber : rightNumber - leftNumber;
                            }
                            return nextOrder === 'asc'
                                ? String(leftValue).localeCompare(String(rightValue))
                                : String(rightValue).localeCompare(String(leftValue));
                        });
                        rows.forEach(row => tbody.appendChild(row));
                    });
                });
            }

            function applyFilters(){
                const gwasQuery = (document.getElementById('gwasSearch')?.value || '').toLowerCase();
                const gwasMinImpact = parseFloat(document.getElementById('gwasMinImpact')?.value || '0') || 0;
                const gwasCategory = (document.getElementById('gwasCategory')?.value || 'ALL').toLowerCase();

                document.querySelectorAll('#gwasTable tbody tr').forEach(row => {
                    const text = row.innerText.toLowerCase();
                    const impact = parseFloat(row.dataset.impact || '0') || 0;
                    const category = (row.dataset.category || '').toLowerCase();
                    const matchesQuery = !gwasQuery || text.includes(gwasQuery);
                    const matchesImpact = impact >= gwasMinImpact;
                    const matchesCategory = gwasCategory === 'all' || category === gwasCategory;
                    row.style.display = matchesQuery && matchesImpact && matchesCategory ? '' : 'none';
                });

                const polyQuery = (document.getElementById('polySearch')?.value || '').toLowerCase();
                const polyMinCoverage = parseFloat(document.getElementById('polyMinCoverage')?.value || '0') || 0;
                const polyRisk = (document.getElementById('polyRiskCategory')?.value || 'ALL').toLowerCase();
                const polyTraitCategory = (document.getElementById('polyTraitCategory')?.value || 'ALL').toLowerCase();

                document.querySelectorAll('#polygenicTable tbody tr').forEach(row => {
                    const text = row.innerText.toLowerCase();
                    const coverage = parseFloat(row.dataset.coverage || '0') || 0;
                    const risk = (row.dataset.riskCategory || '').toLowerCase();
                    const traitCategory = (row.dataset.traitCategory || '').toLowerCase();
                    const matchesQuery = !polyQuery || text.includes(polyQuery);
                    const matchesCoverage = coverage >= polyMinCoverage;
                    const matchesRisk = polyRisk === 'all' || risk === polyRisk;
                    const matchesTraitCategory = polyTraitCategory === 'all' || traitCategory === polyTraitCategory;
                    row.style.display = matchesQuery && matchesCoverage && matchesRisk && matchesTraitCategory ? '' : 'none';
                });
            }

            document.addEventListener('DOMContentLoaded', () => {
                const gwasSearch = document.getElementById('gwasSearch');
                const gwasMinImpact = document.getElementById('gwasMinImpact');
                const gwasCategory = document.getElementById('gwasCategory');
                const polySearch = document.getElementById('polySearch');
                const polyMinCoverage = document.getElementById('polyMinCoverage');
                const polyRiskCategory = document.getElementById('polyRiskCategory');
                const polyTraitCategory = document.getElementById('polyTraitCategory');

                [gwasSearch, gwasMinImpact, gwasCategory, polySearch, polyMinCoverage, polyRiskCategory, polyTraitCategory].forEach(control => {
                    if (control) control.addEventListener('input', applyFilters);
                    if (control) control.addEventListener('change', applyFilters);
                });

                const gwasTable = document.getElementById('gwasTable');
                const polygenicTable = document.getElementById('polygenicTable');
                if (gwasTable) sortTable(gwasTable);
                if (polygenicTable) sortTable(polygenicTable);
                applyFilters();
            });
            """

            def _gwas_impact(match):
                impact = getattr(match, 'impact_score', None)
                return impact if impact is not None else getattr(match, 'impact', None)

            def _fmt_risk_label(value):
                if isinstance(value, Enum):
                    return value.value
                if value is None or value == '':
                    return 'N/A'
                return str(value)

            parts = []
            parts.append('<!doctype html>')
            parts.append('<html lang="en">')
            parts.append('<head>')
            parts.append('<meta charset="utf-8">')
            parts.append('<meta name="viewport" content="width=device-width,initial-scale=1">')
            parts.append(f'<title>{html_lib.escape(APP_NAME)} Report</title>')
            parts.append(f'<style>{css}</style>')
            parts.append('</head>')
            parts.append('<body>')
            parts.append('<div class="shell">')
            parts.append('<header>')
            parts.append('<div>')
            parts.append(f'<h1>{html_lib.escape(APP_NAME)} — Genetic Exploration Report</h1>')
            parts.append(f'<div class="muted">Generated on {html_lib.escape(_fmt_datetime(datetime.now()))}</div>')
            parts.append('</div>')
            parts.append('<div class="actions">')
            parts.append(f'<span class="badge">{len(matches):,} GWAS matches</span>')
            parts.append(f' <span class="badge">{len(poly_results):,} polygenic scores</span>')
            parts.append('</div>')
            parts.append('</header>')

            parts.append('<section class="card wide" id="index">')
            parts.append('<h2>Index</h2>')
            parts.append('<div class="index-links">')
            for anchor, label in [
                ('#sample-info', 'Sample Information'),
                ('#summary', 'Summary'),
                ('#gwas-section', 'Monogenic Analysis'),
                ('#polygenic-section', 'Polygenic Analysis'),
                ('#methodology', 'Methodology & Limitations'),
                ('#details', 'Details'),
            ]:
                parts.append(f'<a href="{anchor}">{html_lib.escape(label)}</a>')
            parts.append('</div>')
            parts.append('</section>')

            parts.append('<div class="grid">')
            parts.append('<section class="card">')
            parts.append('<h2>Sample Information</h2>')
            parts.append('<div class="detail-meta">')
            parts.append(f'<div><strong>File</strong>{html_lib.escape(_fmt_text(source_filename, "Unknown"))}</div>')
            parts.append(f'<div><strong>Path</strong>{html_lib.escape(_fmt_text(source_path, "Unknown"))}</div>')
            parts.append(f'<div><strong>Detected format</strong>{html_lib.escape(_fmt_text(detected_format, "Unknown"))}</div>')
            parts.append(f'<div><strong>File size</strong>{html_lib.escape(_human_file_size(file_size_bytes))}</div>')
            parts.append(f'<div><strong>Analysis started</strong>{html_lib.escape(_fmt_datetime(analysis_started))}</div>')
            parts.append(f'<div><strong>Analysis completed</strong>{html_lib.escape(_fmt_datetime(analysis_completed))}</div>')
            parts.append('</div>')
            parts.append('</section>')

            parts.append('<section class="card">')
            parts.append('<h2>Summary</h2>')
            parts.append('<div class="detail-meta">')
            if gwas_ver:
                parts.append(f'<div><strong>GWAS DB version</strong>{html_lib.escape(_fmt_text(gwas_ver.version))}</div>')
                parts.append(f'<div><strong>GWAS records</strong>{html_lib.escape(f"{gwas_ver.record_count:,}")}</div>')
            if pgs_ver:
                parts.append(f'<div><strong>PGS DB version</strong>{html_lib.escape(_fmt_text(pgs_ver.version))}</div>')
                parts.append(f'<div><strong>PGS records</strong>{html_lib.escape(f"{pgs_ver.record_count:,}")}</div>')
            if gwas_stats.get('variants') is not None:
                parts.append(f'<div><strong>GWAS variants</strong>{html_lib.escape(f"{gwas_stats.get('variants'):,}")}</div>')
            parts.append(f'<div><strong>Monogenic matches</strong>{html_lib.escape(f"{len(matches):,}")}</div>')
            parts.append(f'<div><strong>Polygenic scores</strong>{html_lib.escape(f"{len(poly_results):,}")}</div>')
            parts.append('</div>')
            parts.append('</section>')
            parts.append('</div>')

            parts.append('<section class="section" id="gwas-section">')
            parts.append('<h2>GWAS Matches</h2>')
            parts.append(f'<p class="small">Found {len(matches):,} monogenic matches.</p>')
            parts.append('<section class="toolbar">')
            parts.append('<label>Search <input id="gwasSearch" type="search" placeholder="Search traits, genes, SNP IDs"></label>')
            parts.append('<label>Min. Impact Score <input id="gwasMinImpact" type="number" min="0" max="10" step="0.1" value="0"></label>')
            parts.append('<label>Category <select id="gwasCategory"><option value="ALL">All</option>')
            for category_label in TRAIT_CATEGORIES:
                parts.append(f'<option value="{html_lib.escape(category_label)}">{html_lib.escape(category_label)}</option>')
            parts.append('</select></label>')
            parts.append('</section>')
            parts.append('<table id="gwasTable">')
            parts.append('<thead><tr>')
            for label, type_name in [
                ('SNP ID', 'text'), ('Gene', 'text'), ('Trait', 'text'), ('User Genotype', 'text'),
                ('Risk Allele', 'text'), ('P-value', 'numeric'), ('Category', 'text'), ('Impact Score', 'numeric'),
                ('Interpretation', 'text'), ('Details', 'text')
            ]:
                parts.append(f'<th class="sortable" data-type="{type_name}">{html_lib.escape(label)}</th>')
            parts.append('</tr></thead><tbody>')

            for match in matches:
                rsid = getattr(match, 'rsid', '')
                detail_id = f'gwas-detail-{_slugify(rsid or match.trait)}'
                impact_score = _gwas_impact(match)
                interpretation = get_score_interpretation(impact_score) if impact_score is not None else 'N/A'
                p_value = getattr(match, 'p_value', None)
                p_value_text = f'{p_value:.2e}' if p_value is not None else 'N/A'
                category_text = _fmt_text(getattr(match, 'category', None))
                parts.append(f'<tr data-impact="{html_lib.escape(str(impact_score if impact_score is not None else 0))}" data-category="{html_lib.escape(category_text.lower())}">')
                parts.append(f'<td>{html_lib.escape(_fmt_text(rsid))}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(getattr(match, "gene", None), "-"))}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(getattr(match, "trait", None)))}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(getattr(match, "user_genotype", None)))}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(getattr(match, "risk_allele", None)))}</td>')
                parts.append(f'<td data-value="{html_lib.escape(str(p_value if p_value is not None else 0))}">{html_lib.escape(p_value_text)}</td>')
                parts.append(f'<td>{html_lib.escape(category_text)}</td>')
                parts.append(f'<td data-value="{html_lib.escape(str(impact_score if impact_score is not None else 0))}">{html_lib.escape(f"{impact_score:.2f}" if impact_score is not None else "N/A")}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(interpretation))}</td>')
                parts.append(f'<td><a class="action-link" href="#{detail_id}">Explain</a></td>')
                parts.append('</tr>')

            parts.append('</tbody></table>')
            parts.append('</section>')

            parts.append('<section class="section" id="polygenic-section">')
            parts.append('<h2>Polygenic Risk Scores</h2>')
            parts.append(f'<p class="small">Found {len(poly_results):,} computed scores.</p>')
            parts.append('<section class="toolbar">')
            parts.append('<label>Search <input id="polySearch" type="search" placeholder="Search traits, PGS IDs"></label>')
            parts.append('<label>Min. Coverage <input id="polyMinCoverage" type="number" min="0" max="100" step="0.1" value="0"></label>')
            parts.append('<label>Risk Category <select id="polyRiskCategory"><option value="ALL">All</option><option value="LOW">Low</option><option value="INTERMEDIATE">Intermediate</option><option value="HIGH">High</option></select></label>')
            parts.append('<label>Trait Category <select id="polyTraitCategory"><option value="ALL">All</option>')
            for category_label in TRAIT_CATEGORIES:
                parts.append(f'<option value="{html_lib.escape(category_label)}">{html_lib.escape(category_label)}</option>')
            parts.append('</select></label>')
            parts.append('</section>')
            parts.append('<table id="polygenicTable">')
            parts.append('<thead><tr>')
            for label, type_name in [
                ('PGS ID', 'text'), ('Trait', 'text'), ('Category', 'text'), ('Variants', 'numeric'),
                ('Coverage', 'numeric'), ('Percentile', 'numeric'), ('Risk', 'text'), ('Details', 'text')
            ]:
                parts.append(f'<th class="sortable" data-type="{type_name}">{html_lib.escape(label)}</th>')
            parts.append('</tr></thead><tbody>')

            for result in poly_results:
                score = getattr(result, 'score', None) or poly_score_lookup.get(getattr(result, 'pgs_id', None))
                pgs_id = getattr(result, 'pgs_id', None) or getattr(score, 'pgs_id', None)
                trait_name = getattr(result, 'trait_name', None) or getattr(score, 'trait_name', None) or getattr(score, 'trait_reported', None)
                trait_category = getattr(result, 'trait_category', None) or getattr(score, 'trait_category', None)
                variants_found = getattr(result, 'variants_found', None)
                variants_total = getattr(result, 'variants_total', None)
                coverage_percent = getattr(result, 'coverage_percent', None)
                percentile = getattr(result, 'percentile', None)
                risk_category = getattr(result, 'risk_category', None)
                detail_id = f'poly-detail-{_slugify(pgs_id or trait_name)}'

                category_text = _fmt_text(trait_category)
                risk_text = _fmt_risk_label(risk_category)
                variants_text = f'{variants_found:,} / {variants_total:,}' if variants_found is not None and variants_total is not None else _fmt_text(variants_found)
                coverage_value = coverage_percent if coverage_percent is not None else 0
                percentile_value = percentile if percentile is not None else 0

                parts.append(f'<tr data-coverage="{html_lib.escape(str(coverage_value))}" data-risk-category="{html_lib.escape(risk_text.lower())}" data-trait-category="{html_lib.escape(category_text.lower())}">')
                parts.append(f'<td>{html_lib.escape(_fmt_text(pgs_id))}</td>')
                parts.append(f'<td>{html_lib.escape(_fmt_text(trait_name))}</td>')
                parts.append(f'<td>{html_lib.escape(category_text)}</td>')
                parts.append(f'<td data-value="{html_lib.escape(str(variants_found if variants_found is not None else 0))}">{html_lib.escape(variants_text)}</td>')
                parts.append(f'<td data-value="{html_lib.escape(str(coverage_value))}">{html_lib.escape(f"{coverage_value:.0f}%" if coverage_percent is not None else "N/A")}</td>')
                parts.append(f'<td data-percentile="{html_lib.escape(str(percentile_value))}" data-value="{html_lib.escape(str(percentile_value))}">{html_lib.escape(_safe_float_text(percentile))}</td>')
                parts.append(f'<td>{html_lib.escape(risk_text)}</td>')
                parts.append(f'<td><a class="details-link secondary" href="#{detail_id}">View</a></td>')
                parts.append('</tr>')

            parts.append('</tbody></table>')
            parts.append('</section>')

            parts.append('<section class="section" id="methodology">')
            parts.append('<h2>Methodology & Limitations</h2>')
            parts.append('<div class="card">')
            parts.append('<p class="small">')
            parts.append('This report combines your genotype file with local GWAS and PGS catalog data. ')
            parts.append('Monogenic matches summarize individual variant associations, while polygenic scores are population-based estimates.')
            parts.append('</p>')
            parts.append('<p class="small">')
            parts.append('Low coverage reduces confidence in a score. Risk categories are derived from percentile thresholds and may not transfer across ancestries.')
            parts.append('</p>')
            parts.append('</div>')
            parts.append('</section>')

            total_pgs_variants = sum(
                getattr(score, 'num_variants', 0) or 0
                for score in getattr(self.polygenic_widget, 'scores', []) or []
            )
            pgs_scores_loaded = len(getattr(self.polygenic_widget, 'scores', []) or [])

            parts.append('<section class="section" id="database-stats">')
            parts.append('<h2>Database Statistics</h2>')
            parts.append('<div class="grid">')
            parts.append('<div class="card">')
            parts.append('<h3>GWAS Catalog</h3>')
            parts.append('<div class="detail-meta">')
            parts.append(f'<div><strong>Version</strong>{html_lib.escape(_fmt_text(getattr(gwas_ver, "version", None), "N/A"))}</div>')
            parts.append(f'<div><strong>Records</strong>{html_lib.escape(f"{getattr(gwas_ver, 'record_count', 0):,}" if gwas_ver else "N/A")}</div>')
            parts.append(f'<div><strong>Variants</strong>{html_lib.escape(f"{gwas_stats.get('variants', 0):,}" if gwas_stats.get('variants') is not None else "N/A")}</div>')
            parts.append(f'<div><strong>Traits</strong>{html_lib.escape(f"{gwas_stats.get('traits', 0):,}" if gwas_stats.get('traits') is not None else "N/A")}</div>')
            parts.append(f'<div><strong>Genes</strong>{html_lib.escape(f"{gwas_stats.get('genes', 0):,}" if gwas_stats.get('genes') is not None else "N/A")}</div>')
            parts.append(f'<div><strong>Associations</strong>{html_lib.escape(f"{gwas_stats.get('associations', 0):,}" if gwas_stats.get('associations') is not None else "N/A")}</div>')
            parts.append('</div>')
            parts.append('</div>')
            parts.append('<div class="card">')
            parts.append('<h3>PGS Catalog</h3>')
            parts.append('<div class="detail-meta">')
            parts.append(f'<div><strong>Version</strong>{html_lib.escape(_fmt_text(getattr(pgs_ver, "version", None), "N/A"))}</div>')
            parts.append(f'<div><strong>Scores loaded</strong>{html_lib.escape(f"{pgs_scores_loaded:,}" if pgs_scores_loaded else "0")}</div>')
            parts.append(f'<div><strong>Variant definitions</strong>{html_lib.escape(f"{total_pgs_variants:,}" if total_pgs_variants else "0")}</div>')
            parts.append(f'<div><strong>Sample source</strong>{html_lib.escape(_fmt_text(getattr(pgs_ver, "source_type", None), "catalog"))}</div>')
            parts.append('</div>')
            parts.append('</div>')
            parts.append('</div>')
            parts.append('</section>')

            parts.append('<section class="section" id="details">')
            parts.append('<h2>Details</h2>')
            parts.append('<p class="small">Use the links above to jump to each card. These anchors restore the original Explain/View workflow.</p>')

            for match in matches:
                rsid = getattr(match, 'rsid', '')
                detail_id = f'gwas-detail-{_slugify(rsid or match.trait)}'
                impact_score = _gwas_impact(match)
                parts.append(f'<div class="detail-block detail-anchor" id="{detail_id}">')
                parts.append(f'<h3>Explain: {html_lib.escape(_fmt_text(rsid))}</h3>')
                detail_links = []
                if rsid:
                    detail_links.append(f'<a class="details-link" href="https://www.ncbi.nlm.nih.gov/snp/{html_lib.escape(str(rsid))}" target="_blank" rel="noopener noreferrer">dbSNP</a>')
                    detail_links.append(f'<a class="details-link" href="https://www.ebi.ac.uk/gwas/variants/{html_lib.escape(str(rsid))}" target="_blank" rel="noopener noreferrer">GWAS Catalog</a>')
                gene_name = getattr(match, 'gene', None)
                if gene_name:
                    detail_links.append(f'<a class="details-link secondary" href="https://www.genecards.org/Search/Keyword?queryString={html_lib.escape(str(gene_name))}" target="_blank" rel="noopener noreferrer">GeneCards</a>')
                if detail_links:
                    parts.append('<div class="detail-links">' + ' '.join(detail_links) + '</div>')
                parts.append('<div class="detail-meta">')
                for label, value in [
                    ('Gene', getattr(match, 'gene', None)),
                    ('Trait', getattr(match, 'trait', None)),
                    ('Chromosome', getattr(match, 'chromosome', None)),
                    ('Position', f'{getattr(match, "position", None):,}' if getattr(match, 'position', None) is not None else None),
                    ('User genotype', getattr(match, 'user_genotype', None)),
                    ('Risk allele', getattr(match, 'risk_allele', None)),
                    ('P-value', f'{getattr(match, "p_value", 0):.2e}' if getattr(match, 'p_value', None) is not None else None),
                    ('Category', getattr(match, 'category', None)),
                    ('Impact score', f'{impact_score:.2f}' if impact_score is not None else None),
                    ('Interpretation', get_score_interpretation(impact_score) if impact_score is not None else None),
                    ('Odds ratio', f'{getattr(match, "odds_ratio", None):.2f}' if getattr(match, 'odds_ratio', None) is not None else None),
                    ('Sample size', f'{getattr(match, "sample_size", None):,}' if getattr(match, 'sample_size', None) is not None else None),
                    ('Allele frequency', f'{getattr(match, "allele_frequency", None):.1%}' if getattr(match, 'allele_frequency', None) is not None else None),
                ]:
                    parts.append(f'<div><strong>{html_lib.escape(label)}</strong>{html_lib.escape(_fmt_text(value))}</div>')
                parts.append('</div>')
                parts.append('<div class="detail-note">')
                if getattr(match, 'has_risk_allele', lambda: False)():
                    parts.append('You carry the reported risk allele for this association. This is not diagnostic and should be interpreted in context.')
                else:
                    parts.append('You do not carry the reported risk allele for this association.')
                parts.append('</div>')
                parts.append('<p><a class="action-link" href="#gwas-section">Back to GWAS table</a></p>')
                parts.append('</div>')

            for result in poly_results:
                score = getattr(result, 'score', None) or poly_score_lookup.get(getattr(result, 'pgs_id', None))
                pgs_id = getattr(result, 'pgs_id', None) or getattr(score, 'pgs_id', None)
                trait_name = getattr(result, 'trait_name', None) or getattr(score, 'trait_name', None) or getattr(score, 'trait_reported', None)
                detail_id = f'poly-detail-{_slugify(pgs_id or trait_name)}'
                parts.append(f'<div class="detail-block detail-anchor" id="{detail_id}">')
                parts.append(f'<h3>View: {html_lib.escape(_fmt_text(trait_name))}</h3>')
                detail_links = []
                publication_doi = getattr(result, 'publication_doi', None) or (getattr(score, 'publication_doi', None) if score else None)
                score_pgs_id = pgs_id or (getattr(score, 'pgs_id', None) if score else None)
                if publication_doi:
                    detail_links.append(f'<a class="details-link" href="https://doi.org/{html_lib.escape(str(publication_doi))}" target="_blank" rel="noopener noreferrer">DOI</a>')
                if score_pgs_id:
                    detail_links.append(f'<a class="details-link secondary" href="https://www.pgscatalog.org/score/{html_lib.escape(str(score_pgs_id))}/" target="_blank" rel="noopener noreferrer">PGS Catalog</a>')
                if trait_name:
                    trait_search = str(trait_name).replace(' ', '+')
                    detail_links.append(f'<a class="details-link" href="https://pubmed.ncbi.nlm.nih.gov/?term={html_lib.escape(trait_search)}+polygenic" target="_blank" rel="noopener noreferrer">PubMed</a>')
                if detail_links:
                    parts.append('<div class="detail-links">' + ' '.join(detail_links) + '</div>')
                parts.append('<div class="detail-meta">')
                for label, value in [
                    ('PGS ID', pgs_id),
                    ('Category', getattr(result, 'trait_category', None) or getattr(score, 'trait_category', None)),
                    ('Raw score', f'{getattr(result, "raw_score", 0):.4f}' if getattr(result, 'raw_score', None) is not None else None),
                    ('Normalized score', f'{getattr(result, "normalized_score", 0):.2f}' if getattr(result, 'normalized_score', None) is not None else None),
                    ('Percentile', _safe_float_text(getattr(result, 'percentile', None))),
                    ('Risk category', _fmt_risk_label(getattr(result, 'risk_category', None))),
                    ('Variants found', f'{getattr(result, "variants_found", 0):,}' if getattr(result, 'variants_found', None) is not None else None),
                    ('Variants total', f'{getattr(result, "variants_total", 0):,}' if getattr(result, 'variants_total', None) is not None else None),
                    ('Coverage', f'{getattr(result, "coverage_percent", 0):.0f}%' if getattr(result, 'coverage_percent', None) is not None else None),
                    ('Population reference', getattr(result, 'population_reference', None)),
                    ('Computation time', f'{getattr(result, "computation_time_ms", 0):.0f} ms' if getattr(result, 'computation_time_ms', None) is not None else None),
                ]:
                    parts.append(f'<div><strong>{html_lib.escape(label)}</strong>{html_lib.escape(_fmt_text(value))}</div>')
                parts.append('</div>')

                if score:
                    parts.append('<h4>Scientific Context</h4>')
                    parts.append('<div class="detail-meta">')
                    for label, value in [
                        ('Trait name', getattr(score, 'trait_name', None)),
                        ('Publication year', getattr(score, 'publication_year', None)),
                        ('Population', getattr(score, 'study_population', None)),
                        ('Sample size', f'{getattr(score, "sample_size", 0):,}' if getattr(score, 'sample_size', None) is not None else None),
                        ('Variant count', f'{getattr(score, "num_variants", 0):,}' if getattr(score, 'num_variants', None) is not None else None),
                        ('Description', getattr(score, 'description', None)),
                    ]:
                        parts.append(f'<div><strong>{html_lib.escape(label)}</strong>{html_lib.escape(_fmt_text(value))}</div>')
                    parts.append('</div>')

                if getattr(result, 'variant_contributions', None):
                    parts.append('<h4>Top Contributing Variants</h4>')
                    parts.append('<table class="contributors">')
                    parts.append('<thead><tr><th>Variant</th><th>Contribution</th></tr></thead><tbody>')
                    for rsid, contribution in result.get_top_contributors(10):
                        parts.append(f'<tr><td>{html_lib.escape(_fmt_text(rsid))}</td><td>{html_lib.escape(f"{contribution:.4f}")}</td></tr>')
                    parts.append('</tbody></table>')

                parts.append('<p><a class="action-link" href="#polygenic-section">Back to polygenic table</a></p>')
                parts.append('</div>')

            parts.append('</section>')
            parts.append('<footer>')
            parts.append(f'Report generated by {html_lib.escape(APP_NAME)} v{html_lib.escape(str(APP_VERSION))}.')
            parts.append(f' <a href="https://github.com/pmpfe/genexplore" target="_blank" rel="noopener noreferrer">Project on GitHub</a>')
            parts.append('</footer>')
            parts.append('</div>')
            parts.append(f'<script>{js}</script>')
            parts.append('</body></html>')

            html = '\n'.join(parts)
            with open(filepath, 'w', encoding='utf-8') as fh:
                fh.write(html)

            self.status_bar.showMessage(f'Report exported to {filepath}')
            QMessageBox.information(self, 'Export Complete', f'Report exported to:\n{filepath}')

        except Exception as e:
            logger.error(f'Failed to export report: {e}')
            QMessageBox.critical(self, 'Export Failed', f'Failed to export report:\n{str(e)}')
    
    def _on_prev_page(self) -> None:
        """Go to previous page."""
        if self.current_page > 0:
            self.current_page -= 1
            self._update_table()
    
    def _on_next_page(self) -> None:
        """Go to next page."""
        total_pages = (len(self.filtered_matches) + RESULTS_PER_PAGE - 1) // RESULTS_PER_PAGE
        if self.current_page < total_pages - 1:
            self.current_page += 1
            self._update_table()
    
    def closeEvent(self, event) -> None:
        """Handle window close event."""
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.worker.wait()
        
        # Also stop polygenic worker
        if (hasattr(self.polygenic_widget, 'worker') and 
            self.polygenic_widget.worker and 
            self.polygenic_widget.worker.isRunning()):
            self.polygenic_widget.worker.cancel()
            self.polygenic_widget.worker.wait()
        
        event.accept()
