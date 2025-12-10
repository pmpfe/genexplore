"""
Main window UI for the Genetic Analysis Application.

Implements PyQt6 interface with tabbed layout for monogenic and polygenic analysis.
"""

from typing import List, Optional
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QTableWidget, QTableWidgetItem, QFileDialog, QLabel, QLineEdit,
    QComboBox, QSlider, QProgressBar, QStatusBar, QMessageBox,
    QHeaderView, QGroupBox, QSpinBox, QFrame, QSplitter, QApplication,
    QDialog, QTextEdit, QScrollArea, QDialogButtonBox, QTabWidget
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QUrl
from PyQt6.QtGui import QFont, QColor, QDesktopServices

from models.data_models import SNPRecord, GWASMatch, FilterCriteria
from models.polygenic_models import PolygenicResult
from backend.parsers import Parser23andMe, ParseError
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
    Background worker thread for processing 23andMe files.
    
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
            
            parser = Parser23andMe()
            
            # Set up progress callback for file parsing
            def parsing_progress(lines_processed: int, total_lines: int) -> None:
                if total_lines > 0:
                    progress_percent = int((lines_processed / total_lines) * 100)
                    self.file_progress.emit(progress_percent, 
                                           f"{lines_processed:,}/{total_lines:,} lines")
            
            parser.set_progress_callback(parsing_progress)
            snp_records = parser.parse_file(self.filepath)
            stats = parser.get_parse_stats()
            
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
        
        # Scrollable text area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        content = QWidget()
        content_layout = QVBoxLayout(content)
        
        help_text = QTextEdit()
        help_text.setReadOnly(True)
        help_text.setHtml(self._get_help_content())
        help_text.setMinimumHeight(500)
        
        content_layout.addWidget(help_text)
        scroll.setWidget(content)
        layout.addWidget(scroll)
        
        # Close button
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)
    
    def _get_help_content(self) -> str:
        return """
        <h1>Genetic Analysis Application</h1>
        
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
            <li><b>Upload your data:</b> Click "Upload 23andMe File" and select your raw data file</li>
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
        self.filtered_matches: List[GWASMatch] = []
        self.current_page = 0
        self.worker: Optional[ProcessingWorker] = None
        self.snp_records: List[SNPRecord] = []
        self.polygenic_results: List[PolygenicResult] = []
        
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
        
        self.load_btn = QPushButton("📂 Load")
        self.load_btn.setMinimumWidth(80)
        self.load_btn.setMinimumHeight(40)
        self.load_btn.setFont(QFont('Arial', 11))
        self.load_btn.setToolTip("Load previously saved analysis")
        header_layout.addWidget(self.load_btn)
        
        self.upload_btn = QPushButton("📁 Upload 23andMe File")
        self.upload_btn.setMinimumWidth(200)
        self.upload_btn.setMinimumHeight(40)
        self.upload_btn.setFont(QFont('Arial', 11))
        header_layout.addWidget(self.upload_btn)
        
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
        self.status_bar.showMessage("Ready - Upload a 23andMe file to begin")
        
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
        self.pvalue_slider.setValue(100)
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
            "Select 23andMe Raw Data File",
            "",
            "Text Files (*.txt);;All Files (*)"
        )
        
        if not filepath:
            return
        
        self._start_processing(filepath)
    
    def _start_processing(self, filepath: str) -> None:
        """Start processing a 23andMe file in background."""
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
        self.snp_records = snp_records
        self.current_page = 0
        
        logger.info(f"Processing complete: {stats}")
        
        self._reset_filters()
        
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
                "No analysis data to save. Please upload a 23andMe file first."
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
        
        try:
            self.status_bar.showMessage("Saving session...")
            QApplication.processEvents()
            
            SessionManager.save_session(
                filepath=filepath,
                snp_records=self.snp_records,
                gwas_matches=self.all_matches,
                polygenic_results=self.polygenic_results,
                metadata={
                    "app_version": APP_VERSION,
                    "snp_file": getattr(self, '_last_loaded_file', None)
                }
            )
            
            self.status_bar.showMessage(f"Session saved to {filepath}")
            QMessageBox.information(
                self, "Session Saved",
                f"Analysis session saved successfully.\n\n"
                f"File: {filepath}\n"
                f"SNPs: {len(self.snp_records):,}\n"
                f"GWAS matches: {len(self.all_matches):,}\n"
                f"Polygenic scores: {len(self.polygenic_results)}"
            )
            
        except Exception as e:
            logger.error(f"Failed to save session: {e}")
            QMessageBox.critical(
                self, "Save Error",
                f"Failed to save session:\n{str(e)}"
            )
            self.status_bar.showMessage("Save failed")
    
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
        
        try:
            self.status_bar.showMessage("Loading session...")
            QApplication.processEvents()
            
            snp_records, gwas_matches, polygenic_results, metadata = \
                SessionManager.load_session(filepath)
            
            # Update internal state
            self.snp_records = snp_records
            self.all_matches = gwas_matches
            self.polygenic_results = polygenic_results
            self.current_page = 0
            
            # Update UI
            self.genotype_status.setText(f"✓ {len(snp_records):,} SNPs loaded")
            self.genotype_status.setStyleSheet("color: #006000; font-weight: bold;")
            
            # Enable filters and save button
            self.filters_group.setEnabled(True)
            self.save_btn.setEnabled(True)
            
            # Refresh monogenic results table
            self._reset_filters()
            
            # Update polygenic widget with loaded data
            self.polygenic_widget.set_genotype_data(snp_records)
            self.polygenic_widget.display_loaded_results(polygenic_results)
            
            # Show summary
            created_at = metadata.get('created_at', 'Unknown')
            self.status_bar.showMessage(
                f"✓ Session loaded | {len(snp_records):,} SNPs | "
                f"{len(gwas_matches):,} matches | {len(polygenic_results)} scores"
            )
            
            QMessageBox.information(
                self, "Session Loaded",
                f"Analysis session loaded successfully.\n\n"
                f"Created: {created_at}\n"
                f"SNPs: {len(snp_records):,}\n"
                f"GWAS matches: {len(gwas_matches):,}\n"
                f"Polygenic scores: {len(polygenic_results)}"
            )
            
        except ValueError as e:
            QMessageBox.critical(
                self, "Invalid File",
                f"The file is not a valid GenExplore session:\n{str(e)}"
            )
            self.status_bar.showMessage("Load failed - invalid file")
        except Exception as e:
            logger.error(f"Failed to load session: {e}")
            QMessageBox.critical(
                self, "Load Error",
                f"Failed to load session:\n{str(e)}"
            )
            self.status_bar.showMessage("Load failed")
    
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
