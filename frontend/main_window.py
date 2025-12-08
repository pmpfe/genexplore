"""
Main window UI for the Genetic Analysis Application.

Implements PyQt6 interface with file upload, results table, filters, and search.
"""

from typing import List, Optional
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QTableWidget, QTableWidgetItem, QFileDialog, QLabel, QLineEdit,
    QComboBox, QSlider, QProgressBar, QStatusBar, QMessageBox,
    QHeaderView, QGroupBox, QSpinBox, QFrame, QSplitter, QApplication,
    QDialog, QTextEdit, QScrollArea, QDialogButtonBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QUrl
from PyQt6.QtGui import QFont, QColor, QDesktopServices

from models.data_models import SNPRecord, GWASMatch, FilterCriteria
from backend.parsers import Parser23andMe, ParseError
from backend.search_engine import SearchEngine, DatabaseError
from backend.scoring import get_score_interpretation
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
        finished: Emitted when processing is complete with results.
        error: Emitted when an error occurs.
        progress: Emitted with progress updates (0-100).
    """
    finished = pyqtSignal(list, dict)
    error = pyqtSignal(str)
    progress = pyqtSignal(int, str)
    
    def __init__(self, filepath: str, search_engine: SearchEngine) -> None:
        super().__init__()
        self.filepath = filepath
        self.search_engine = search_engine
        self._is_cancelled = False
    
    def run(self) -> None:
        """Execute the processing pipeline."""
        try:
            self.progress.emit(10, "Parsing 23andMe file...")
            
            parser = Parser23andMe()
            
            # Set up progress callback
            def parsing_progress(lines_processed: int, total_lines: int) -> None:
                if total_lines > 0:
                    # Use 10-40% of progress for parsing
                    progress_percent = 10 + int((lines_processed / total_lines) * 30)
                    self.progress.emit(progress_percent, 
                                      f"Parsing 23andMe file... ({lines_processed}/{total_lines} lines)")
            
            parser.set_progress_callback(parsing_progress)
            snp_records = parser.parse_file(self.filepath)
            stats = parser.get_parse_stats()
            
            if self._is_cancelled:
                return
            
            self.progress.emit(50, f"Matching {len(snp_records)} SNPs against GWAS database...")
            
            matches = self.search_engine.match_user_snps(snp_records)
            
            if self._is_cancelled:
                return
            
            self.progress.emit(100, "Processing complete!")
            
            stats['matches_found'] = len(matches)
            self.finished.emit(matches, stats)
            
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
        <p>This application analyzes your raw genetic data from 23andMe and compares it against 
        the <b>GWAS Catalog</b> (Genome-Wide Association Studies Catalog) - a curated database 
        of genetic variants associated with various traits and diseases.</p>
        
        <h2>How does it work?</h2>
        <ol>
            <li><b>Upload your data:</b> Click "Upload 23andMe File" and select your raw data file 
            (usually named something like "genome_Your_Name.txt")</li>
            <li><b>Parsing:</b> The application reads your genetic variants (SNPs - Single Nucleotide Polymorphisms)</li>
            <li><b>Matching:</b> Each of your SNPs is compared against the GWAS database to find associations</li>
            <li><b>Scoring:</b> An Impact Score (0-10) is calculated for each match based on statistical significance</li>
            <li><b>Results:</b> Matches are displayed in a sortable, filterable table</li>
        </ol>
        
        <h2>Table Columns Explained</h2>
        <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%;">
            <tr style="background-color: #e0e0e0;">
                <th>Column</th>
                <th>Description</th>
            </tr>
            <tr>
                <td><b>SNP ID</b></td>
                <td>The unique identifier for the genetic variant (e.g., rs6983267). 
                The "rs" prefix stands for "Reference SNP".</td>
            </tr>
            <tr>
                <td><b>Gene</b></td>
                <td>The gene where this variant is located or nearest to. 
                Genes are segments of DNA that encode proteins.</td>
            </tr>
            <tr>
                <td><b>Trait</b></td>
                <td>The disease, condition, or characteristic associated with this variant 
                according to published GWAS studies.</td>
            </tr>
            <tr>
                <td><b>User Genotype</b></td>
                <td>YOUR specific genotype at this position. You have two copies (alleles) - 
                one from each parent. E.g., "AG" means you have one A and one G allele.</td>
            </tr>
            <tr>
                <td><b>Risk Allele</b></td>
                <td>The allele (A, T, C, or G) associated with increased risk or effect 
                for the trait. If your genotype contains this allele, it's highlighted in red.</td>
            </tr>
            <tr>
                <td><b>P-value</b></td>
                <td>Statistical significance of the association. Lower = more significant. 
                Values like 1e-10 mean 0.0000000001. Generally, p < 5e-8 is considered genome-wide significant.</td>
            </tr>
            <tr>
                <td><b>Category</b></td>
                <td>Classification of the trait: Metabolic, Cardiovascular, Neuropsychiatric, 
                Physical Trait, Oncology, Immune, Infectious, or Other.</td>
            </tr>
            <tr>
                <td><b>Impact Score</b></td>
                <td>A calculated score from 0-10 combining p-value significance and allele rarity. 
                Higher scores indicate potentially more impactful variants. 
                <br><br>Formula: Score = (P-value component × 7) + (Rarity component × 3)
                <br>- P-value component: Based on -log10(p-value)
                <br>- Rarity component: Based on how rare the variant is in the population</td>
            </tr>
            <tr>
                <td><b>Interpretation</b></td>
                <td>Human-readable interpretation of the impact score:
                <br>- Very High Impact (≥8): Strong statistical association
                <br>- High Impact (6-8): Significant association  
                <br>- Moderate Impact (4-6): Notable association
                <br>- Low Impact (2-4): Weak association
                <br>- Minimal Impact (<2): Very weak association</td>
            </tr>
            <tr>
                <td><b>Explain</b></td>
                <td>Click this button to see detailed information about that specific result, 
                including how the score was calculated and links for further research.</td>
            </tr>
        </table>
        
        <h2>Using Filters</h2>
        <ul>
            <li><b>Min Impact Score:</b> Show only results above this threshold</li>
            <li><b>Max P-value:</b> Show only results with statistical significance below this value</li>
            <li><b>Category:</b> Filter by trait category</li>
            <li><b>Search:</b> Free text search in traits, genes, and SNP IDs</li>
        </ul>
        
        <h2>Important Disclaimer</h2>
        <p style="color: #b00000;"><b>⚠️ This application is for educational and informational purposes only.</b></p>
        <p>The results should NOT be used for medical diagnosis or treatment decisions. 
        Genetic associations are complex and influenced by many factors including:</p>
        <ul>
            <li>Environment and lifestyle</li>
            <li>Gene-gene interactions</li>
            <li>Population-specific effects</li>
            <li>Incomplete scientific understanding</li>
        </ul>
        <p>Always consult with qualified healthcare professionals and genetic counselors 
        for interpretation of genetic data.</p>
        
        <h2>Data Sources</h2>
        <ul>
            <li><b>GWAS Catalog:</b> <a href="https://www.ebi.ac.uk/gwas/">https://www.ebi.ac.uk/gwas/</a></li>
            <li><b>dbSNP:</b> <a href="https://www.ncbi.nlm.nih.gov/snp/">https://www.ncbi.nlm.nih.gov/snp/</a></li>
            <li><b>gnomAD (allele frequencies):</b> <a href="https://gnomad.broadinstitute.org/">https://gnomad.broadinstitute.org/</a></li>
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
    - File upload for 23andMe data
    - Results table with sorting and pagination
    - Real-time filtering by score, p-value, category, and text search
    """
    
    def __init__(self) -> None:
        super().__init__()
        
        self.all_matches: List[GWASMatch] = []
        self.filtered_matches: List[GWASMatch] = []
        self.current_page = 0
        self.worker: Optional[ProcessingWorker] = None
        
        self.search_engine = SearchEngine(DATABASE_PATH)
        
        self._init_ui()
        self._setup_connections()
        self._verify_database()
    
    def _init_ui(self) -> None:
        """Initialize the user interface."""
        self.setWindowTitle(f"{APP_NAME} v{APP_VERSION}")
        self.setMinimumSize(1200, 800)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(10)
        main_layout.setContentsMargins(15, 15, 15, 15)
        
        # Header section
        header_layout = QHBoxLayout()
        
        title_label = QLabel(APP_NAME)
        title_label.setFont(QFont('Arial', 18, QFont.Weight.Bold))
        header_layout.addWidget(title_label)
        
        header_layout.addStretch()
        
        self.help_btn = QPushButton("❓ Help")
        self.help_btn.setMinimumWidth(80)
        self.help_btn.setMinimumHeight(40)
        self.help_btn.setFont(QFont('Arial', 11))
        header_layout.addWidget(self.help_btn)
        
        self.upload_btn = QPushButton("📁 Upload 23andMe File")
        self.upload_btn.setMinimumWidth(200)
        self.upload_btn.setMinimumHeight(40)
        self.upload_btn.setFont(QFont('Arial', 11))
        header_layout.addWidget(self.upload_btn)
        
        main_layout.addLayout(header_layout)
        
        # Progress section
        self.progress_frame = QFrame()
        self.progress_frame.setVisible(False)
        progress_layout = QVBoxLayout(self.progress_frame)
        
        self.progress_label = QLabel("Processing...")
        progress_layout.addWidget(self.progress_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        progress_layout.addWidget(self.progress_bar)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setMaximumWidth(100)
        progress_layout.addWidget(self.cancel_btn)
        
        main_layout.addWidget(self.progress_frame)
        
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
        
        main_layout.addWidget(self.filters_group)
        self.filters_group.setEnabled(False)
        
        # Results count label
        self.results_label = QLabel("No results to display")
        self.results_label.setFont(QFont('Arial', 11))
        main_layout.addWidget(self.results_label)
        
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
        
        main_layout.addWidget(self.results_table, stretch=1)
        
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
        
        main_layout.addLayout(pagination_layout)
        
        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready - Upload a 23andMe file to begin")
        
        self._apply_styles()
    
    def _apply_styles(self) -> None:
        """Apply stylesheet to the application."""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f5;
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
        self.progress_bar.setValue(0)
        
        self.worker = ProcessingWorker(filepath, self.search_engine)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_processing_finished)
        self.worker.error.connect(self._on_processing_error)
        self.worker.start()
    
    def _on_progress(self, value: int, message: str) -> None:
        """Handle progress updates from worker."""
        self.progress_bar.setValue(value)
        self.progress_label.setText(message)
        self.status_bar.showMessage(message)
    
    def _on_processing_finished(self, matches: List[GWASMatch], stats: dict) -> None:
        """Handle successful processing completion."""
        self.progress_frame.setVisible(False)
        self.upload_btn.setEnabled(True)
        self.filters_group.setEnabled(True)
        
        self.all_matches = matches
        self.current_page = 0
        
        logger.info(f"Processing complete: {stats}")
        
        self._reset_filters()
        
        self.status_bar.showMessage(
            f"Loaded {stats.get('valid_snps', 0)} SNPs | "
            f"Found {stats.get('matches_found', 0)} GWAS matches | "
            f"Skipped {stats.get('skipped_undetermined', 0)} undetermined"
        )
    
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
        event.accept()
