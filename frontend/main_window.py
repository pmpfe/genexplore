"""
Main window UI for the Genetic Analysis Application.

Implements PyQt6 interface with file upload, results table, filters, and search.
"""

from typing import List, Optional
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QTableWidget, QTableWidgetItem, QFileDialog, QLabel, QLineEdit,
    QComboBox, QSlider, QProgressBar, QStatusBar, QMessageBox,
    QHeaderView, QGroupBox, QSpinBox, QFrame, QSplitter, QApplication
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QFont, QColor

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
        self.results_table.setColumnCount(9)
        self.results_table.setHorizontalHeaderLabels([
            "SNP ID", "Gene", "Trait", "User Genotype", 
            "Risk Allele", "P-value", "Category", "Impact Score", "Interpretation"
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
        self.cancel_btn.clicked.connect(self._on_cancel_clicked)
        self.prev_btn.clicked.connect(self._on_prev_page)
        self.next_btn.clicked.connect(self._on_next_page)
        self.reset_filters_btn.clicked.connect(self._reset_filters)
        
        # Real-time filter connections
        self.score_slider.valueChanged.connect(self._on_score_changed)
        self.pvalue_slider.valueChanged.connect(self._on_pvalue_changed)
        self.category_combo.currentTextChanged.connect(self._apply_filters)
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
