"""
Polygenic analysis UI components.

Contains widgets for:
- Polygenic score browser table
- Score detail view with distribution visualization
- Database management settings
"""

import sys
from pathlib import Path
from typing import List, Optional, Dict
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox,
    QLineEdit, QProgressBar, QFrame, QGroupBox, QSplitter,
    QDialog, QTextEdit, QTextBrowser, QScrollArea, QDialogButtonBox,
    QMessageBox, QSizePolicy, QSlider, QPlainTextEdit, QFormLayout
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QUrl, QProcess
from PyQt6.QtGui import QFont, QColor, QPainter, QPen, QBrush, QDesktopServices
import gc
import re

from models.polygenic_models import (
    PolygenicScore, PolygenicResult, PopulationDistribution,
    RiskCategory, TraitCategory, DatabaseVersion
)
from models.data_models import SNPRecord
from backend.polygenic_scoring import PolygenicScorer, get_risk_interpretation, format_score_summary
from database.polygenic_database import (
    PolygenicDatabase, DatabaseVersionManager, get_gwas_database_stats
)
from config import BASE_DIR
from utils.logging_config import get_logger

logger = get_logger(__name__)


def _coverage_confidence_label(coverage_percent: float) -> tuple[str, str, str]:
    """
    Return a (label, color_hex, explanation) tuple for a coverage percentage.
    """
    if coverage_percent >= 80:
        return (
            "High",
            "#c8e6c9",  # light green background
            "High coverage (>=80%): Results are based on most variants in this score and are more reliable."
        )
    if coverage_percent >= 50:
        return (
            "Moderate",
            "#fff3cd",  # light amber
            "Moderate coverage (50–79%): Useful but interpret with caution; some variants are missing."
        )
    if coverage_percent >= 20:
        return (
            "Low",
            "#ffe0b2",  # light orange
            "Low coverage (20–49%): Limited reliability; consider reprocessing with denser data or imputation."
        )
    return (
        "Very low",
        "#f8d7da",  # light red
        "Very low coverage (<20%): Not informative for clinical interpretation. Consider using WGS/VCF or imputation."
    )


class PolygenicComputeWorker(QThread):
    """
    Background worker for computing all polygenic scores.
    
    Loads score variants incrementally from database to avoid UI freeze.
    
    Signals:
        progress: (current, total, message) progress updates
        score_computed: (pgs_id, result) individual score completed
        finished: (results) all scores computed
        error: (message) error occurred
    """
    progress = pyqtSignal(int, int, str)
    score_computed = pyqtSignal(str, object)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    
    def __init__(
        self,
        snp_records: List[SNPRecord],
        score_ids: List[str],  # Just IDs, not full scores
        distributions: Dict[str, PopulationDistribution],
        db_path: str  # Database path for loading variants in thread
    ) -> None:
        super().__init__()
        self.snp_records = snp_records
        self.score_ids = score_ids
        self.distributions = distributions
        self.db_path = db_path
        self._is_cancelled = False
    
    def run(self) -> None:
        """Execute polygenic score computation."""
        try:
            # Create database connection in this thread
            pgs_db = PolygenicDatabase(self.db_path)
            
            scorer = PolygenicScorer()
            scorer.load_genotypes(self.snp_records)
            # Release worker's reference to the full SNP list to reduce peak
            # memory usage. The main thread retains its own reference.
            try:
                self.snp_records = []
                del self.snp_records
            except Exception:
                pass
            gc.collect()

            results = []
            total = len(self.score_ids)
            
            for i, pgs_id in enumerate(self.score_ids):
                if self._is_cancelled:
                    return
                
                self.progress.emit(i + 1, total, f"Loading & computing score {i+1}/{total}...")
                
                # Load score with variants in background thread
                score = pgs_db.get_score_with_variants(pgs_id)
                if not score:
                    continue
                
                pop_dist = self.distributions.get(pgs_id)
                result = scorer.compute_score(score, pop_dist)
                results.append(result)
                
                self.score_computed.emit(pgs_id, result)
            
            self.finished.emit(results)
            
        except Exception as e:
            logger.exception(f"Error computing polygenic scores: {e}")
            self.error.emit(str(e))
    
    def cancel(self) -> None:
        """Cancel computation."""
        self._is_cancelled = True


class DistributionWidget(QWidget):
    """
    Widget to display a population distribution with user's position marked.
    """
    
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumSize(400, 150)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        
        self.distribution: Optional[PopulationDistribution] = None
        self.user_score: Optional[float] = None
        self.user_percentile: Optional[float] = None
    
    def set_data(
        self,
        distribution: PopulationDistribution,
        user_score: float,
        user_percentile: float
    ) -> None:
        """Set the distribution and user data to display."""
        self.distribution = distribution
        self.user_score = user_score
        self.user_percentile = user_percentile
        self.update()
    
    def clear_data(self) -> None:
        """Clear the display."""
        self.distribution = None
        self.user_score = None
        self.user_percentile = None
        self.update()
    
    def paintEvent(self, event) -> None:
        """Paint the distribution visualization."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        width = self.width()
        height = self.height()
        margin = 40
        graph_width = width - 2 * margin
        graph_height = height - 2 * margin
        
        # Background
        painter.fillRect(0, 0, width, height, QColor(255, 255, 255))
        
        if not self.distribution:
            painter.setPen(QPen(QColor(128, 128, 128)))
            painter.drawText(
                margin, margin, graph_width, graph_height,
                Qt.AlignmentFlag.AlignCenter,
                "No distribution data available"
            )
            return
        
        # Draw bell curve approximation
        import math
        mean = self.distribution.mean
        std = self.distribution.std
        
        if std <= 0:
            return
        
        # Calculate curve points
        points = []
        min_x = mean - 3 * std
        max_x = mean + 3 * std
        
        for i in range(101):
            x_val = min_x + (max_x - min_x) * i / 100
            z = (x_val - mean) / std
            y_val = math.exp(-0.5 * z * z) / (std * math.sqrt(2 * math.pi))
            points.append((x_val, y_val))
        
        max_y = max(p[1] for p in points)
        
        # Draw axes
        painter.setPen(QPen(QColor(100, 100, 100), 1))
        painter.drawLine(margin, height - margin, width - margin, height - margin)  # X axis
        painter.drawLine(margin, margin, margin, height - margin)  # Y axis
        
        # Draw curve
        painter.setPen(QPen(QColor(70, 130, 180), 2))
        for i in range(len(points) - 1):
            x1 = margin + (points[i][0] - min_x) / (max_x - min_x) * graph_width
            y1 = height - margin - (points[i][1] / max_y) * graph_height
            x2 = margin + (points[i + 1][0] - min_x) / (max_x - min_x) * graph_width
            y2 = height - margin - (points[i + 1][1] / max_y) * graph_height
            painter.drawLine(int(x1), int(y1), int(x2), int(y2))
        
        # Draw user position
        if self.user_score is not None:
            user_x = margin + (self.user_score - min_x) / (max_x - min_x) * graph_width
            user_x = max(margin, min(width - margin, user_x))
            
            # Vertical line for user position
            painter.setPen(QPen(QColor(220, 50, 50), 2, Qt.PenStyle.DashLine))
            painter.drawLine(int(user_x), margin, int(user_x), height - margin)
            
            # User marker
            painter.setBrush(QBrush(QColor(220, 50, 50)))
            painter.drawEllipse(int(user_x) - 6, height - margin - 6, 12, 12)
            
            # Label
            painter.setPen(QPen(QColor(0, 0, 0)))
            painter.setFont(QFont('Arial', 9, QFont.Weight.Bold))
            label = f"You: {self.user_percentile:.0f}th percentile"
            painter.drawText(int(user_x) - 50, margin - 5, label)
        
        # Draw percentile labels
        painter.setPen(QPen(QColor(100, 100, 100)))
        painter.setFont(QFont('Arial', 8))
        for pct in [5, 25, 50, 75, 95]:
            if pct in self.distribution.percentiles:
                pct_x = margin + (self.distribution.percentiles[pct] - min_x) / (max_x - min_x) * graph_width
                if margin < pct_x < width - margin:
                    painter.drawText(int(pct_x) - 10, height - margin + 15, f"{pct}%")


class ScoreDetailDialog(QDialog):
    """
    Dialog showing detailed analysis of a polygenic score result.
    """
    
    def __init__(
        self,
        result: PolygenicResult,
        score: PolygenicScore,
        distribution: Optional[PopulationDistribution],
        parent=None
    ) -> None:
        super().__init__(parent)
        self.result = result
        self.score = score
        self.distribution = distribution
        
        self.setWindowTitle(f"Polygenic Score Details: {result.trait_name}")
        self.setMinimumSize(900, 750)
        self._init_ui()
    
    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        
        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(8)
        
        # Row 1: Your Results (left) + Population Distribution (right)
        row1_layout = QHBoxLayout()
        row1_layout.setSpacing(10)
        
        # Your Results group (left)
        results_group = QGroupBox("Your Results")
        results_layout = QVBoxLayout(results_group)
        results_layout.setContentsMargins(8, 12, 8, 8)
        results_layout.setSpacing(4)
        
        # Results shown as individual labeled fields with tooltips
        results_form = QFormLayout()
        results_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        pct_label = QLabel(f"{self.result.percentile:.1f}th")
        pct_label.setToolTip("Percentile: estimated position in the reference population (0–100). Higher percentiles indicate greater genetic predisposition for the trait.")
        results_form.addRow(QLabel("<b>Percentile:</b>"), pct_label)

        risk_color = {
            RiskCategory.LOW: "#006000",
            RiskCategory.INTERMEDIATE: "#806000",
            RiskCategory.HIGH: "#b00000"
        }.get(self.result.risk_category, "#000000")
        risk_label = QLabel(self.result.risk_category.value)
        risk_label.setStyleSheet(f"color: {risk_color}; font-weight: bold;")
        risk_label.setToolTip("Risk category: classification based on percentile thresholds (Low/Intermediate/High).")
        results_form.addRow(QLabel("<b>Risk Category:</b>"), risk_label)

        raw_label = QLabel(f"{self.result.raw_score:.4f}")
        raw_label.setToolTip("Raw score: weighted sum of variant contributions for this individual.")
        results_form.addRow(QLabel("<b>Raw Score:</b>"), raw_label)

        z_label = QLabel(f"{self.result.normalized_score:.2f}")
        z_label.setToolTip("Z-score: raw score normalized to the reference population (z = (raw - mean)/sd).")
        results_form.addRow(QLabel("<b>Z-Score:</b>"), z_label)

        pop_label = QLabel(self.result.population_reference)
        pop_label.setToolTip("Reference population used for normalization/interpolation.")
        results_form.addRow(QLabel("<b>Population:</b>"), pop_label)

        results_widget = QWidget()
        results_widget.setLayout(results_form)
        results_widget.setMinimumHeight(150)
        results_widget.setMaximumHeight(180)
        results_layout.addWidget(results_widget)
        row1_layout.addWidget(results_group, stretch=1)
        
        # Population Distribution group (right)
        dist_group = QGroupBox("Population Distribution")
        dist_layout = QVBoxLayout(dist_group)
        dist_layout.setContentsMargins(8, 12, 8, 8)
        
        self.dist_widget = DistributionWidget()
        self.dist_widget.setMinimumHeight(150)
        self.dist_widget.setMaximumHeight(180)
        if self.distribution:
            self.dist_widget.set_data(
                self.distribution,
                self.result.raw_score,
                self.result.percentile
            )
        dist_layout.addWidget(self.dist_widget)
        row1_layout.addWidget(dist_group, stretch=1)
        
        content_layout.addLayout(row1_layout)
        
        # Row 2: Scientific Context (left) + Quality Information (right)
        row2_layout = QHBoxLayout()
        row2_layout.setSpacing(10)
        
        # Scientific context (left)
        context_group = QGroupBox("Scientific Context")
        context_layout = QVBoxLayout(context_group)
        context_layout.setContentsMargins(8, 12, 8, 8)
        context_layout.setSpacing(2)
        
        # Scientific context shown as labeled fields with tooltips
        context_form = QFormLayout()

        id_label = QLabel(self.score.pgs_id)
        id_label.setToolTip("PGS identifier in the PGS Catalog.")
        context_form.addRow(QLabel("<b>Score ID:</b>"), id_label)

        year_label = QLabel(str(self.score.publication_year or 'N/A'))
        year_label.setToolTip("Year of publication/source for this PGS.")
        context_form.addRow(QLabel("<b>Year:</b>"), year_label)

        study_pop_label = QLabel(self.score.study_population or 'N/A')
        study_pop_label.setToolTip("Population used to develop the score (may affect transferability).")
        context_form.addRow(QLabel("<b>Population:</b>"), study_pop_label)

        sample_size_label = QLabel(f"{self.score.sample_size:,}" if self.score.sample_size else 'N/A')
        sample_size_label.setToolTip("Sample size of the discovery study used to derive the PGS.")
        context_form.addRow(QLabel("<b>Sample Size:</b>"), sample_size_label)

        num_variants_label = QLabel(f"{self.score.num_variants:,}")
        num_variants_label.setToolTip("Total number of variants in the score definition.")
        context_form.addRow(QLabel("<b>Variants:</b>"), num_variants_label)

        links_label = QLabel()
        links_html = []
        if self.score.publication_doi:
            links_html.append(f'<a href="https://doi.org/{self.score.publication_doi}">DOI</a>')
        links_html.append(f'<a href="https://www.pgscatalog.org/score/{self.score.pgs_id}/">PGS Catalog</a>')
        trait_search = self.score.trait_name.replace(' ', '+')
        links_html.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/?term={trait_search}+polygenic">PubMed</a>')
        links_label.setText(' | '.join(links_html))
        links_label.setOpenExternalLinks(True)
        links_label.setToolTip('Links to external resources for this score (DOI, PGS Catalog, PubMed).')
        context_form.addRow(QLabel("<b>Links:</b>"), links_label)

        context_widget = QWidget()
        context_widget.setLayout(context_form)
        context_widget.setMinimumHeight(120)
        context_widget.setMaximumHeight(150)
        context_layout.addWidget(context_widget)
        row2_layout.addWidget(context_group, stretch=1)
        
        # Quality information (right)
        quality_group = QGroupBox("Quality Information")
        quality_layout = QVBoxLayout(quality_group)
        quality_layout.setContentsMargins(8, 12, 8, 8)
        quality_layout.setSpacing(2)
        
        # Quality information shown with per-field tooltips
        quality_form = QFormLayout()

        cov_label = QLabel(f"{self.result.coverage_percent:.0f}%")
        cov_label.setToolTip(
            "Coverage = percentage of variants in the PGS that were found in your genotype file (variants_found / variants_total). Low coverage reduces reliability."
        )
        quality_form.addRow(QLabel("<b>Coverage:</b>"), cov_label)

        found_label = QLabel(f"{self.result.variants_found:,}")
        found_label.setToolTip("Number of variants from this PGS that were matched in your genotype file.")
        quality_form.addRow(QLabel("<b>Variants Found:</b>"), found_label)

        total_label = QLabel(f"{self.result.variants_total:,}")
        total_label.setToolTip("Total number of variants defined for this PGS.")
        quality_form.addRow(QLabel("<b>Variants Total:</b>"), total_label)

        warning_label = QLabel("Yes" if self.result.is_low_coverage() else "No")
        warning_label.setToolTip("Indicates whether coverage is below the recommended threshold (80%).")
        quality_form.addRow(QLabel("<b>Low Coverage Warning:</b>"), warning_label)

        quality_widget = QWidget()
        quality_widget.setLayout(quality_form)
        quality_widget.setMinimumHeight(120)
        quality_widget.setMaximumHeight(150)
        quality_layout.addWidget(quality_widget)
        row2_layout.addWidget(quality_group, stretch=1)
        
        content_layout.addLayout(row2_layout)
        
        # Top contributors - takes remaining space
        if self.result.variant_contributions:
            contrib_group = QGroupBox("Top Contributing Variants")
            contrib_layout = QVBoxLayout(contrib_group)
            contrib_layout.setContentsMargins(8, 12, 8, 8)
            
            contrib_table = QTableWidget()
            contrib_table.setColumnCount(5)
            contrib_table.setHorizontalHeaderLabels([
                "Variant", "Contribution", "Direction", "Gene Info", "Resources"
            ])
            
            header = contrib_table.horizontalHeader()
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
            header.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
            
            top_contribs = self.result.get_top_contributors(15)
            contrib_table.setRowCount(len(top_contribs))

            # Add header tooltips explaining each column
            headers = ["Variant (rsID)", "Absolute contribution to score", "Direction of effect", "Quick gene lookup", "External resources"]
            for hi in range(contrib_table.columnCount()):
                hi_item = contrib_table.horizontalHeaderItem(hi)
                if hi_item:
                    hi_item.setToolTip(headers[hi])
            
            for row, (rsid, contrib) in enumerate(top_contribs):
                # Variant ID
                rsid_item = QTableWidgetItem(rsid)
                rsid_item.setFlags(rsid_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                contrib_table.setItem(row, 0, rsid_item)
                
                # Contribution value
                contrib_item = QTableWidgetItem(f"{abs(contrib):.4f}")
                contrib_item.setFlags(contrib_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                contrib_table.setItem(row, 1, contrib_item)
                
                # Direction
                direction = "↑ Increases" if contrib > 0 else "↓ Decreases"
                dir_item = QTableWidgetItem(direction)
                dir_item.setFlags(dir_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                dir_item.setForeground(QColor(180, 0, 0) if contrib > 0 else QColor(0, 100, 0))
                contrib_table.setItem(row, 2, dir_item)
                
                # Gene Info button (links to multiple resources)
                gene_btn = QPushButton("🧬 Gene Info")
                gene_btn.setStyleSheet("""
                    QPushButton { background-color: #5cb85c; color: white; border: none; 
                                  padding: 3px 6px; border-radius: 3px; font-size: 10px; }
                    QPushButton:hover { background-color: #449d44; }
                """)
                gene_btn.clicked.connect(lambda _, r=rsid: self._open_gene_info(r))
                contrib_table.setCellWidget(row, 3, gene_btn)
                
                # Resources button (links to dbSNP, ClinVar, etc.)
                resources_widget = QWidget()
                resources_layout = QHBoxLayout(resources_widget)
                resources_layout.setContentsMargins(2, 2, 2, 2)
                resources_layout.setSpacing(4)
                
                dbsnp_btn = QPushButton("dbSNP")
                dbsnp_btn.setStyleSheet("""
                    QPushButton { background-color: #337ab7; color: white; border: none; 
                                  padding: 2px 5px; border-radius: 2px; font-size: 9px; }
                    QPushButton:hover { background-color: #286090; }
                """)
                dbsnp_btn.clicked.connect(lambda _, r=rsid: self._open_dbsnp(r))
                resources_layout.addWidget(dbsnp_btn)
                
                clinvar_btn = QPushButton("ClinVar")
                clinvar_btn.setStyleSheet("""
                    QPushButton { background-color: #d9534f; color: white; border: none; 
                                  padding: 2px 5px; border-radius: 2px; font-size: 9px; }
                    QPushButton:hover { background-color: #c9302c; }
                """)
                clinvar_btn.clicked.connect(lambda _, r=rsid: self._open_clinvar(r))
                resources_layout.addWidget(clinvar_btn)
                
                ensembl_btn = QPushButton("Ensembl")
                ensembl_btn.setStyleSheet("""
                    QPushButton { background-color: #f0ad4e; color: white; border: none; 
                                  padding: 2px 5px; border-radius: 2px; font-size: 9px; }
                    QPushButton:hover { background-color: #ec971f; }
                """)
                ensembl_btn.clicked.connect(lambda _, r=rsid: self._open_ensembl(r))
                resources_layout.addWidget(ensembl_btn)
                
                resources_layout.addStretch()
                contrib_table.setCellWidget(row, 4, resources_widget)
            
            contrib_table.setMinimumHeight(200)
            contrib_layout.addWidget(contrib_table)
            content_layout.addWidget(contrib_group, stretch=1)
        
        scroll.setWidget(content)
        layout.addWidget(scroll, stretch=1)
        
        # Disclaimer
        disclaimer = QLabel(
            "⚠️ Disclaimer: This is a research tool for educational purposes only. "
            "Results should NOT be used for medical diagnosis or treatment decisions."
        )
        disclaimer.setWordWrap(True)
        disclaimer.setStyleSheet("color: #b00000; font-weight: bold; padding: 5px;")
        layout.addWidget(disclaimer)
        
        # Close button
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)
    
    def _open_gene_info(self, rsid: str) -> None:
        """Open GeneCards search for the variant."""
        url = f"https://www.genecards.org/Search/Keyword?queryString={rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _open_dbsnp(self, rsid: str) -> None:
        """Open dbSNP page for the variant."""
        url = f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _open_clinvar(self, rsid: str) -> None:
        """Open ClinVar search for the variant."""
        url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _open_ensembl(self, rsid: str) -> None:
        """Open Ensembl page for the variant."""
        url = f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={rsid}"
        QDesktopServices.openUrl(QUrl(url))
    
    def _get_results_html(self) -> str:
        r = self.result
        risk_color = {
            RiskCategory.LOW: "#006000",
            RiskCategory.INTERMEDIATE: "#806000",
            RiskCategory.HIGH: "#b00000"
        }.get(r.risk_category, "#000000")
        
        return f"""
        <style>td {{ padding: 2px 5px; }} p {{ margin: 5px 0; }}</style>
        <table border="0" style="width: 100%;">
            <tr>
                <td><b>Percentile:</b></td>
                <td><span style="font-size: 14pt; font-weight: bold;">{r.percentile:.1f}th</span></td>
            </tr>
            <tr>
                <td><b>Risk Category:</b></td>
                <td><span style="font-size: 12pt; font-weight: bold; color: {risk_color};">{r.risk_category.value}</span></td>
            </tr>
            <tr>
                <td><b>Raw Score:</b></td>
                <td>{r.raw_score:.4f}</td>
            </tr>
            <tr>
                <td><b>Z-Score:</b></td>
                <td>{r.normalized_score:.2f}</td>
            </tr>
            <tr>
                <td><b>Population:</b></td>
                <td>{r.population_reference}</td>
            </tr>
        </table>
        """
    
    def _get_context_html(self) -> str:
        s = self.score
        
        # Build links
        links = []
        if s.publication_doi:
            links.append(f'<a href="https://doi.org/{s.publication_doi}">DOI: {s.publication_doi}</a>')
        
        # Add PGS Catalog link
        links.append(f'<a href="https://www.pgscatalog.org/score/{s.pgs_id}/">PGS Catalog</a>')
        
        # Add PubMed search link
        trait_search = s.trait_name.replace(" ", "+")
        links.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/?term={trait_search}+polygenic">PubMed</a>')
        
        links_html = " | ".join(links)
        
        sample_size_str = f"{s.sample_size:,}" if s.sample_size else "N/A"
        num_variants_str = f"{s.num_variants:,}" if s.num_variants else "N/A"
        
        return f"""
        <style>td {{ padding: 1px 4px; }} p {{ margin: 3px 0; }}</style>
        <table border="0" style="width: 100%;">
            <tr><td><b>Score ID:</b></td><td>{s.pgs_id}</td></tr>
            <tr><td><b>Year:</b></td><td>{s.publication_year or 'N/A'}</td></tr>
            <tr><td><b>Population:</b></td><td>{s.study_population or 'N/A'}</td></tr>
            <tr><td><b>Sample Size:</b></td><td>{sample_size_str}</td></tr>
            <tr><td><b>Variants:</b></td><td>{num_variants_str}</td></tr>
        </table>
        <p><b>Links:</b> {links_html}</p>
        """
    
    def _get_quality_html(self) -> str:
        r = self.result
        coverage_color = "#006000" if r.coverage_percent >= 80 else "#b00000"
        
        warning = ""
        if r.is_low_coverage():
            warning = '<p style="color: #b00000; margin: 3px 0;"><b>⚠️ Low coverage</b> - accuracy may be reduced</p>'
        
        pop_warning = ""
        if r.population_reference == "EUR":
            pop_warning = '<p style="color: #666; margin: 3px 0; font-size: 10px;">Note: Score optimized for European ancestry</p>'
        
        return f"""
        <style>td {{ padding: 1px 4px; }}</style>
        <table border="0" style="width: 100%;">
            <tr><td><b>Variants Found:</b></td><td>{r.variants_found} / {r.variants_total}</td></tr>
            <tr><td><b>Coverage:</b></td><td><span style="color: {coverage_color}; font-weight: bold;">{r.coverage_percent:.0f}%</span></td></tr>
        </table>
        {warning}{pop_warning}
        """


class PolygenicBrowserWidget(QWidget):
    """
    Widget for browsing and computing polygenic scores.
    """
    
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        
        self.pgs_db = PolygenicDatabase()
        self.version_manager = DatabaseVersionManager()
        self.scores: List[PolygenicScore] = []
        self.results: Dict[str, PolygenicResult] = {}
        self.distributions: Dict[str, PopulationDistribution] = {}
        self.snp_records: List[SNPRecord] = []
        self.worker: Optional[PolygenicComputeWorker] = None
        
        self._init_ui()
        self._load_scores()
    
    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Info banner - dynamic based on database content
        self.info_frame = QFrame()
        self.info_frame.setStyleSheet("""
            QFrame { background-color: #e7f3fe; border: 1px solid #2196F3; border-radius: 4px; padding: 5px; }
            QLabel { color: #1565C0; }
        """)
        info_layout = QHBoxLayout(self.info_frame)
        info_layout.setContentsMargins(10, 5, 10, 5)
        
        info_icon = QLabel("ℹ️")
        info_layout.addWidget(info_icon)
        
        self.info_text = QLabel()
        self.info_text.setWordWrap(True)
        self.info_text.setFont(QFont('Arial', 9))
        info_layout.addWidget(self.info_text, stretch=1)
        
        layout.addWidget(self.info_frame)
        
        # Status bar
        status_layout = QHBoxLayout()
        
        self.status_label = QLabel("Load genetic data to compute polygenic scores")
        self.status_label.setFont(QFont('Arial', 10))
        status_layout.addWidget(self.status_label)
        
        status_layout.addStretch()
        
        self.compute_btn = QPushButton("🧬 Compute All Scores")
        self.compute_btn.setEnabled(False)
        self.compute_btn.clicked.connect(self._start_computation)
        status_layout.addWidget(self.compute_btn)
        
        layout.addLayout(status_layout)
        
        # Progress bar (hidden by default)
        self.progress_frame = QFrame()
        self.progress_frame.setVisible(False)
        progress_layout = QHBoxLayout(self.progress_frame)
        progress_layout.setContentsMargins(0, 5, 0, 5)
        
        self.progress_label = QLabel("Computing...")
        progress_layout.addWidget(self.progress_label)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        progress_layout.addWidget(self.progress_bar, stretch=1)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self._cancel_computation)
        progress_layout.addWidget(self.cancel_btn)
        
        layout.addWidget(self.progress_frame)
        
        # Filters
        filter_layout = QHBoxLayout()
        
        filter_layout.addWidget(QLabel("Category:"))
        self.category_filter = QComboBox()
        self.category_filter.addItem("All Categories")
        for cat in TraitCategory:
            self.category_filter.addItem(cat.value)
        self.category_filter.currentTextChanged.connect(self._apply_filters)
        filter_layout.addWidget(self.category_filter)
        
        filter_layout.addWidget(QLabel("Search:"))
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search traits...")
        self.search_input.setClearButtonEnabled(True)
        self.search_input.textChanged.connect(self._on_search_changed)
        filter_layout.addWidget(self.search_input, stretch=1)
        
        layout.addLayout(filter_layout)
        
        # Coverage filter slider
        coverage_layout = QHBoxLayout()
        coverage_layout.addWidget(QLabel("Min. Coverage:"))
        
        self.coverage_slider = QSlider(Qt.Orientation.Horizontal)
        self.coverage_slider.setRange(0, 100)
        self.coverage_slider.setValue(0)
        self.coverage_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.coverage_slider.setTickInterval(10)
        self.coverage_slider.setMaximumWidth(200)
        self.coverage_slider.valueChanged.connect(self._on_coverage_changed)
        coverage_layout.addWidget(self.coverage_slider)
        
        self.coverage_label = QLabel("0%")
        self.coverage_label.setMinimumWidth(40)
        coverage_layout.addWidget(self.coverage_label)
        
        coverage_layout.addStretch()
        layout.addLayout(coverage_layout)
        
        # Results table
        self.table = QTableWidget()
        self.table.setColumnCount(8)
        self.table.setHorizontalHeaderLabels([
            "Trait", "Category", "Variants", "Cover%", "Score", "Percentile", "Risk", "Details"
        ])
        
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(6, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(7, QHeaderView.ResizeMode.ResizeToContents)
        
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        
        layout.addWidget(self.table)
        
        # Search debounce timer
        self._search_timer = QTimer()
        self._search_timer.setSingleShot(True)
        self._search_timer.timeout.connect(self._apply_filters)
        
        self._apply_styles()
    
    def _apply_styles(self) -> None:
        self.setStyleSheet("""
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
        """)
    
    def _load_scores(self) -> None:
        """Load polygenic scores from database."""
        try:
            self.scores = self.pgs_db.get_all_scores()
            self.distributions = self.pgs_db.get_all_distributions()
            self._update_table()
            self._update_info_banner()
            logger.info(f"Loaded {len(self.scores)} polygenic scores")
        except Exception as e:
            logger.error(f"Error loading polygenic scores: {e}")
            self._update_info_banner()
    
    def _update_info_banner(self) -> None:
        """Update the info banner based on database status."""
        num_scores = len(self.scores)
        is_sample_db = self.version_manager.pgs_db.is_sample_pgs_database()
        
        if num_scores == 0:
            self.info_text.setText(
                "No PGS scores are available. The app ships with example scores only when the bundled "
                "database is present; otherwise run 'python database/update_databases.py --pgs' "
                "from the update tab or terminal to download the real PGS Catalog."
            )
            self.info_frame.setStyleSheet("""
                QFrame { background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; padding: 5px; }
                QLabel { color: #856404; }
            """)
        elif is_sample_db:
            self.info_text.setText(
                f"Using the bundled example scores ({num_scores} scores). "
                "Run 'python database/update_databases.py --pgs' to download the real PGS Catalog."
            )
            self.info_frame.setStyleSheet("""
                QFrame { background-color: #e7f3fe; border: 1px solid #2196F3; border-radius: 4px; padding: 5px; }
                QLabel { color: #1565C0; }
            """)
        else:
            self.info_text.setText(
                f"Loaded {num_scores} polygenic scores from the PGS Catalog database."
            )
            self.info_frame.setStyleSheet("""
                QFrame { background-color: #d4edda; border: 1px solid #28a745; border-radius: 4px; padding: 5px; }
                QLabel { color: #155724; }
            """)
    
    def set_genotype_data(self, snp_records: List[SNPRecord]) -> None:
        """
        Set the user's genotype data for computation.
        
        Args:
            snp_records: List of SNP records from user's file.
        """
        self.snp_records = snp_records
        self.compute_btn.setEnabled(True)
        self.status_label.setText(f"{len(snp_records)} SNPs loaded - Ready to compute scores")
    
    def display_loaded_results(self, results: List[PolygenicResult]) -> None:
        """
        Display previously loaded/saved polygenic results.
        
        Used when loading a saved session to show results without recomputing.
        
        Args:
            results: List of PolygenicResult objects from saved session.
        """
        # Store results in dictionary by pgs_id
        self.results = {r.pgs_id: r for r in results}
        
        # Update table to show loaded results
        self._update_table()
        
        # Update status
        self.status_label.setText(f"Loaded {len(results)} polygenic scores from session")
        self.compute_btn.setEnabled(True)
        # Re-apply any active filters (e.g., coverage slider) after loading
        self._apply_filters()
    
    def _start_computation(self) -> None:
        """Start computing all polygenic scores."""
        if not self.snp_records:
            return
        
        # Just get score IDs - variants will be loaded in background thread
        score_ids = [score.pgs_id for score in self.scores]
        
        if not score_ids:
            self.status_label.setText("No scores available to compute")
            return
        
        self.compute_btn.setEnabled(False)
        self.progress_frame.setVisible(True)
        self.progress_bar.setValue(0)
        self.progress_label.setText(f"Starting computation of {len(score_ids)} scores...")
        
        self.worker = PolygenicComputeWorker(
            self.snp_records,
            score_ids,
            self.distributions,
            self.pgs_db.db_path  # Pass db path for thread-safe access
        )
        self.worker.progress.connect(self._on_progress)
        self.worker.score_computed.connect(self._on_score_computed)
        self.worker.finished.connect(self._on_computation_finished)
        self.worker.error.connect(self._on_computation_error)
        self.worker.start()
    
    def _cancel_computation(self) -> None:
        """Cancel ongoing computation."""
        if self.worker and self.worker.isRunning():
            self.worker.cancel()
            self.worker.wait()
        
        self.progress_frame.setVisible(False)
        self.compute_btn.setEnabled(True)
        self.status_label.setText("Computation cancelled")
    
    def _on_progress(self, current: int, total: int, message: str) -> None:
        """Handle progress updates."""
        self.progress_bar.setValue(int(current / total * 100))
        self.progress_label.setText(message)
    
    def _on_score_computed(self, pgs_id: str, result: PolygenicResult) -> None:
        """Handle individual score completion."""
        self.results[pgs_id] = result
        self._update_table_row(pgs_id, result)
    
    def _on_computation_finished(self, results: List[PolygenicResult]) -> None:
        """Handle computation completion."""
        self.progress_frame.setVisible(False)
        self.compute_btn.setEnabled(True)
        self.status_label.setText(f"Computed {len(results)} polygenic scores")
        self._update_table()
    
    def _on_computation_error(self, error_message: str) -> None:
        """Handle computation error."""
        self.progress_frame.setVisible(False)
        self.compute_btn.setEnabled(True)
        self.status_label.setText(f"Error: {error_message}")
        QMessageBox.critical(self, "Computation Error", error_message)
    
    def _on_search_changed(self) -> None:
        """Handle search input change with debounce."""
        self._search_timer.start(300)
    
    def _on_coverage_changed(self, value: int) -> None:
        """Handle coverage slider change."""
        self.coverage_label.setText(f"{value}%")
        logger.debug("Coverage slider changed: %s", value)
        self._apply_filters()
    
    def _apply_filters(self) -> None:
        """Apply current filters to the table."""
        category = self.category_filter.currentText()
        search = self.search_input.text().lower().strip()
        min_coverage = self.coverage_slider.value()
        
        for row in range(self.table.rowCount()):
            show = True
            
            # Category filter
            if category != "All Categories":
                row_category = self.table.item(row, 1)
                if row_category and row_category.text() != category:
                    show = False
            
            # Search filter
            if search and show:
                trait = self.table.item(row, 0)
                if trait and search not in trait.text().lower():
                    show = False

            # Coverage filter - apply regardless of search presence
            if min_coverage > 0 and show:
                coverage_item = self.table.item(row, 3)
                coverage_val = None

                if coverage_item:
                    # Try EditRole first - handle QVariant-like values robustly
                    try:
                        v = coverage_item.data(Qt.ItemDataRole.EditRole)
                    except Exception:
                        v = None

                    if v is not None:
                        try:
                            # v may be int/float or a QVariant-like object
                            if isinstance(v, (int, float)):
                                coverage_val = float(v)
                            else:
                                coverage_val = float(str(v))
                        except Exception:
                            coverage_val = None

                        # Some code may store fractions (0..1); convert to percentage
                        if coverage_val is not None and 0 <= coverage_val <= 1:
                            coverage_val = coverage_val * 100.0

                        # Treat negative placeholders as missing
                        if coverage_val is not None and coverage_val < 0:
                            coverage_val = None

                    # Fallback: parse displayed text like '12% (Moderate)'
                    if coverage_val is None:
                        try:
                            m = re.search(r'([0-9]+(?:\.[0-9]+)?)%', coverage_item.text())
                            if m:
                                coverage_val = float(m.group(1))
                        except Exception:
                            coverage_val = None

                if coverage_val is None:
                    # Coverage not computed yet or malformed -> hide when filtering
                    show = False
                else:
                    if coverage_val < float(min_coverage):
                        show = False
            
            self.table.setRowHidden(row, not show)
    
    def _update_table(self) -> None:
        """Update the entire table."""
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(self.scores))
        
        for row, score in enumerate(self.scores):
            result = self.results.get(score.pgs_id)
            self._set_table_row(row, score, result)
        
        self.table.setSortingEnabled(True)
    
    def _update_table_row(self, pgs_id: str, result: PolygenicResult) -> None:
        """Update a single table row with computed result."""
        for row in range(self.table.rowCount()):
            # Find the row for this pgs_id
            item = self.table.item(row, 0)
            if item and item.data(Qt.ItemDataRole.UserRole) == pgs_id:
                score = next((s for s in self.scores if s.pgs_id == pgs_id), None)
                if score:
                    self._set_table_row(row, score, result)
                break
    
    def _set_table_row(
        self,
        row: int,
        score: PolygenicScore,
        result: Optional[PolygenicResult]
    ) -> None:
        """Set data for a table row."""
        # Trait name (with pgs_id stored in user data)
        trait_item = QTableWidgetItem(score.trait_name)
        trait_item.setData(Qt.ItemDataRole.UserRole, score.pgs_id)
        trait_item.setFlags(trait_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.table.setItem(row, 0, trait_item)
        
        # Category
        cat_item = QTableWidgetItem(score.trait_category.value)
        cat_item.setFlags(cat_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.table.setItem(row, 1, cat_item)
        
        # Variants (show found/total when available)
        var_item = QTableWidgetItem(str(score.num_variants))
        var_item.setFlags(var_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        # Provide a numeric sort key (total variants) for column sorting
        var_item.setData(Qt.ItemDataRole.EditRole, int(score.num_variants or 0))
        self.table.setItem(row, 2, var_item)
        
        if result:
            # Variants found/total (display) and numeric sort key (found)
            var_item.setText(f"{result.variants_found}/{result.variants_total}")
            var_item.setToolTip(
                f"Variants matched in genotype: {result.variants_found} of {result.variants_total}."
            )
            var_item.setData(Qt.ItemDataRole.EditRole, int(result.variants_found))

            # Coverage with badge label and tooltip
            label, color_hex, explanation = _coverage_confidence_label(result.coverage_percent)
            coverage_text = f"{result.coverage_percent:.0f}% ({label})"
            coverage_item = QTableWidgetItem(coverage_text)
            coverage_item.setFlags(coverage_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            # Numeric value for sorting
            coverage_item.setData(Qt.ItemDataRole.EditRole, float(result.coverage_percent))
            # Set background color for badge feel
            coverage_item.setBackground(QColor(color_hex))
            coverage_item.setToolTip(
                f"Coverage: {result.variants_found}/{result.variants_total} = {result.coverage_percent:.0f}%\n\n{explanation}"
            )
            self.table.setItem(row, 3, coverage_item)
            
            # Score
            score_item = QTableWidgetItem(f"{result.raw_score:.3f}")
            score_item.setFlags(score_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            score_item.setData(Qt.ItemDataRole.EditRole, float(result.raw_score))
            self.table.setItem(row, 4, score_item)
            
            # Percentile
            pct_item = QTableWidgetItem(f"{result.percentile:.0f}%")
            pct_item.setFlags(pct_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            pct_item.setData(Qt.ItemDataRole.EditRole, float(result.percentile))
            self.table.setItem(row, 5, pct_item)
            
            # Risk category with color
            risk_item = QTableWidgetItem(result.risk_category.value)
            risk_item.setFlags(risk_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            # Provide numeric ranking for risk categories to enable proper sorting
            risk_rank_map = {
                RiskCategory.HIGH: 3,
                RiskCategory.INTERMEDIATE: 2,
                RiskCategory.LOW: 1,
            }
            risk_item.setData(Qt.ItemDataRole.EditRole, int(risk_rank_map.get(result.risk_category, 0)))
            
            if result.risk_category == RiskCategory.HIGH:
                risk_item.setBackground(QColor(255, 200, 200))
                risk_item.setForeground(QColor(139, 0, 0))
            elif result.risk_category == RiskCategory.LOW:
                risk_item.setBackground(QColor(200, 255, 200))
                risk_item.setForeground(QColor(0, 100, 0))
            else:
                risk_item.setBackground(QColor(255, 255, 200))
            
            self.table.setItem(row, 6, risk_item)
            
            # Details button
            details_btn = QPushButton("📊 View")
            details_btn.setStyleSheet("""
                QPushButton {
                    background-color: #5cb85c;
                    color: white;
                    border: none;
                    padding: 4px 8px;
                    border-radius: 3px;
                }
                QPushButton:hover {
                    background-color: #449d44;
                }
            """)
            details_btn.clicked.connect(
                lambda checked, pid=score.pgs_id: self._show_details(pid)
            )
            self.table.setCellWidget(row, 7, details_btn)
        else:
            # Not computed yet
            for col in range(3, 7):
                item = QTableWidgetItem("-")
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                item.setForeground(QColor(150, 150, 150))
                # Provide a numeric placeholder so sorting remains stable
                item.setData(Qt.ItemDataRole.EditRole, -1)
                self.table.setItem(row, col, item)
            
            self.table.setCellWidget(row, 7, None)
    
    def _show_details(self, pgs_id: str) -> None:
        """Show detailed view for a score."""
        result = self.results.get(pgs_id)
        if not result:
            return
        
        score = self.pgs_db.get_score_with_variants(pgs_id)
        if not score:
            return
        
        # Re-compute with contributions tracking
        scorer = PolygenicScorer()
        scorer.load_genotypes(self.snp_records)
        distribution = self.distributions.get(pgs_id)
        result = scorer.compute_score(score, distribution, track_contributions=True)
        
        dialog = ScoreDetailDialog(result, score, distribution, self)
        dialog.exec()


class DatabaseSettingsWidget(QWidget):
    """
    Widget for managing database versions and updates.
    """
    
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        
        self.version_manager = DatabaseVersionManager()
        self.update_process: Optional[QProcess] = None
        self.update_output: Optional[QPlainTextEdit] = None
        self._init_ui()
        self._refresh_versions()
    
    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)
        
        # GWAS Database section
        gwas_group = QGroupBox("GWAS Catalog Database")
        gwas_layout = QVBoxLayout(gwas_group)
        
        self.gwas_version_label = QLabel("Version: Loading...")
        gwas_layout.addWidget(self.gwas_version_label)
        
        self.gwas_date_label = QLabel("Last Updated: Loading...")
        gwas_layout.addWidget(self.gwas_date_label)
        
        self.gwas_count_label = QLabel("Records: Loading...")
        gwas_layout.addWidget(self.gwas_count_label)
        
        gwas_link_layout = QHBoxLayout()
        gwas_link = QPushButton("🔗 Visit GWAS Catalog")
        gwas_link.setStyleSheet("QPushButton { text-align: left; color: #337ab7; background: none; border: none; }")
        gwas_link.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://www.ebi.ac.uk/gwas/")))
        gwas_link_layout.addWidget(gwas_link)
        gwas_link_layout.addStretch()
        gwas_layout.addLayout(gwas_link_layout)
        
        layout.addWidget(gwas_group)
        
        # PGS Database section
        pgs_group = QGroupBox("Polygenic Score Catalog")
        pgs_layout = QVBoxLayout(pgs_group)
        
        self.pgs_version_label = QLabel("Version: Loading...")
        pgs_layout.addWidget(self.pgs_version_label)
        
        self.pgs_date_label = QLabel("Last Updated: Loading...")
        pgs_layout.addWidget(self.pgs_date_label)
        
        self.pgs_count_label = QLabel("Scores: Loading...")
        pgs_layout.addWidget(self.pgs_count_label)
        
        self.pgs_note_label = QLabel()
        self.pgs_note_label.setWordWrap(True)
        self.pgs_note_label.setStyleSheet("color: #1565C0; background-color: #e7f3fe; padding: 8px; border-radius: 4px;")
        pgs_layout.addWidget(self.pgs_note_label)
        
        pgs_link_layout = QHBoxLayout()
        pgs_link = QPushButton("🔗 Visit PGS Catalog")
        pgs_link.setStyleSheet("QPushButton { text-align: left; color: #337ab7; background: none; border: none; }")
        pgs_link.clicked.connect(lambda: QDesktopServices.openUrl(QUrl("https://www.pgscatalog.org/")))
        pgs_link_layout.addWidget(pgs_link)
        pgs_link_layout.addStretch()
        pgs_layout.addLayout(pgs_link_layout)
        
        layout.addWidget(pgs_group)
        
        # Backup section
        backup_group = QGroupBox("Backups")
        backup_layout = QVBoxLayout(backup_group)
        
        self.backup_list_label = QLabel("Available backups: 0")
        backup_layout.addWidget(self.backup_list_label)
        
        backup_btn_layout = QHBoxLayout()
        
        self.create_backup_btn = QPushButton("💾 Create Backup")
        self.create_backup_btn.clicked.connect(self._create_backup)
        backup_btn_layout.addWidget(self.create_backup_btn)
        
        backup_btn_layout.addStretch()
        backup_layout.addLayout(backup_btn_layout)
        
        layout.addWidget(backup_group)
        
        layout.addStretch()
        
        # Update instructions
        update_group = QGroupBox("How to Update")
        update_layout = QVBoxLayout(update_group)
        
        update_text = QLabel(
            "Use the buttons below to run database updates from inside the app. "
            "The output area will stream the updater logs while the process is running."
        )
        update_text.setWordWrap(True)
        update_layout.addWidget(update_text)

        update_button_layout = QHBoxLayout()
        self.gwas_update_btn = QPushButton("🔄 Update GWAS")
        self.gwas_update_btn.clicked.connect(lambda: self._start_database_update("gwas"))
        update_button_layout.addWidget(self.gwas_update_btn)

        self.pgs_update_btn = QPushButton("🔄 Update PGS")
        self.pgs_update_btn.clicked.connect(lambda: self._start_database_update("pgs"))
        update_button_layout.addWidget(self.pgs_update_btn)

        update_button_layout.addStretch()
        update_layout.addLayout(update_button_layout)

        self.update_output = QPlainTextEdit()
        self.update_output.setReadOnly(True)
        self.update_output.setMinimumHeight(180)
        self.update_output.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.update_output.setFont(QFont("Monospace"))
        self._reset_update_output()
        update_layout.addWidget(self.update_output)
        
        layout.addWidget(update_group)

    def _reset_update_output(self) -> None:
        """Populate the update panel with current instructions."""
        if not self.update_output:
            return

        self.update_output.setPlainText(
            "Update commands:\n"
            "  python database/update_databases.py --gwas\n"
            "  python database/update_databases.py --pgs\n"
            "  python database/update_databases.py --pgs --limit 50\n\n"
            "Live output from the last update appears below when you press a button."
        )

    def _set_update_buttons_enabled(self, enabled: bool) -> None:
        """Enable or disable the update buttons while a process is running."""
        self.gwas_update_btn.setEnabled(enabled)
        self.pgs_update_btn.setEnabled(enabled)

    def _append_update_output(self, text: str) -> None:
        """Append text to the update console."""
        if not self.update_output or not text:
            return

        self.update_output.appendPlainText(text.rstrip('\n'))

    def _start_database_update(self, command: str, limit: Optional[int] = None) -> None:
        """Run the database updater script in a background process."""
        if self.update_process and self.update_process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Update in progress", "A database update is already running.")
            return

        script_path = Path(BASE_DIR) / "database" / "update_databases.py"
        if not script_path.exists():
            QMessageBox.warning(self, "Updater not found", f"Could not locate {script_path}")
            return

        args = [str(script_path), f"--{command}"]
        if limit is not None:
            args.extend(["--limit", str(limit)])

        self.update_process = QProcess(self)
        self.update_process.setWorkingDirectory(str(BASE_DIR))
        self.update_process.readyReadStandardOutput.connect(self._on_update_stdout)
        self.update_process.readyReadStandardError.connect(self._on_update_stderr)
        self.update_process.finished.connect(self._on_update_finished)

        self._set_update_buttons_enabled(False)
        if self.update_output:
            self.update_output.appendPlainText("")
            self.update_output.appendPlainText("=" * 72)
            self.update_output.appendPlainText(f"Starting update: python database/update_databases.py --{command}")
            if limit is not None:
                self.update_output.appendPlainText(f"Limit: {limit}")
            self.update_output.appendPlainText("=" * 72)

        self.update_process.start(sys.executable, args)

    def _on_update_stdout(self) -> None:
        """Forward updater stdout to the live output console."""
        if not self.update_process or not self.update_output:
            return

        data = bytes(self.update_process.readAllStandardOutput()).decode('utf-8', errors='replace')
        if data:
            self._append_update_output(data)

    def _on_update_stderr(self) -> None:
        """Forward updater stderr to the live output console."""
        if not self.update_process or not self.update_output:
            return

        data = bytes(self.update_process.readAllStandardError()).decode('utf-8', errors='replace')
        if data:
            self._append_update_output(data)

    def _on_update_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        """Finalize the update run and refresh version labels."""
        if self.update_output:
            if exit_status == QProcess.ExitStatus.NormalExit and exit_code == 0:
                self.update_output.appendPlainText("")
                self.update_output.appendPlainText("Update completed successfully.")
            else:
                self.update_output.appendPlainText("")
                self.update_output.appendPlainText(f"Update finished with exit code {exit_code}.")

        self._set_update_buttons_enabled(True)
        self._refresh_versions()
        self.update_process = None
    
    def _refresh_versions(self) -> None:
        """Refresh version information display."""
        # GWAS version
        gwas_version = self.version_manager.get_gwas_version()
        if gwas_version:
            self.gwas_version_label.setText(f"Version: {gwas_version.version}")
            self.gwas_date_label.setText(
                f"Last Updated: {gwas_version.download_date.strftime('%Y-%m-%d %H:%M')}"
            )
            self.gwas_count_label.setText(f"Records: {gwas_version.record_count:,}")
        else:
            self.gwas_version_label.setText("Version: Not installed")
            self.gwas_date_label.setText("Last Updated: N/A")
            self.gwas_count_label.setText("Records: 0")
        
        # PGS version
        pgs_version = self.version_manager.get_pgs_version()
        if pgs_version:
            self.pgs_version_label.setText(f"Version: {pgs_version.version}")
            self.pgs_date_label.setText(
                f"Last Updated: {pgs_version.download_date.strftime('%Y-%m-%d %H:%M')}"
            )
            pgs_db = PolygenicDatabase()
            pgs_count = pgs_db.get_score_count()
            self.pgs_count_label.setText(f"Scores: {pgs_count}")

            if pgs_db.is_sample_pgs_database():
                self.pgs_note_label.setText(
                    "ℹ️ The database currently contains the bundled example scores that ship with the app. "
                    "Use the update button to download the full PGS Catalog."
                )
                self.pgs_note_label.setStyleSheet("color: #1565C0; background-color: #e7f3fe; padding: 8px; border-radius: 4px;")
            else:
                self.pgs_note_label.setText(
                    f"ℹ️ The database currently contains {pgs_count:,} scores from the PGS Catalog."
                )
                self.pgs_note_label.setStyleSheet("color: #155724; background-color: #d4edda; padding: 8px; border-radius: 4px;")
        else:
            self.pgs_version_label.setText("Version: Not installed")
            self.pgs_date_label.setText("Last Updated: N/A")
            self.pgs_count_label.setText("Scores: 0")
            self.pgs_note_label.setText(
                "ℹ️ No PGS database is installed yet. The app ships with example scores only when the bundled database is present."
            )
            self.pgs_note_label.setStyleSheet("color: #856404; background-color: #fff3cd; padding: 8px; border-radius: 4px;")
        
        # Backups
        backups = self.version_manager.list_backups()
        self.backup_list_label.setText(f"Available backups: {len(backups)}")
    
    def _create_backup(self) -> None:
        """Create database backups."""
        from config import DATABASE_PATH
        from database.polygenic_database import PGS_DATABASE_PATH
        
        gwas_backup = self.version_manager.create_backup(DATABASE_PATH)
        pgs_backup = self.version_manager.create_backup(PGS_DATABASE_PATH)
        
        if gwas_backup or pgs_backup:
            QMessageBox.information(
                self,
                "Backup Created",
                "Database backups have been created successfully."
            )
            self._refresh_versions()
        else:
            QMessageBox.warning(
                self,
                "Backup Failed",
                "Failed to create database backups."
            )
