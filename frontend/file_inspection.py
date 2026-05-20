"""UI helpers for pre-analysis genetic file inspection."""

from typing import List, Optional, Tuple

from PyQt6.QtCore import QObject, QThread, pyqtSignal, Qt
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from backend.file_inspector import FileInspectionResult, inspect_genetic_file
from utils.logging_config import get_logger

logger = get_logger(__name__)


class FileInspectionWorker(QThread):
    """Background worker that inspects a genetic file before analysis."""

    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    progress = pyqtSignal(int, str)

    def __init__(self, filepath: str) -> None:
        super().__init__()
        self.filepath = filepath

    def run(self) -> None:
        try:
            result = inspect_genetic_file(self.filepath, progress_callback=self._on_progress)
            self.finished.emit(result)
        except Exception as exc:
            logger.exception("File inspection failed: %s", exc)
            self.error.emit(str(exc))

    def _on_progress(self, percent: int, message: str) -> None:
        self.progress.emit(percent, message)


class FileStatsDialog(QDialog):
    """Modal dialog that displays file inspection statistics."""

    def __init__(self, inspection: FileInspectionResult, parent=None) -> None:
        super().__init__(parent)
        self.inspection = inspection
        self.setWindowTitle(f"File Statistics - {inspection.filename}")
        self.setMinimumSize(860, 720)
        self._init_ui()

    def _init_ui(self) -> None:
        layout = QVBoxLayout(self)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(8)

        self._build_widget_content(content_layout)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        button_box = QDialogButtonBox()
        if self.inspection.analysis_supported:
            button_box.addButton("Analyze", QDialogButtonBox.ButtonRole.AcceptRole)
            button_box.addButton("Cancel", QDialogButtonBox.ButtonRole.RejectRole)
        else:
            button_box.addButton("Close", QDialogButtonBox.ButtonRole.AcceptRole)

        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _build_widget_content(self, layout: QVBoxLayout) -> None:
        stats = self.inspection.stats or {}
        genetic_records = stats.get("genetic_records")

        self._add_intro(layout)
        self._add_section(layout, "Common overview", "Common fields shared by 23andMe and VCF files.")
        common_rows: List[Tuple[str, str, str, str]] = [
            (
                "📄",
                "File",
                self.inspection.filename,
                "Name of the selected input file.",
            ),
            (
                "🏷️",
                "Detected format",
                self.inspection.detected_format,
                "Format inferred from the file contents and extension.",
            ),
            (
                "📏",
                "File size",
                self.inspection.file_size_human,
                "Size on disk before any decompression or parsing.",
            ),
            (
                "🗜️",
                "Compressed",
                "Yes" if self.inspection.is_compressed else "No",
                "Indicates whether the input is gzipped and read as a compressed stream.",
            ),
            (
                "🔢",
                "Total lines scanned",
                self._fmt(stats.get("total_lines")),
                "All lines read from the file, including comments and headers.",
            ),
            (
                "🧬",
                "Genetic records read",
                self._fmt(genetic_records),
                "Validated genetic records found in the file: SNPs for 23andMe or variant records for VCF.",
            ),
            (
                "⚠️",
                "Warnings reported",
                self._fmt(stats.get("warnings_count")),
                "Non-fatal issues encountered while reading the file.",
            ),
        ]

        self._add_rows(layout, common_rows)

        if self.inspection.analysis_supported:
            self._add_section(layout, "23andMe-specific metrics", "Quality and genotype breakdown for 23andMe raw files.")
            self._add_rows(layout, [
                ("💬", "Comment lines", self._fmt(stats.get("comment_lines")), "Lines starting with # that were skipped as comments."),
                ("▫️", "Blank lines", self._fmt(stats.get("blank_lines")), "Empty lines that contained no data."),
                ("✖️", "Invalid lines", self._fmt(stats.get("invalid_lines")), "Lines that looked like data but failed validation."),
                ("—", "Undetermined genotypes", self._fmt(stats.get("skipped_undetermined")), "Rows where the genotype was reported as -- and skipped."),
                ("🧬", "Homozygous SNPs", self._fmt(stats.get("homozygous_snps")), "Validated SNPs where both alleles were the same."),
                ("🧬", "Heterozygous SNPs", self._fmt(stats.get("heterozygous_snps")), "Validated SNPs where the two alleles differed."),
            ])

            self._add_section(layout, "Shared distributions", "Distribution tables common to both formats.")
            self._add_distribution_block(
                layout,
                "Chromosome distribution",
                stats.get("chromosome_counts", {}),
                item_tooltip="Chromosome name or number for each observed record.",
                value_tooltip="Number of validated SNPs on this chromosome.",
            )
            self._add_distribution_block(
                layout,
                "Genotype distribution",
                stats.get("genotype_counts", {}),
                item_tooltip="Observed two-allele genotype string.",
                value_tooltip="Number of times this genotype appeared in the file.",
            )
            self._add_distribution_block(
                layout,
                "Allele counts",
                stats.get("allele_counts", {}),
                item_tooltip="Individual allele symbol observed across validated genotypes.",
                value_tooltip="Total allele observations across all valid genotypes.",
            )
            warnings = stats.get("warnings", [])
        else:
            self._notice_widget(
                layout,
                "VCF analysis is supported for biallelic SNPs. Structural, multiallelic, and non-SNP records are kept in the inspection summary and skipped from downstream matching.",
            )
            self._add_section(layout, "VCF-specific metrics", "Structure and genotype summary for VCF files.")
            self._add_rows(layout, [
                ("ℹ️", "Metadata lines", self._fmt(stats.get("metadata_lines")), "Lines beginning with ## that store VCF metadata."),
                ("🧾", "Header lines", self._fmt(stats.get("header_lines")), "The #CHROM header line that defines the VCF columns and samples."),
                ("↔️", "Column count", self._fmt(stats.get("column_count")), "Number of columns found in the VCF header line."),
                ("👥", "Samples detected", self._fmt(stats.get("sample_count")), "Number of sample columns after the FORMAT field."),
                ("🏷️", "VCF version", self._escape(stats.get("vcf_version")), "Version string reported by the fileformat header."),
                ("✅", "Called genotypes", self._fmt(stats.get("called_genotypes")), "Sample genotype calls with an observed GT value."),
                ("∅", "Missing genotypes", self._fmt(stats.get("missing_genotypes")), "Sample genotype calls with missing or undefined GT values."),
                ("🔒", "Phased genotypes", self._fmt(stats.get("phased_genotypes")), "GT calls that preserve haplotype phase using |."),
                ("🔓", "Unphased genotypes", self._fmt(stats.get("unphased_genotypes")), "GT calls that do not preserve haplotype phase and use /."),
                ("🔁", "SNP variants", self._fmt(stats.get("snv_count")), "VCF records where REF and ALT are single-base substitutions."),
                ("➕", "Indels", self._fmt(stats.get("indel_count")), "VCF records representing insertions or deletions."),
                ("☰", "Multiallelic variants", self._fmt(stats.get("multiallelic_count")), "VCF records with more than one alternate allele."),
            ])

            self._add_section(layout, "Shared distributions", "Distribution tables common to both formats.")
            self._add_distribution_block(
                layout,
                "Variant classes",
                stats.get("variant_types", {}),
                item_tooltip="Class of the VCF record based on REF/ALT contents.",
                value_tooltip="Number of records in this class.",
            )
            self._add_distribution_block(
                layout,
                "Chromosome distribution",
                stats.get("chromosome_counts", {}),
                item_tooltip="Chromosome name or number for each VCF record.",
                value_tooltip="Number of VCF records on this chromosome.",
            )
            self._add_distribution_block(
                layout,
                "Genotype distribution",
                stats.get("genotype_counts", {}),
                item_tooltip="Observed GT call for a sample at each record.",
                value_tooltip="Number of sample genotype calls with this GT value.",
            )
            self._add_distribution_block(
                layout,
                "Filter distribution",
                stats.get("filter_counts", {}),
                item_tooltip="Value in the FILTER column for each record.",
                value_tooltip="Number of records with this FILTER value.",
            )

            self._add_section(layout, "VCF samples", "Sample column names declared in the VCF header.")
            self._add_sample_rows(layout, stats.get("samples", []), tooltip="Sample column names declared in the VCF header.")
            warnings = stats.get("warnings", [])

        if warnings:
            self._add_section(layout, "Warnings", "Non-fatal issues encountered while reading the file.")
            self._add_sample_rows(layout, [str(item) for item in warnings], tooltip="Non-fatal issues encountered while reading the file.", icon="⚠️")

        layout.addStretch()

    def _add_intro(self, layout: QVBoxLayout) -> None:
        title = QLabel("<h2>File inspection</h2>")
        title.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(title)

        description = QLabel("This popup shows the file structure and genetic content before any database matching runs.")
        description.setWordWrap(True)
        description.setStyleSheet("color: #555;")
        layout.addWidget(description)

    def _add_section(self, layout: QVBoxLayout, title: str, tooltip: str) -> None:
        section = QLabel(title)
        section.setToolTip(tooltip)
        section.setStyleSheet("font-size: 13pt; font-weight: 600; margin-top: 8px; color: #222;")
        layout.addWidget(section)

    def _add_rows(self, layout: QVBoxLayout, rows: List[Tuple[str, str, str, str]]) -> None:
        for icon, label, value, tooltip in rows:
            layout.addWidget(self._create_row(icon, label, value, tooltip))

    def _add_distribution_block(
        self,
        layout: QVBoxLayout,
        title: str,
        distribution: dict,
        item_tooltip: str,
        value_tooltip: str,
    ) -> None:
        if not distribution:
            empty = QLabel("Not available.")
            empty.setStyleSheet("color: #777; margin-left: 8px;")
            empty.setToolTip(title)
            layout.addWidget(empty)
            return

        sub_title = QLabel(title)
        sub_title.setToolTip(title)
        sub_title.setStyleSheet("font-weight: 600; margin-top: 4px; color: #333;")
        layout.addWidget(sub_title)

        items = sorted(distribution.items(), key=lambda item: (-item[1], str(item[0])))
        for key, value in items:
            layout.addWidget(
                self._create_row(
                    "🔹",
                    str(key),
                    self._fmt(value),
                    f"{item_tooltip} {value_tooltip}".strip(),
                    value_tooltip=value_tooltip,
                )
            )

    def _add_sample_rows(self, layout: QVBoxLayout, samples: List[str], tooltip: str = "", icon: str = "👤") -> None:
        if not samples:
            empty = QLabel("No sample columns were detected.")
            empty.setStyleSheet("color: #777; margin-left: 8px;")
            empty.setToolTip(tooltip)
            layout.addWidget(empty)
            return

        for sample in samples[:10]:
            layout.addWidget(self._create_row(icon, sample, "", tooltip))

        if len(samples) > 10:
            more = QLabel(f"and {len(samples) - 10} more sample(s).")
            more.setStyleSheet("color: #777; font-style: italic; margin-left: 8px;")
            more.setToolTip(tooltip)
            layout.addWidget(more)

    def _notice_widget(self, layout: QVBoxLayout, text: str) -> None:
        notice = QFrame()
        notice.setStyleSheet("background:#fff4d6; border:1px solid #e0c27a; border-radius:8px;")
        notice_layout = QVBoxLayout(notice)
        notice_layout.setContentsMargins(12, 10, 12, 10)
        label = QLabel(text)
        label.setWordWrap(True)
        label.setToolTip(text)
        notice_layout.addWidget(label)
        layout.addWidget(notice)

    def _create_row(self, icon: str, label: str, value: str, tooltip: str, value_tooltip: str = "") -> QFrame:
        row = QFrame()
        row.setToolTip(tooltip)
        row.setStyleSheet(
            "QFrame { background: white; border: 1px solid #e5e5e5; border-radius: 6px; }"
        )
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(10, 6, 10, 6)
        row_layout.setSpacing(10)

        icon_label = QLabel(icon)
        icon_label.setToolTip(tooltip)
        icon_label.setFixedWidth(24)
        row_layout.addWidget(icon_label)

        field_label = QLabel(label)
        field_label.setToolTip(tooltip)
        field_label.setStyleSheet("font-weight: 600; color: #222;")
        field_label.setWordWrap(True)
        row_layout.addWidget(field_label, stretch=1)

        value_label = QLabel(value)
        value_label.setToolTip(value_tooltip or tooltip)
        value_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        value_label.setStyleSheet("color: #333;")
        value_label.setWordWrap(True)
        row_layout.addWidget(value_label)

        return row

    def _fmt(self, value) -> str:
        if value is None:
            return "-"
        if isinstance(value, float):
            return f"{value:.2f}"
        try:
            return f"{int(value):,}"
        except (TypeError, ValueError):
            return self._escape(str(value))

    def _escape(self, text) -> str:
        if text is None:
            return "-"
        return (
            str(text)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
        )


class FileInspectionCoordinator(QObject):
    """Coordinates file inspection and the pre-analysis dialog flow."""

    def __init__(self, main_window) -> None:
        super().__init__(main_window)
        self.main_window = main_window
        self.worker: Optional[FileInspectionWorker] = None
        self._upload_was_enabled = False

    def start(self, filepath: str) -> None:
        """Inspect a file before any analysis is launched."""
        if self.worker and self.worker.isRunning():
            return

        self._upload_was_enabled = self.main_window.upload_btn.isEnabled()
        self.main_window.upload_btn.setEnabled(False)
        self.main_window.status_bar.showMessage("Inspecting genetic file...")

        self.worker = FileInspectionWorker(filepath)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.progress.connect(self._on_progress)
        self.worker.start()

    def _on_progress(self, percent: int, message: str) -> None:
        logger.info("File inspection progress: %s%% - %s", percent, message)
        self.main_window.status_bar.showMessage(f"Inspecting genetic file: {percent}% - {message}")

    def _on_finished(self, inspection: FileInspectionResult) -> None:
        self.main_window._current_inspection_result = inspection
        dialog = FileStatsDialog(inspection, self.main_window)
        dialog_result = dialog.exec()

        if dialog_result == QDialog.DialogCode.Accepted:
            self.main_window._start_processing(
                inspection.filepath,
                snp_records=inspection.records,
                inspection_stats=inspection.stats,
            )
            return

        self._restore_upload_button()
        self.main_window.status_bar.showMessage("File inspection cancelled")

    def _on_error(self, error_message: str) -> None:
        self._restore_upload_button()
        QMessageBox.critical(self.main_window, "Inspection Error", error_message)
        self.main_window.status_bar.showMessage(f"Error: {error_message}")

    def _restore_upload_button(self) -> None:
        self.main_window.upload_btn.setEnabled(self._upload_was_enabled)
