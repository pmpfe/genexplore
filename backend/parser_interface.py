"""Common parser interface for supported genetic file formats."""

from abc import ABC, abstractmethod
from typing import Callable, List, Optional

from models.data_models import SNPRecord


class GeneticFileParser(ABC):
    """Base interface implemented by genetic file parsers."""

    format_name: str = 'unknown'

    @classmethod
    @abstractmethod
    def can_parse(cls, filepath: str) -> bool:
        """Return True when the parser can likely handle the file."""

    @abstractmethod
    def inspect_file(
        self,
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ):
        """Inspect a file and return statistics plus optional normalized records."""

    @abstractmethod
    def parse_file(
        self,
        filepath: str,
        progress_callback: Optional[Callable[[int, str], None]] = None,
    ) -> List[SNPRecord]:
        """Return normalized variant records from the file."""