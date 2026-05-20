"""Common parsing primitives for genetic file formats."""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional


class ParseError(Exception):
    """Exception raised when parsing or inspection fails."""

    pass


@dataclass
class ParserIssue:
    """Structured parser warning or error message."""

    severity: str
    code: str
    message: str
    line_number: Optional[int] = None
    details: Dict[str, Any] = field(default_factory=dict)
