"""Compatibility wrapper for legacy parser imports."""

from backend.parser_common import ParseError
from backend.parsers_23andme import Parser23andMe, parse_23andme_file

__all__ = ["ParseError", "Parser23andMe", "parse_23andme_file"]
