"""Quality checks for taxonomic classification output."""

from .grader import (
    CRITICAL,
    NICHE_EXPECTATIONS,
    WARNING,
    Flag,
    grade_classification,
)

__version__ = "1.0.0"

__all__ = [
    "Flag",
    "CRITICAL",
    "WARNING",
    "NICHE_EXPECTATIONS",
    "grade_classification",
]
