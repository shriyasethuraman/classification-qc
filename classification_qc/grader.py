"""Quality checks for taxonomic classification output.

A metagenomic classification pipeline will happily return a confident-looking
answer from a run that went wrong: host depletion that failed, a reference
database missing the organism actually present, a reagent contaminant that
swamped a low-biomass sample. The classifier does not know any of this. It
reports what matched.

This module inspects classification output and flags the patterns that mean the
result should not be trusted, before anyone acts on it.
"""

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set

__all__ = [
    "Flag",
    "CRITICAL",
    "WARNING",
    "NICHE_EXPECTATIONS",
    "grade_classification",
]

# ---------------------------------------------------------------------------
# Severity levels
# ---------------------------------------------------------------------------

CRITICAL = "CRITICAL"
"""A specific upstream failure is implicated. The result should not be used."""

WARNING = "WARNING"
"""The output is suspicious and needs a human to look before it is trusted."""


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_UNCLASSIFIED_THRESHOLD = 0.50
"""Fraction of unclassified reads above which the run itself is suspect."""

DEFAULT_DOMINANCE_THRESHOLD = 0.80
"""Fraction of classified reads from one taxon above which diversity is suspect."""

NICHE_EXPECTATIONS: Dict[str, Set[str]] = {
    "gut": {"Bacteroides", "Firmicutes", "Escherichia", "Bifidobacterium"},
    "skin": {"Staphylococcus", "Cutibacterium", "Corynebacterium", "Streptococcus"},
    "soil": {"Streptomyces", "Bradyrhizobium", "Bacillus", "Acidobacteria"},
}
"""Genera plausibly expected per sample type.

Deliberately a small, legible starting point rather than an attempt at a
complete microbiome reference. Its purpose is to catch the obviously
out-of-place, not to adjudicate what belongs in a niche. Pass your own mapping
via ``niche_expectations`` to extend or replace it. A sample type absent from
this mapping disables the niche check rather than flagging everything in it.
"""

_UNCLASSIFIED = {"NA", "", None}
"""Values a pipeline may use to mean 'no classification assigned'."""


# ---------------------------------------------------------------------------
# Flag
# ---------------------------------------------------------------------------


@dataclass
class Flag:
    """A single quality concern raised about the classification output.

    ``read_id`` is ``None`` for run-level findings. The distinction matters:
    a per-read flag points at one observation, while a run-level flag says
    something about the pipeline or the sample as a whole, and the two warrant
    different responses.
    """

    read_id: Optional[str]
    flag_type: str
    severity: str
    message: str

    @property
    def is_run_level(self) -> bool:
        return self.read_id is None


# ---------------------------------------------------------------------------
# Grader
# ---------------------------------------------------------------------------


def grade_classification(
    sample_type: str,
    reads: Sequence[Mapping[str, object]],
    *,
    unclassified_threshold: float = DEFAULT_UNCLASSIFIED_THRESHOLD,
    dominance_threshold: float = DEFAULT_DOMINANCE_THRESHOLD,
    niche_expectations: Optional[Mapping[str, Iterable[str]]] = None,
) -> List[Flag]:
    """Inspect classified reads and report what makes the result untrustworthy.

    Args:
        sample_type: Where the sample came from, e.g. ``"gut"``, ``"skin"``,
            ``"soil"``. Used only by the niche-plausibility check; an unknown
            value disables that check rather than flagging every read.
        reads: One mapping per read, with keys ``read_id``, ``seq`` and
            ``tax_name``. ``tax_name`` is the assigned taxon, or ``"NA"`` /
            ``None`` where the classifier assigned nothing.
        unclassified_threshold: Unclassified fraction above which the run is
            flagged.
        dominance_threshold: Single-taxon fraction of classified reads above
            which diversity is flagged.
        niche_expectations: Override for :data:`NICHE_EXPECTATIONS`.

    Returns:
        Flags in the order found: per-read flags first, then run-level flags.
        An empty list means nothing suspicious was detected, which is not the
        same as the classification being correct.
    """
    flags: List[Flag] = []

    total_reads = len(reads)
    if total_reads == 0:
        # No reads is a caller-side condition, not a classification fault.
        return flags

    expectations = NICHE_EXPECTATIONS if niche_expectations is None else niche_expectations
    expected_genera = set(expectations.get(sample_type.lower(), set()))

    taxon_counts: Counter = Counter()
    unclassified_count = 0

    for read in reads:
        read_id = read.get("read_id")
        tax = read.get("tax_name")

        if tax in _UNCLASSIFIED:
            unclassified_count += 1
            continue

        tax = str(tax)
        taxon_counts[tax] += 1
        genus = tax.split()[0] if tax.split() else tax

        # Host DNA is checked first and short-circuits. It is the more specific
        # and more serious finding, and it applies whatever the sample type, so
        # reporting it alongside a generic "unexpected taxon" would only add
        # noise to the more important signal.
        if genus.lower() == "homo" or "sapiens" in tax.lower():
            flags.append(
                Flag(
                    read_id=read_id,
                    flag_type="HOST_DNA_LEAK",
                    severity=CRITICAL,
                    message=(
                        f"Read classified as human ({tax}); "
                        f"host depletion may have failed."
                    ),
                )
            )
            continue

        if expected_genera and genus not in expected_genera:
            flags.append(
                Flag(
                    read_id=read_id,
                    flag_type="UNEXPECTED_NICHE_TAXON",
                    severity=WARNING,
                    message=(
                        f"Taxon '{tax}' is not among the genera expected for "
                        f"sample type '{sample_type}'."
                    ),
                )
            )

    # --- Run-level checks ---------------------------------------------------

    unclassified_ratio = unclassified_count / total_reads
    if unclassified_ratio > unclassified_threshold:
        flags.append(
            Flag(
                read_id=None,
                flag_type="HIGH_UNCLASSIFIED_RATE",
                severity=CRITICAL,
                message=(
                    f"{unclassified_ratio:.1%} of reads are unclassified; "
                    f"the reference database may not cover what is in this "
                    f"sample, or the run quality may be poor."
                ),
            )
        )

    classified_total = total_reads - unclassified_count
    if classified_total > 0:
        top_taxon, top_count = taxon_counts.most_common(1)[0]
        top_fraction = top_count / classified_total
        if top_fraction > dominance_threshold:
            flags.append(
                Flag(
                    read_id=None,
                    flag_type="DOMINANT_TAXON",
                    severity=WARNING,
                    message=(
                        f"'{top_taxon}' accounts for {top_fraction:.1%} of "
                        f"classified reads; consider a contamination spike or "
                        f"processing artifact rather than true community "
                        f"structure."
                    ),
                )
            )

    return flags
