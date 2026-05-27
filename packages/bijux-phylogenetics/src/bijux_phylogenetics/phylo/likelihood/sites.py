from __future__ import annotations

from collections.abc import Callable
import math

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord

from .patterns import CompressedAlignmentSitePatterns
from .patterns import iter_uncompressed_alignment_sites


def sum_alignment_site_log_likelihoods(
    records: list[AlignmentRecord],
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> float:
    """Sum one per-site log likelihood across every uncompressed alignment column."""
    total = 0.0
    for states in iter_uncompressed_alignment_sites(records):
        total += _validated_site_log_likelihood(site_log_likelihood(states))
    return total


def sum_compressed_site_pattern_log_likelihoods(
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> float:
    """Sum one per-pattern log likelihood using each pattern's integer weight."""
    total = 0.0
    for pattern in compressed_patterns.patterns:
        total += pattern.weight * _validated_site_log_likelihood(
            site_log_likelihood(pattern.states)
        )
    return total


def _validated_site_log_likelihood(log_likelihood: float) -> float:
    if not math.isfinite(log_likelihood):
        raise ValueError("site log likelihood must be finite")
    return log_likelihood
