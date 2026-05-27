from __future__ import annotations

from collections.abc import Callable
import math
from pathlib import Path
import csv

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord

from .models import FixedTopologySiteLogLikelihoodReport
from .models import SiteLogLikelihoodRow
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


def expanded_site_log_likelihood_rows_from_patterns(
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> tuple[list[SiteLogLikelihoodRow], float]:
    """Expand compressed-pattern log likelihoods back to one stable row per site."""
    rows: list[SiteLogLikelihoodRow] = []
    total = 0.0
    for pattern in compressed_patterns.patterns:
        log_likelihood = _validated_site_log_likelihood(
            site_log_likelihood(pattern.states)
        )
        total += pattern.weight * log_likelihood
        for site_position in pattern.site_positions:
            rows.append(
                SiteLogLikelihoodRow(
                    pattern_id=pattern.pattern_id,
                    pattern_weight=pattern.weight,
                    site_position=site_position,
                    site_states=pattern.states,
                    log_likelihood=log_likelihood,
                )
            )
    return rows, total


def write_site_log_likelihood_table(
    path: Path,
    report: FixedTopologySiteLogLikelihoodReport,
) -> Path:
    """Write one expanded per-site log-likelihood TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model_name",
                "taxon_order",
                "pattern_id",
                "pattern_weight",
                "site_position",
                "site_states",
                "log_likelihood",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        taxon_order = "|".join(report.taxa)
        for row in report.site_log_likelihoods:
            writer.writerow(
                {
                    "model_name": report.model_name,
                    "taxon_order": taxon_order,
                    "pattern_id": row.pattern_id,
                    "pattern_weight": row.pattern_weight,
                    "site_position": row.site_position,
                    "site_states": "|".join(row.site_states),
                    "log_likelihood": repr(row.log_likelihood),
                }
            )
    return path


def _validated_site_log_likelihood(log_likelihood: float) -> float:
    if not math.isfinite(log_likelihood):
        raise ValueError("site log likelihood must be finite")
    return log_likelihood
