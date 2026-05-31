from __future__ import annotations

from collections.abc import Callable
import csv
import math
from pathlib import Path

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord

from .models import (
    FixedTopologySiteLogLikelihoodReport,
    LocalClockLikelihoodReport,
    ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport,
    ProteinEmpiricalDiscreteGammaTreeLikelihoodReport,
    ProteinEmpiricalInvariantMixtureTreeLikelihoodReport,
    ProteinEmpiricalMatrixTreeLikelihoodReport,
    ProteinPoissonTreeLikelihoodReport,
    SiteLogLikelihoodRow,
    StrictClockLikelihoodReport,
)
from .patterns import (
    CompressedAlignmentSitePatterns,
    iter_pattern_sites_in_alignment_order,
    iter_uncompressed_alignment_sites,
)

SequenceSiteLogLikelihoodReport = (
    FixedTopologySiteLogLikelihoodReport
    | LocalClockLikelihoodReport
    | ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport
    | ProteinEmpiricalDiscreteGammaTreeLikelihoodReport
    | ProteinEmpiricalInvariantMixtureTreeLikelihoodReport
    | ProteinEmpiricalMatrixTreeLikelihoodReport
    | ProteinPoissonTreeLikelihoodReport
    | StrictClockLikelihoodReport
)


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
    log_likelihood_by_pattern_id: dict[str, float] = {}
    total = 0.0
    for pattern in compressed_patterns.patterns:
        log_likelihood = _validated_site_log_likelihood(
            site_log_likelihood(pattern.states)
        )
        log_likelihood_by_pattern_id[pattern.pattern_id] = log_likelihood
        total += pattern.weight * log_likelihood

    rows: list[SiteLogLikelihoodRow] = []
    for site_position, pattern in iter_pattern_sites_in_alignment_order(
        compressed_patterns
    ):
        rows.append(
            SiteLogLikelihoodRow(
                pattern_id=pattern.pattern_id,
                pattern_weight=pattern.weight,
                site_position=site_position,
                site_states=pattern.states,
                log_likelihood=log_likelihood_by_pattern_id[pattern.pattern_id],
            )
        )
    return rows, total


def write_site_log_likelihood_table(
    path: Path,
    report: SequenceSiteLogLikelihoodReport,
) -> Path:
    """Write one expanded per-site log-likelihood TSV for one sequence model."""
    path.parent.mkdir(parents=True, exist_ok=True)
    model_name = _site_log_likelihood_report_model_name(report)
    taxon_order = "|".join(report.taxa)
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
        for row in _iter_sequence_site_log_likelihood_rows(report):
            writer.writerow(
                {
                    "model_name": model_name,
                    "taxon_order": taxon_order,
                    "pattern_id": row.pattern_id,
                    "pattern_weight": row.pattern_weight,
                    "site_position": row.site_position,
                    "site_states": "|".join(row.site_states),
                    "log_likelihood": repr(row.log_likelihood),
                }
            )
    return path


def _iter_sequence_site_log_likelihood_rows(
    report: SequenceSiteLogLikelihoodReport,
) -> list[SiteLogLikelihoodRow]:
    if hasattr(report, "site_log_likelihoods"):
        return list(report.site_log_likelihoods)
    if isinstance(report, ProteinEmpiricalDiscreteGammaTreeLikelihoodReport):
        return [
            SiteLogLikelihoodRow(
                pattern_id=row.pattern_id,
                pattern_weight=row.pattern_weight,
                site_position=row.site_position,
                site_states=row.site_states,
                log_likelihood=row.log_likelihood,
            )
            for row in report.site_likelihoods
        ]
    if isinstance(report, ProteinEmpiricalInvariantMixtureTreeLikelihoodReport):
        return [
            SiteLogLikelihoodRow(
                pattern_id=row.pattern_id,
                pattern_weight=row.pattern_weight,
                site_position=row.site_position,
                site_states=row.site_states,
                log_likelihood=row.log_likelihood,
            )
            for row in report.site_likelihoods
        ]
    if isinstance(report, ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport):
        return [
            SiteLogLikelihoodRow(
                pattern_id=row.pattern_id,
                pattern_weight=row.pattern_weight,
                site_position=row.site_position,
                site_states=row.site_states,
                log_likelihood=row.log_likelihood,
            )
            for row in report.site_likelihoods
        ]
    raise TypeError(
        "sequence site log-likelihood export requires a supported sequence likelihood report"
    )


def _site_log_likelihood_report_model_name(
    report: SequenceSiteLogLikelihoodReport,
) -> str:
    if hasattr(report, "model_name"):
        return report.model_name
    if isinstance(report, ProteinPoissonTreeLikelihoodReport):
        return "protein Poisson"
    if isinstance(report, ProteinEmpiricalMatrixTreeLikelihoodReport):
        return "empirical protein matrix"
    if isinstance(report, ProteinEmpiricalDiscreteGammaTreeLikelihoodReport):
        return "empirical protein matrix +G"
    if isinstance(report, ProteinEmpiricalInvariantMixtureTreeLikelihoodReport):
        return "empirical protein matrix +I"
    if isinstance(report, ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport):
        return "empirical protein matrix +G+I"
    raise TypeError(
        "sequence site log-likelihood export requires a supported sequence likelihood report"
    )


def _validated_site_log_likelihood(log_likelihood: float) -> float:
    if not math.isfinite(log_likelihood):
        raise ValueError("site log likelihood must be finite")
    return log_likelihood
