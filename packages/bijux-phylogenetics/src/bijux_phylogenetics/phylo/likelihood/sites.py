from __future__ import annotations

from collections.abc import Callable
import csv
import math
from pathlib import Path
from typing import Protocol

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


class PatternedSiteLogLikelihoodRow(Protocol):
    """Duck-typed site-likelihood row contract shared across sequence reports."""

    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    log_likelihood: float


def sum_alignment_site_log_likelihoods(
    records: list[AlignmentRecord],
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> float:
    """Sum one per-site log likelihood across every uncompressed alignment column."""
    return math.fsum(
        _validated_site_log_likelihood(site_log_likelihood(states))
        for states in iter_uncompressed_alignment_sites(records)
    )


def sum_compressed_site_pattern_log_likelihoods(
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> float:
    """Sum one per-pattern log likelihood using each pattern's integer weight."""
    return math.fsum(
        pattern.weight
        * _validated_site_log_likelihood(site_log_likelihood(pattern.states))
        for pattern in compressed_patterns.patterns
    )


def expanded_site_log_likelihood_rows_from_patterns(
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    site_log_likelihood: Callable[[tuple[str, ...]], float],
) -> tuple[list[SiteLogLikelihoodRow], float]:
    """Expand compressed-pattern log likelihoods back to one stable row per site."""
    log_likelihood_by_pattern_id: dict[str, float] = {}
    weighted_pattern_terms: list[float] = []
    for pattern in compressed_patterns.patterns:
        log_likelihood = _validated_site_log_likelihood(
            site_log_likelihood(pattern.states)
        )
        log_likelihood_by_pattern_id[pattern.pattern_id] = log_likelihood
        weighted_pattern_terms.append(pattern.weight * log_likelihood)

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
    return rows, math.fsum(weighted_pattern_terms)


def sum_site_log_likelihood_rows(
    rows: list[PatternedSiteLogLikelihoodRow],
) -> float:
    """Sum expanded site rows exactly once per emitted site."""
    return math.fsum(_validated_site_log_likelihood(row.log_likelihood) for row in rows)


def sum_weighted_site_pattern_log_likelihood_rows(
    rows: list[PatternedSiteLogLikelihoodRow],
) -> float:
    """Sum one weighted term per unique pattern after validating repeated rows."""
    pattern_terms: dict[str, tuple[int, tuple[str, ...], float]] = {}
    for row in rows:
        log_likelihood = _validated_site_log_likelihood(row.log_likelihood)
        prior = pattern_terms.get(row.pattern_id)
        current = (row.pattern_weight, row.site_states, log_likelihood)
        if prior is None:
            pattern_terms[row.pattern_id] = current
            continue
        if prior != current:
            raise ValueError(
                "site log-likelihood rows require stable pattern weights, states, and log likelihoods"
            )
    return math.fsum(
        pattern_weight * log_likelihood
        for pattern_weight, _states, log_likelihood in pattern_terms.values()
    )


def validate_site_log_likelihood_reconstruction(
    rows: list[PatternedSiteLogLikelihoodRow],
    *,
    expected_total_log_likelihood: float,
    expected_site_count: int,
    expected_pattern_count: int,
    owner_name: str,
) -> None:
    """Reject emitted site rows that do not reconstruct their declared total."""
    if len(rows) != expected_site_count:
        raise ValueError(
            f"{owner_name} emitted {len(rows)} site rows for {expected_site_count} declared sites"
        )
    unique_pattern_count = len({row.pattern_id for row in rows})
    if unique_pattern_count != expected_pattern_count:
        raise ValueError(
            f"{owner_name} emitted {unique_pattern_count} patterns for {expected_pattern_count} declared patterns"
        )
    expanded_total = sum_site_log_likelihood_rows(rows)
    weighted_pattern_total = sum_weighted_site_pattern_log_likelihood_rows(rows)
    if not math.isclose(
        expanded_total,
        expected_total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            f"{owner_name} expanded site rows did not reconstruct the declared total likelihood"
        )
    if not math.isclose(
        weighted_pattern_total,
        expected_total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            f"{owner_name} weighted pattern rows did not reconstruct the declared total likelihood"
        )


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
