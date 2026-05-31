from __future__ import annotations

from pathlib import Path

from .genetic_distance_matrix import (
    compute_pairwise_genetic_distance_matrix,
)
from .models import (
    AmbiguityPolicy,
    DistanceModel,
    DistanceSaturationDiagnosticsReport,
    DistanceSaturationWarning,
    GapHandlingMode,
    GeneticDistanceMatrix,
    SaturatedDistancePair,
)


def _warning_kind(reason: str, *, distance: float | None) -> str:
    lowered = reason.lower()
    if "tends to infinity" in lowered or "correction limit" in lowered:
        return "infinite-corrected-distance"
    if distance is None:
        return "undefined-corrected-distance"
    if "correction range" in lowered or "undefined" in lowered:
        return "undefined-corrected-distance"
    return "saturation-risk"


def diagnose_distance_saturation_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
) -> DistanceSaturationDiagnosticsReport:
    """Summarize pair-level saturation warnings before distance-tree inference."""
    warning_rows: list[DistanceSaturationWarning] = []
    for pair in report.pairs:
        if pair.left_identifier == pair.right_identifier or not pair.saturated:
            continue
        reason = pair.saturation_reason or "distance enters a saturation regime"
        kind = _warning_kind(reason, distance=pair.distance)
        blocks_tree_inference = pair.distance is None
        warning_rows.append(
            DistanceSaturationWarning(
                left_identifier=pair.left_identifier,
                right_identifier=pair.right_identifier,
                distance=pair.distance,
                comparable_sites=pair.comparable_sites,
                warning_kind=kind,
                reason=reason,
                blocks_tree_inference=blocks_tree_inference,
            )
        )
    warning_rows.sort(
        key=lambda row: (
            not row.blocks_tree_inference,
            row.left_identifier,
            row.right_identifier,
            row.warning_kind,
        )
    )
    warnings = list(report.warnings)
    if warning_rows:
        warnings.append(
            "one or more pairwise distances enter a saturation regime before tree inference"
        )
    if any(row.blocks_tree_inference for row in warning_rows):
        warnings.append(
            "one or more corrected distances are undefined or infinite, so distance-tree inference is blocked"
        )
    return DistanceSaturationDiagnosticsReport(
        alignment_path=report.path,
        model=report.model,
        gap_handling=report.gap_handling,
        ambiguity_policy=report.ambiguity_policy,
        taxon_count=len(report.identifiers),
        pair_count=len(report.pairs),
        blocking_warning_count=sum(
            1 for row in warning_rows if row.blocks_tree_inference
        ),
        warning_rows=warning_rows,
        warnings=warnings,
        blocks_tree_inference=any(row.blocks_tree_inference for row in warning_rows),
    )


def diagnose_distance_saturation(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceSaturationDiagnosticsReport:
    """Summarize pair-level saturation warnings from one alignment."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    return diagnose_distance_saturation_from_genetic_distance_matrix(matrix)


def saturated_pairs_from_diagnostics(
    report: DistanceSaturationDiagnosticsReport,
) -> list[SaturatedDistancePair]:
    """Project saturation diagnostics onto the existing saturated-pair summary rows."""
    return [
        SaturatedDistancePair(
            left_identifier=row.left_identifier,
            right_identifier=row.right_identifier,
            distance=row.distance,
            comparable_sites=row.comparable_sites,
            reason=row.reason,
        )
        for row in report.warning_rows
    ]
