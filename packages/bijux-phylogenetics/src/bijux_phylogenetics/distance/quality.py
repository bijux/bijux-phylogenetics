from __future__ import annotations

import math
from pathlib import Path
from statistics import median

from .genetic_distance_matrix import (
    _build_alignment_distance_lookup,
    compute_pairwise_genetic_distance_matrix,
)
from .models import (
    AmbiguityPolicy,
    DistanceMatrixQualityReport,
    DistanceMethodAssessment,
    DistanceMethodAssumptionReport,
    DistanceModel,
    DistanceOutlierPair,
    GapHandlingMode,
    GeneticDistanceMatrix,
    LowInformationPair,
)
from .shared import _iter_ultrametric_violations
from .saturation import (
    diagnose_distance_saturation_from_genetic_distance_matrix,
    saturated_pairs_from_diagnostics,
)


def assess_distance_method_assumptions(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit clock-like compatibility and core distance-tree assumptions for an alignment."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    distances = _build_alignment_distance_lookup(matrix)
    expected_pairs = (len(matrix.identifiers) * (len(matrix.identifiers) - 1)) // 2
    violations = _iter_ultrametric_violations(
        matrix.identifiers,
        distances,
        tolerance=ultrametric_tolerance,
    )
    warnings: list[str] = []
    if len(distances) < expected_pairs:
        warnings.append(
            "UPGMA clock-assumption auditing is incomplete because one or more pairwise distances are undefined under the selected model"
        )
    if violations:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=path,
        source_kind="alignment",
        taxon_count=len(matrix.identifiers),
        pair_count=len(distances),
        nj_assumptions=[
            "neighbor-joining treats the matrix as an additive approximation and does not require a strict molecular clock",
            "neighbor-joining still becomes unreliable when pairwise distances are heavily saturated or estimated from too few comparable sites",
        ],
        upgma_assumptions=[
            "UPGMA assumes the pairwise distances are ultrametric and therefore consistent with a clock-like process",
            "UPGMA can mis-cluster taxa when rates vary among lineages even if the matrix remains symmetric and complete",
        ],
        ultrametric_compatible=not violations,
        ultrametric_tolerance=ultrametric_tolerance,
        upgma_ultrametric_violations=violations,
        warnings=warnings,
    )


def assess_distance_method_assumptions_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit core distance-tree assumptions from one in-memory distance matrix."""
    distances = _build_alignment_distance_lookup(report)
    expected_pairs = (len(report.identifiers) * (len(report.identifiers) - 1)) // 2
    violations = _iter_ultrametric_violations(
        report.identifiers,
        distances,
        tolerance=ultrametric_tolerance,
    )
    warnings = list(report.warnings)
    if len(distances) < expected_pairs:
        warnings.append(
            "UPGMA clock-assumption auditing is incomplete because one or more pairwise distances are undefined under the selected model"
        )
    if violations:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=report.path,
        source_kind="alignment",
        taxon_count=len(report.identifiers),
        pair_count=len(distances),
        nj_assumptions=[
            "neighbor-joining treats the matrix as an additive approximation and does not require a strict molecular clock",
            "neighbor-joining still becomes unreliable when pairwise distances are heavily saturated or estimated from too few comparable sites",
        ],
        upgma_assumptions=[
            "UPGMA assumes the pairwise distances are ultrametric and therefore consistent with a clock-like process",
            "UPGMA can mis-cluster taxa when rates vary among lineages even if the matrix remains symmetric and complete",
        ],
        ultrametric_compatible=not violations,
        ultrametric_tolerance=ultrametric_tolerance,
        upgma_ultrametric_violations=violations,
        warnings=warnings,
    )


def inspect_distance_matrix_quality(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceMatrixQualityReport:
    """Report saturation, outlier, and low-information risks for one computed matrix."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    saturation_report = diagnose_distance_saturation_from_genetic_distance_matrix(
        matrix
    )
    off_diagonal = [
        pair for pair in matrix.pairs if pair.left_identifier != pair.right_identifier
    ]
    defined_pairs = [pair for pair in off_diagonal if pair.distance is not None]
    saturated_pairs = saturated_pairs_from_diagnostics(saturation_report)

    high_distance_outliers: list[DistanceOutlierPair] = []
    if defined_pairs:
        distances = sorted(
            float(pair.distance) for pair in defined_pairs if pair.distance is not None
        )
        threshold = 0.75
        if len(distances) >= 4:
            q1 = distances[len(distances) // 4]
            q3 = distances[(3 * len(distances)) // 4]
            threshold = max(threshold, q3 + (1.5 * (q3 - q1)))
        for pair in defined_pairs:
            if pair.distance is not None and pair.distance >= threshold:
                high_distance_outliers.append(
                    DistanceOutlierPair(
                        left_identifier=pair.left_identifier,
                        right_identifier=pair.right_identifier,
                        distance=pair.distance,
                        note="pairwise distance is unusually large relative to the dataset baseline",
                    )
                )

    comparable_baseline = (
        median(pair.comparable_sites for pair in off_diagonal) if off_diagonal else 0
    )
    low_information_cutoff = max(10, int(math.floor(comparable_baseline * 0.5)))
    low_information_pairs = [
        LowInformationPair(
            left_identifier=pair.left_identifier,
            right_identifier=pair.right_identifier,
            comparable_sites=pair.comparable_sites,
            note="too few comparable sites remain for a stable distance estimate",
        )
        for pair in off_diagonal
        if pair.comparable_sites < low_information_cutoff
    ]
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )

    warnings: list[str] = []
    blocker_reasons: list[str] = []
    if saturated_pairs:
        warnings.append(
            "one or more pairwise distances are saturated or undefined under the selected model"
        )
        if len(saturated_pairs) / max(1, len(off_diagonal)) > 0.25:
            warnings.append(
                "many pairwise distances sit in a saturation regime that weakens distance-method assumptions"
            )
    if low_information_pairs:
        warnings.append("one or more taxon pairs retain too few comparable sites")
        if len(low_information_pairs) / max(1, len(off_diagonal)) > 0.5:
            warnings.append(
                "many taxon pairs retain too few comparable sites after filtering"
            )
    if not defined_pairs:
        blocker_reasons.append(
            "no off-diagonal distances could be computed from the selected alignment and policies"
        )
    if matrix.alignment_length < 10:
        warnings.append("alignment is short for robust distance-based tree building")
    if high_distance_outliers:
        warnings.append("one or more taxon pairs are unusually divergent")
    warnings.extend(assumptions.warnings)
    decision = "blocked" if blocker_reasons else ("risky" if warnings else "allowed")
    assessment_reasons = blocker_reasons if blocker_reasons else warnings
    return DistanceMatrixQualityReport(
        path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=matrix.inferred_alphabet,
        taxon_count=len(matrix.identifiers),
        pair_count=len(matrix.pairs),
        saturated_pairs=saturated_pairs,
        high_distance_outliers=high_distance_outliers,
        low_information_pairs=low_information_pairs,
        assumptions=assumptions,
        warnings=warnings,
        method_assessment=DistanceMethodAssessment(
            decision=decision,
            reasons=assessment_reasons,
        ),
    )
