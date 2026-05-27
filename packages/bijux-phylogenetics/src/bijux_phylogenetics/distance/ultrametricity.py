from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .genetic_distance_matrix import (
    _build_alignment_distance_lookup,
    compute_pairwise_genetic_distance_matrix,
)
from .models import (
    AmbiguityPolicy,
    DistanceModel,
    DistanceUltrametricityDiagnosticsReport,
    DistanceUltrametricityViolation,
    GapHandlingMode,
    GeneticDistanceMatrix,
)
from .shared import _pair_key


def diagnose_distance_ultrametricity(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    tolerance: float = 1e-6,
) -> DistanceUltrametricityDiagnosticsReport:
    """Test the three-point condition across all comparable alignment taxon triples."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    return diagnose_distance_ultrametricity_from_genetic_distance_matrix(
        matrix,
        tolerance=tolerance,
    )


def diagnose_distance_ultrametricity_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    tolerance: float = 1e-6,
) -> DistanceUltrametricityDiagnosticsReport:
    """Test one in-memory genetic distance matrix for ultrametric three-point violations."""
    warnings = list(report.warnings)
    diagnostics = _build_ultrametricity_diagnostics(
        identifiers=report.identifiers,
        pair_distances=_build_alignment_distance_lookup(report),
        source_path=report.path,
        source_kind="alignment",
        tolerance=tolerance,
        warnings=warnings,
    )
    if diagnostics.skipped_triple_count > 0:
        diagnostics.warnings.append(
            "one or more taxon triples could not be tested because one or more pairwise distances were undefined under the selected model or policies"
        )
    return diagnostics


def diagnose_imported_distance_matrix_ultrametricity(
    path: Path,
    *,
    tolerance: float = 1e-6,
) -> DistanceUltrametricityDiagnosticsReport:
    """Test one imported distance matrix for ultrametric three-point violations."""
    from .imported import (
        _symmetric_distances,
        load_imported_distance_matrix,
        validate_imported_distance_matrix,
    )

    validation = validate_imported_distance_matrix(path)
    if not validation.zero_diagonal:
        raise InvalidDistanceMatrixError("distance matrix has nonzero diagonal entries")
    if not validation.symmetric:
        raise InvalidDistanceMatrixError(
            "distance matrix contains asymmetric directional entries"
        )
    if not validation.nonnegative:
        raise InvalidDistanceMatrixError("distance matrix contains negative distances")
    diagnostics = _build_ultrametricity_diagnostics(
        identifiers=validation.identifiers,
        pair_distances=_symmetric_distances(load_imported_distance_matrix(path)),
        source_path=path,
        source_kind="imported-distance-matrix",
        tolerance=tolerance,
        warnings=list(validation.warnings),
    )
    if diagnostics.skipped_triple_count > 0:
        diagnostics.warnings.append(
            "one or more taxon triples could not be tested because the imported matrix omitted one or more required pairwise distances"
        )
    return diagnostics


def _build_ultrametricity_diagnostics(
    *,
    identifiers: list[str],
    pair_distances: dict[tuple[str, str], float],
    source_path: Path | None,
    source_kind: str,
    tolerance: float,
    warnings: list[str],
) -> DistanceUltrametricityDiagnosticsReport:
    violations: list[DistanceUltrametricityViolation] = []
    tested_triple_count = 0
    skipped_triple_count = 0
    max_violation = 0.0
    for left_index, left_identifier in enumerate(identifiers):
        for middle_index in range(left_index + 1, len(identifiers)):
            middle_identifier = identifiers[middle_index]
            for right_index in range(middle_index + 1, len(identifiers)):
                right_identifier = identifiers[right_index]
                left_middle_key = _pair_key(left_identifier, middle_identifier)
                left_right_key = _pair_key(left_identifier, right_identifier)
                middle_right_key = _pair_key(middle_identifier, right_identifier)
                if (
                    left_middle_key not in pair_distances
                    or left_right_key not in pair_distances
                    or middle_right_key not in pair_distances
                ):
                    skipped_triple_count += 1
                    continue
                tested_triple_count += 1
                left_middle_distance = pair_distances[left_middle_key]
                left_right_distance = pair_distances[left_right_key]
                middle_right_distance = pair_distances[middle_right_key]
                ordered = sorted(
                    (
                        left_middle_distance,
                        left_right_distance,
                        middle_right_distance,
                    )
                )
                violation = abs(ordered[2] - ordered[1])
                max_violation = max(max_violation, violation)
                if violation > tolerance:
                    violations.append(
                        DistanceUltrametricityViolation(
                            left_identifier=left_identifier,
                            middle_identifier=middle_identifier,
                            right_identifier=right_identifier,
                            left_middle_distance=left_middle_distance,
                            left_right_distance=left_right_distance,
                            middle_right_distance=middle_right_distance,
                            second_largest_distance=ordered[1],
                            largest_distance=ordered[2],
                            violation=violation,
                        )
                    )
    sorted_violations = sorted(
        violations,
        key=lambda row: (
            row.left_identifier,
            row.middle_identifier,
            row.right_identifier,
        ),
    )
    return DistanceUltrametricityDiagnosticsReport(
        source_path=source_path,
        source_kind=source_kind,
        taxon_count=len(identifiers),
        defined_pair_count=len(pair_distances),
        tested_triple_count=tested_triple_count,
        skipped_triple_count=skipped_triple_count,
        tolerance=tolerance,
        ultrametric=not sorted_violations,
        max_violation=max_violation,
        violating_triples=sorted_violations,
        warnings=warnings,
    )
