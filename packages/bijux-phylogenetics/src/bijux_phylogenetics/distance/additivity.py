from __future__ import annotations

from dataclasses import asdict
from itertools import combinations
import json
from pathlib import Path

from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .genetic_distance_matrix import (
    _build_alignment_distance_lookup,
    compute_pairwise_genetic_distance_matrix,
)
from .models import (
    AmbiguityPolicy,
    DistanceAdditivityDiagnosticsReport,
    DistanceFourPointViolation,
    DistanceModel,
    GapHandlingMode,
    GeneticDistanceMatrix,
)
from .shared import _pair_key


def diagnose_distance_additivity(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    tolerance: float = 1e-6,
) -> DistanceAdditivityDiagnosticsReport:
    """Test quartet additivity across one alignment-derived distance matrix."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    return diagnose_distance_additivity_from_genetic_distance_matrix(
        matrix,
        tolerance=tolerance,
    )


def diagnose_distance_additivity_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    tolerance: float = 1e-6,
) -> DistanceAdditivityDiagnosticsReport:
    """Test quartet additivity across one in-memory genetic distance matrix."""
    diagnostics = _build_additivity_diagnostics(
        identifiers=report.identifiers,
        pair_distances=_build_alignment_distance_lookup(report),
        source_path=report.path,
        source_kind="alignment",
        tolerance=tolerance,
        warnings=list(report.warnings),
    )
    if diagnostics.skipped_quartet_count > 0:
        diagnostics.warnings.append(
            "one or more quartets could not be tested because one or more pairwise distances were undefined under the selected model or policies"
        )
    return diagnostics


def diagnose_imported_distance_matrix_additivity(
    path: Path,
    *,
    tolerance: float = 1e-6,
) -> DistanceAdditivityDiagnosticsReport:
    """Test quartet additivity across one imported distance matrix."""
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
    diagnostics = _build_additivity_diagnostics(
        identifiers=validation.identifiers,
        pair_distances=_symmetric_distances(load_imported_distance_matrix(path)),
        source_path=path,
        source_kind="imported-distance-matrix",
        tolerance=tolerance,
        warnings=list(validation.warnings),
    )
    if diagnostics.skipped_quartet_count > 0:
        diagnostics.warnings.append(
            "one or more quartets could not be tested because the imported matrix omitted one or more required pairwise distances"
        )
    return diagnostics


def write_distance_additivity_table(
    path: Path,
    report: DistanceAdditivityDiagnosticsReport,
) -> Path:
    """Write the governed four-point violation ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "quartet",
        "split_ab_cd_sum",
        "split_ac_bd_sum",
        "split_ad_bc_sum",
        "best_split",
        "violation_magnitude",
    ]
    lines = ["\t".join(columns)]
    for row in report.violating_quartets:
        lines.append(
            "\t".join(
                [
                    _quartet_label(
                        row.first_identifier,
                        row.second_identifier,
                        row.third_identifier,
                        row.fourth_identifier,
                    ),
                    format(row.split_ab_cd_sum, ".12g"),
                    format(row.split_ac_bd_sum, ".12g"),
                    format(row.split_ad_bc_sum, ".12g"),
                    row.best_split,
                    format(row.violation_magnitude, ".12g"),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_additivity_run_json(
    path: Path,
    report: DistanceAdditivityDiagnosticsReport,
) -> Path:
    """Write the complete four-point additivity payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source_path": None if report.source_path is None else str(report.source_path),
        "source_kind": report.source_kind,
        "taxon_count": report.taxon_count,
        "defined_pair_count": report.defined_pair_count,
        "tested_quartet_count": report.tested_quartet_count,
        "skipped_quartet_count": report.skipped_quartet_count,
        "tolerance": report.tolerance,
        "additive": report.additive,
        "max_violation": report.max_violation,
        "violating_quartets": [asdict(row) for row in report.violating_quartets],
        "warnings": report.warnings,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_distance_additivity_artifacts(
    out_dir: Path,
    report: DistanceAdditivityDiagnosticsReport,
) -> dict[str, Path]:
    """Write governed quartet additivity artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    violation_table_path = write_distance_additivity_table(
        out_dir / "four_point_violations.tsv",
        report,
    )
    run_json_path = write_distance_additivity_run_json(out_dir / "run.json", report)
    return {
        "four_point_violations": violation_table_path,
        "run_json": run_json_path,
    }


def _build_additivity_diagnostics(
    *,
    identifiers: list[str],
    pair_distances: dict[tuple[str, str], float],
    source_path: Path | None,
    source_kind: str,
    tolerance: float,
    warnings: list[str],
) -> DistanceAdditivityDiagnosticsReport:
    violations: list[DistanceFourPointViolation] = []
    tested_quartet_count = 0
    skipped_quartet_count = 0
    max_violation = 0.0
    for (
        first_identifier,
        second_identifier,
        third_identifier,
        fourth_identifier,
    ) in combinations(
        identifiers,
        4,
    ):
        quartet_pairs = [
            _pair_key(first_identifier, second_identifier),
            _pair_key(first_identifier, third_identifier),
            _pair_key(first_identifier, fourth_identifier),
            _pair_key(second_identifier, third_identifier),
            _pair_key(second_identifier, fourth_identifier),
            _pair_key(third_identifier, fourth_identifier),
        ]
        if any(pair_key not in pair_distances for pair_key in quartet_pairs):
            skipped_quartet_count += 1
            continue
        tested_quartet_count += 1
        split_ab_cd_sum = (
            pair_distances[_pair_key(first_identifier, second_identifier)]
            + pair_distances[_pair_key(third_identifier, fourth_identifier)]
        )
        split_ac_bd_sum = (
            pair_distances[_pair_key(first_identifier, third_identifier)]
            + pair_distances[_pair_key(second_identifier, fourth_identifier)]
        )
        split_ad_bc_sum = (
            pair_distances[_pair_key(first_identifier, fourth_identifier)]
            + pair_distances[_pair_key(second_identifier, third_identifier)]
        )
        split_rows = [
            (
                f"{first_identifier},{second_identifier}|{third_identifier},{fourth_identifier}",
                split_ab_cd_sum,
            ),
            (
                f"{first_identifier},{third_identifier}|{second_identifier},{fourth_identifier}",
                split_ac_bd_sum,
            ),
            (
                f"{first_identifier},{fourth_identifier}|{second_identifier},{third_identifier}",
                split_ad_bc_sum,
            ),
        ]
        ordered_sums = sorted(value for _label, value in split_rows)
        violation_magnitude = abs(ordered_sums[2] - ordered_sums[1])
        max_violation = max(max_violation, violation_magnitude)
        if violation_magnitude > tolerance:
            best_split = min(split_rows, key=lambda item: (item[1], item[0]))[0]
            violations.append(
                DistanceFourPointViolation(
                    first_identifier=first_identifier,
                    second_identifier=second_identifier,
                    third_identifier=third_identifier,
                    fourth_identifier=fourth_identifier,
                    split_ab_cd_sum=split_ab_cd_sum,
                    split_ac_bd_sum=split_ac_bd_sum,
                    split_ad_bc_sum=split_ad_bc_sum,
                    best_split=best_split,
                    violation_magnitude=violation_magnitude,
                )
            )
    sorted_violations = sorted(
        violations,
        key=lambda row: (
            row.first_identifier,
            row.second_identifier,
            row.third_identifier,
            row.fourth_identifier,
        ),
    )
    return DistanceAdditivityDiagnosticsReport(
        source_path=source_path,
        source_kind=source_kind,
        taxon_count=len(identifiers),
        defined_pair_count=len(pair_distances),
        tested_quartet_count=tested_quartet_count,
        skipped_quartet_count=skipped_quartet_count,
        tolerance=tolerance,
        additive=not sorted_violations,
        max_violation=max_violation,
        violating_quartets=sorted_violations,
        warnings=warnings,
    )


def _quartet_label(
    first_identifier: str,
    second_identifier: str,
    third_identifier: str,
    fourth_identifier: str,
) -> str:
    return "|".join(
        [
            first_identifier,
            second_identifier,
            third_identifier,
            fourth_identifier,
        ]
    )
