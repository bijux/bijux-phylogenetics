from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import median

from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .missing_distance_policy import apply_missing_distance_policy
from .models import (
    DistanceMatrixAsymmetry,
    DistanceMethodAssumptionReport,
    ImportedDistanceEntry,
    ImportedDistanceMatrixQualityReport,
    ImportedDistanceMatrixReport,
    ImportedDistanceTreeBuildReport,
    LowInformationPair,
    MissingDistancePolicy,
    MissingDistancePolicyReport,
    NonMetricDistanceObservation,
    SaturatedDistancePair,
    UPGMAUltrametricViolation,
)
from .shared import (
    _build_distance_tree_from_lookup,
    _pair_key,
    _require_supported_distance_tree_method,
)
from .ultrametricity import diagnose_imported_distance_matrix_ultrametricity


def _distance_rows(entries: list[ImportedDistanceEntry]) -> list[ImportedDistanceEntry]:
    return [
        entry for entry in entries if entry.left_identifier != entry.right_identifier
    ]


def _unique_distance_rows(
    entries: list[ImportedDistanceEntry],
) -> list[ImportedDistanceEntry]:
    unique: dict[tuple[str, str], ImportedDistanceEntry] = {}
    for entry in _distance_rows(entries):
        unique.setdefault(
            _pair_key(entry.left_identifier, entry.right_identifier), entry
        )
    return [unique[key] for key in sorted(unique)]


def _imported_distance_scale(entries: list[ImportedDistanceEntry]) -> str:
    off_diagonal = _unique_distance_rows(entries)
    if off_diagonal and all(0.0 <= entry.distance <= 1.5 for entry in off_diagonal):
        return "unit-interval-like"
    return "unknown"


def load_imported_distance_matrix(path: Path) -> list[ImportedDistanceEntry]:
    """Load a long-form imported distance matrix table."""
    if not path.exists():
        raise FileNotFoundError(f"distance matrix file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"left_identifier", "right_identifier", "distance"}
        if reader.fieldnames is None or not required_columns <= set(reader.fieldnames):
            raise InvalidDistanceMatrixError(
                "distance matrix must contain left_identifier, right_identifier, and distance columns"
            )
        comparable_sites_column = (
            "comparable_sites" if "comparable_sites" in reader.fieldnames else None
        )
        entries: list[ImportedDistanceEntry] = []
        seen_directional_pairs: set[tuple[str, str]] = set()
        for row_index, row in enumerate(reader, start=2):
            left_identifier = str(row.get("left_identifier", "")).strip()
            right_identifier = str(row.get("right_identifier", "")).strip()
            raw_distance = str(row.get("distance", "")).strip()
            if not left_identifier or not right_identifier:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} is missing a left_identifier or right_identifier value"
                )
            if not raw_distance:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} is missing a distance value"
                )
            try:
                distance = float(raw_distance)
            except ValueError as error:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} contains a non-numeric distance value '{raw_distance}'"
                ) from error
            comparable_sites: int | None = None
            if comparable_sites_column is not None:
                raw_comparable_sites = str(row.get(comparable_sites_column, "")).strip()
                if raw_comparable_sites:
                    try:
                        comparable_sites = int(raw_comparable_sites)
                    except ValueError as error:
                        raise InvalidDistanceMatrixError(
                            f"row {row_index} in {path} contains a non-integer comparable_sites value '{raw_comparable_sites}'"
                        ) from error
            directional_pair = (left_identifier, right_identifier)
            if directional_pair in seen_directional_pairs:
                raise InvalidDistanceMatrixError(
                    f"distance matrix contains duplicate directional entry {left_identifier}/{right_identifier}"
                )
            seen_directional_pairs.add(directional_pair)
            entries.append(
                ImportedDistanceEntry(
                    left_identifier=left_identifier,
                    right_identifier=right_identifier,
                    distance=round(distance, 15),
                    comparable_sites=comparable_sites,
                )
            )

    if not entries:
        raise InvalidDistanceMatrixError(f"distance matrix contains no rows: {path}")
    return entries


def _symmetric_distances(
    entries: list[ImportedDistanceEntry],
) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}
    by_direction = {
        (entry.left_identifier, entry.right_identifier): entry.distance
        for entry in entries
    }
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    for left_identifier in identifiers:
        for right_identifier in identifiers:
            pair_key = _pair_key(left_identifier, right_identifier)
            if pair_key in distances:
                continue
            if left_identifier == right_identifier:
                if (left_identifier, right_identifier) in by_direction:
                    distances[pair_key] = by_direction[
                        (left_identifier, right_identifier)
                    ]
                continue
            if (left_identifier, right_identifier) in by_direction:
                distances[pair_key] = by_direction[(left_identifier, right_identifier)]
            elif (right_identifier, left_identifier) in by_direction:
                distances[pair_key] = by_direction[(right_identifier, left_identifier)]
    return distances


def assess_imported_distance_method_assumptions(
    path: Path,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit clock-like compatibility and core distance-tree assumptions for an imported matrix."""
    validation = validate_imported_distance_matrix(path)
    diagnostics = diagnose_imported_distance_matrix_ultrametricity(
        path,
        tolerance=ultrametric_tolerance,
    )
    distances = _symmetric_distances(load_imported_distance_matrix(path))
    warnings = list(diagnostics.warnings)
    if diagnostics.violating_triples:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=path,
        source_kind="imported-distance-matrix",
        taxon_count=len(validation.identifiers),
        pair_count=len(distances),
        nj_assumptions=[
            "neighbor-joining treats the matrix as an additive approximation and does not require a strict molecular clock",
            "neighbor-joining still inherits any errors or non-metric distortions present in the imported matrix",
        ],
        upgma_assumptions=[
            "UPGMA assumes the imported matrix is ultrametric and therefore compatible with a clock-like process",
            "UPGMA can enforce an ultrametric tree even when the source matrix does not satisfy that assumption",
        ],
        ultrametric_compatible=diagnostics.ultrametric,
        ultrametric_tolerance=ultrametric_tolerance,
        upgma_ultrametric_violations=_upgma_violations_from_diagnostics(diagnostics),
        warnings=warnings,
    )


def _upgma_violations_from_diagnostics(
    report,
) -> list[UPGMAUltrametricViolation]:
    return [
        UPGMAUltrametricViolation(
            left_identifier=row.left_identifier,
            middle_identifier=row.middle_identifier,
            right_identifier=row.right_identifier,
            smallest_distance=min(
                row.left_middle_distance,
                row.left_right_distance,
                row.middle_right_distance,
            ),
            middle_distance=row.second_largest_distance,
            largest_distance=row.largest_distance,
            deviation=row.violation,
        )
        for row in report.violating_triples
    ]


def validate_imported_distance_matrix(path: Path) -> ImportedDistanceMatrixReport:
    """Validate a long-form imported distance matrix table."""
    entries = load_imported_distance_matrix(path)
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    by_direction = {
        (entry.left_identifier, entry.right_identifier): entry for entry in entries
    }

    missing_pairs: list[str] = []
    diagonal_problems: list[str] = []
    negative_distance_pairs: list[str] = []
    asymmetric_pairs: list[DistanceMatrixAsymmetry] = []
    symmetric_distances = _symmetric_distances(entries)

    for left_identifier in identifiers:
        diagonal = by_direction.get((left_identifier, left_identifier))
        if diagonal is None:
            missing_pairs.append(f"{left_identifier}/{left_identifier}")
        elif diagonal.distance != 0.0:
            diagonal_problems.append(
                f"{left_identifier}/{left_identifier} has diagonal distance {diagonal.distance:g}"
            )
        if diagonal is not None and diagonal.distance < 0:
            negative_distance_pairs.append(f"{left_identifier}/{left_identifier}")

        for right_identifier in identifiers:
            if left_identifier >= right_identifier:
                continue
            left_to_right = by_direction.get((left_identifier, right_identifier))
            right_to_left = by_direction.get((right_identifier, left_identifier))
            if left_to_right is None and right_to_left is None:
                missing_pairs.append(f"{left_identifier}/{right_identifier}")
                continue
            if left_to_right is not None and left_to_right.distance < 0:
                negative_distance_pairs.append(f"{left_identifier}/{right_identifier}")
            if right_to_left is not None and right_to_left.distance < 0:
                negative_distance_pairs.append(f"{right_identifier}/{left_identifier}")
            if (
                left_to_right is not None
                and right_to_left is not None
                and not math.isclose(
                    left_to_right.distance,
                    right_to_left.distance,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            ):
                asymmetric_pairs.append(
                    DistanceMatrixAsymmetry(
                        left_identifier=left_identifier,
                        right_identifier=right_identifier,
                        left_to_right_distance=left_to_right.distance,
                        right_to_left_distance=right_to_left.distance,
                    )
                )

    nonmetric_observations: list[NonMetricDistanceObservation] = []
    if not missing_pairs and not diagonal_problems and not negative_distance_pairs:
        for left_index, left_identifier in enumerate(identifiers):
            for middle_index, middle_identifier in enumerate(identifiers):
                if middle_index == left_index:
                    continue
                for right_index, right_identifier in enumerate(identifiers):
                    if len({left_index, middle_index, right_index}) < 3:
                        continue
                    if left_identifier > right_identifier:
                        continue
                    direct_distance = symmetric_distances.get(
                        _pair_key(left_identifier, right_identifier)
                    )
                    left_middle = symmetric_distances.get(
                        _pair_key(left_identifier, middle_identifier)
                    )
                    middle_right = symmetric_distances.get(
                        _pair_key(middle_identifier, right_identifier)
                    )
                    if (
                        direct_distance is None
                        or left_middle is None
                        or middle_right is None
                    ):
                        continue
                    indirect_distance = round(left_middle + middle_right, 15)
                    if direct_distance > indirect_distance + 1e-12:
                        nonmetric_observations.append(
                            NonMetricDistanceObservation(
                                left_identifier=left_identifier,
                                middle_identifier=middle_identifier,
                                right_identifier=right_identifier,
                                direct_distance=direct_distance,
                                indirect_distance=indirect_distance,
                            )
                        )

    warnings: list[str] = []
    if missing_pairs:
        warnings.append("distance matrix is missing one or more required pairs")
    if diagonal_problems:
        warnings.append("distance matrix contains nonzero diagonal entries")
    if negative_distance_pairs:
        warnings.append("distance matrix contains negative distances")
    if asymmetric_pairs:
        warnings.append("distance matrix contains asymmetric directional entries")
    if nonmetric_observations:
        warnings.append(
            "distance matrix violates triangle inequality for one or more taxon triples"
        )

    return ImportedDistanceMatrixReport(
        path=path,
        identifiers=identifiers,
        pair_count=len(entries),
        complete=not missing_pairs,
        zero_diagonal=not diagonal_problems,
        symmetric=not asymmetric_pairs,
        nonnegative=not negative_distance_pairs,
        missing_pairs=missing_pairs,
        diagonal_problems=diagonal_problems,
        negative_distance_pairs=sorted(set(negative_distance_pairs)),
        asymmetric_pairs=sorted(
            asymmetric_pairs,
            key=lambda row: (row.left_identifier, row.right_identifier),
        ),
        nonmetric_observations=sorted(
            nonmetric_observations,
            key=lambda row: (
                row.left_identifier,
                row.middle_identifier,
                row.right_identifier,
            ),
        ),
        warnings=warnings,
    )


def _distance_lookup_from_imported(
    report: ImportedDistanceMatrixReport,
    entries: list[ImportedDistanceEntry],
    *,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> tuple[dict[tuple[str, str], float], MissingDistancePolicyReport]:
    if not report.zero_diagonal:
        raise InvalidDistanceMatrixError("distance matrix has nonzero diagonal entries")
    if not report.symmetric:
        raise InvalidDistanceMatrixError(
            "distance matrix contains asymmetric directional entries"
        )
    if not report.nonnegative:
        raise InvalidDistanceMatrixError("distance matrix contains negative distances")
    return apply_missing_distance_policy(
        report.identifiers,
        _symmetric_distances(entries),
        policy=missing_distance_policy,
    )


def inspect_imported_distance_matrix_quality(
    path: Path,
) -> ImportedDistanceMatrixQualityReport:
    """Report structural and heuristic quality risks for an imported distance matrix."""
    validation = validate_imported_distance_matrix(path)
    entries = load_imported_distance_matrix(path)
    off_diagonal = _unique_distance_rows(entries)
    distance_scale = _imported_distance_scale(entries)

    saturated_pairs: list[SaturatedDistancePair] = []
    if distance_scale == "unit-interval-like":
        saturated_pairs = [
            SaturatedDistancePair(
                left_identifier=entry.left_identifier,
                right_identifier=entry.right_identifier,
                distance=entry.distance,
                comparable_sites=entry.comparable_sites or 0,
                reason="imported distance lies in a high-divergence regime often associated with saturation for raw or lightly corrected genetic distances",
            )
            for entry in off_diagonal
            if entry.distance >= 0.75
        ]

    comparable_site_counts = [
        entry.comparable_sites
        for entry in off_diagonal
        if entry.comparable_sites is not None
    ]
    low_information_pair_cutoff: int | None = None
    low_information_pairs: list[LowInformationPair] = []
    if comparable_site_counts:
        comparable_baseline = median(comparable_site_counts)
        low_information_pair_cutoff = max(
            10, int(math.floor(comparable_baseline * 0.5))
        )
        low_information_pairs = [
            LowInformationPair(
                left_identifier=entry.left_identifier,
                right_identifier=entry.right_identifier,
                comparable_sites=int(entry.comparable_sites or 0),
                note="imported pair retains too few comparable sites for a stable distance interpretation",
            )
            for entry in off_diagonal
            if entry.comparable_sites is not None
            and entry.comparable_sites < low_information_pair_cutoff
        ]

    warnings = list(validation.warnings)
    if distance_scale != "unit-interval-like":
        warnings.append(
            "saturation heuristics were skipped because imported distances do not resemble unit-interval genetic distances"
        )
    elif saturated_pairs:
        warnings.append(
            "one or more imported distances fall in a high-divergence regime that deserves saturation review"
        )
    if comparable_site_counts:
        if low_information_pairs:
            warnings.append(
                "one or more imported distances retain too few comparable sites"
            )
    else:
        warnings.append(
            "low-information pair auditing is unavailable because the imported matrix does not provide comparable_sites"
        )

    return ImportedDistanceMatrixQualityReport(
        validation=validation,
        saturated_pairs=saturated_pairs,
        low_information_pairs=low_information_pairs,
        low_information_pair_cutoff=low_information_pair_cutoff,
        saturation_audit_scale=distance_scale,
        warnings=warnings,
    )


def build_tree_from_imported_distance_matrix(
    path: Path,
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> tuple[PhyloTree, ImportedDistanceTreeBuildReport]:
    """Build a distance-based tree from an imported long-form distance matrix."""
    method_policy = _require_supported_distance_tree_method(method)
    entries = load_imported_distance_matrix(path)
    validation = validate_imported_distance_matrix(path)
    assumptions = assess_imported_distance_method_assumptions(path)
    distance_lookup, missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
        missing_distance_policy=missing_distance_policy,
    )
    tree = _build_distance_tree_from_lookup(
        validation.identifiers,
        distance_lookup,
        method=method_policy.method,
    )
    return tree, ImportedDistanceTreeBuildReport(
        matrix_path=path,
        method=method_policy.method,
        method_policy=method_policy,
        taxon_count=len(validation.identifiers),
        pair_count=validation.pair_count,
        assumptions=assumptions,
        missing_distance_policy_report=missing_distance_policy_report,
    )
