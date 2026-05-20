from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology.branch_lengths import (
    compare_branch_score_distance,
)
from bijux_phylogenetics.compare.topology.comparison import (
    compare_robinson_foulds,
)


@dataclass(slots=True)
class TreeDistanceReferenceObservation:
    """One governed validation observation for a tree-distance hard case."""

    case_id: str
    category: str
    method: str
    reference_tool: str
    reference_version: str
    reference_source: str
    rf_mode: str | None
    taxon_overlap_policy: str
    input_fixtures: list[Path]
    expected_metrics: dict[str, object]
    observed_metrics: dict[str, object]
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class TreeDistanceReferenceValidationReport:
    """Integrated tree-distance validation across governed hard cases."""

    case_count: int
    external_case_count: int
    policy_case_count: int
    all_passed: bool
    observations: list[TreeDistanceReferenceObservation]


def validate_tree_distance_reference_examples(
    *,
    tolerance: float = 1e-12,
) -> TreeDistanceReferenceValidationReport:
    """Validate RF and branch-score behavior against governed hard cases."""
    observations = [
        *_validate_robinson_foulds_reference_rows(tolerance=tolerance),
        *_validate_branch_score_reference_rows(tolerance=tolerance),
        _validate_identical_taxa_policy_for_robinson_foulds(),
        _validate_identical_taxa_policy_for_branch_score(),
    ]
    external_case_count = sum(
        observation.category == "external-reference" for observation in observations
    )
    policy_case_count = sum(
        observation.category == "overlap-policy" for observation in observations
    )
    return TreeDistanceReferenceValidationReport(
        case_count=len(observations),
        external_case_count=external_case_count,
        policy_case_count=policy_case_count,
        all_passed=all(observation.passed for observation in observations),
        observations=observations,
    )


def _validate_robinson_foulds_reference_rows(
    *,
    tolerance: float,
) -> list[TreeDistanceReferenceObservation]:
    observations: list[TreeDistanceReferenceObservation] = []
    for row in _load_robinson_foulds_reference_rows():
        input_paths = [
            _trees_root() / row["left_tree"],
            _trees_root() / row["right_tree"],
        ]
        report = compare_robinson_foulds(
            input_paths[0],
            input_paths[1],
            rf_mode=str(row["rf_mode"]),
            taxon_overlap_policy=str(row["taxon_overlap_policy"]),
        )
        expected_metrics = {
            "left_split_count": int(row["left_split_count"]),
            "right_split_count": int(row["right_split_count"]),
            "robinson_foulds_distance": int(row["robinson_foulds_distance"]),
            "normalized_robinson_foulds": float(row["normalized_robinson_foulds"]),
        }
        observed_metrics = {
            "left_split_count": report.left_split_count,
            "right_split_count": report.right_split_count,
            "robinson_foulds_distance": report.robinson_foulds_distance,
            "normalized_robinson_foulds": report.normalized_robinson_foulds,
            "shared_taxa": report.shared_taxa,
            "left_only_taxa": report.left_only_taxa,
            "right_only_taxa": report.right_only_taxa,
            "topology_equal": report.topology_equal,
        }
        passed = _expected_metrics_match(
            expected_metrics,
            observed_metrics,
            tolerance=tolerance,
        )
        observations.append(
            TreeDistanceReferenceObservation(
                case_id=str(row["case_id"]),
                category="external-reference",
                method="robinson-foulds-distance",
                reference_tool=str(row["reference_tool"]),
                reference_version=str(row["reference_version"]),
                reference_source=str(row["reference_source"]),
                rf_mode=str(row["rf_mode"]),
                taxon_overlap_policy=str(row["taxon_overlap_policy"]),
                input_fixtures=input_paths,
                expected_metrics=expected_metrics,
                observed_metrics=observed_metrics,
                passed=passed,
                notes=[]
                if passed
                else ["one or more RF metrics drifted from the governed reference"],
            )
        )
    return observations


def _validate_branch_score_reference_rows(
    *,
    tolerance: float,
) -> list[TreeDistanceReferenceObservation]:
    observations: list[TreeDistanceReferenceObservation] = []
    for row in _load_branch_score_reference_rows():
        input_paths = [
            _trees_root() / row["left_tree"],
            _trees_root() / row["right_tree"],
        ]
        report = compare_branch_score_distance(
            input_paths[0],
            input_paths[1],
            taxon_overlap_policy=str(row["taxon_overlap_policy"]),
        )
        expected_metrics = {
            "same_taxon_set": row["same_taxon_set"] == "true",
            "branch_score_distance": float(row["branch_score_distance"]),
        }
        observed_metrics = {
            "same_taxon_set": report.same_taxon_set,
            "branch_score_distance": report.branch_score_distance,
            "shared_taxa": report.shared_taxa,
            "left_only_taxa": report.left_only_taxa,
            "right_only_taxa": report.right_only_taxa,
            "split_count": report.split_count,
            "shared_split_count": report.shared_split_count,
            "left_only_split_count": report.left_only_split_count,
            "right_only_split_count": report.right_only_split_count,
            "missing_length_split_count": report.missing_length_split_count,
        }
        passed = _expected_metrics_match(
            expected_metrics,
            observed_metrics,
            tolerance=tolerance,
        )
        observations.append(
            TreeDistanceReferenceObservation(
                case_id=str(row["case_id"]),
                category="external-reference",
                method="branch-score-distance",
                reference_tool=str(row["reference_tool"]),
                reference_version=str(row["reference_version"]),
                reference_source=str(row["reference_source"]),
                rf_mode=None,
                taxon_overlap_policy=str(row["taxon_overlap_policy"]),
                input_fixtures=input_paths,
                expected_metrics=expected_metrics,
                observed_metrics=observed_metrics,
                passed=passed,
                notes=[]
                if passed
                else [
                    "branch-score distance drifted from the governed external reference"
                ],
            )
        )
    return observations


def _validate_identical_taxa_policy_for_robinson_foulds() -> (
    TreeDistanceReferenceObservation
):
    left_path = _trees_root() / "example_tree.nwk"
    right_path = _trees_root() / "example_tree_overlap.nwk"
    expected_metrics = {
        "raised": True,
        "error_substring": "requires identical taxon sets",
    }
    try:
        compare_robinson_foulds(
            left_path,
            right_path,
            taxon_overlap_policy="require-identical",
        )
    except ValueError as error:
        observed_metrics = {"raised": True, "error": str(error)}
        passed = expected_metrics["error_substring"] in str(error)
    else:  # pragma: no cover - defensive assertion
        observed_metrics = {"raised": False, "error": ""}
        passed = False
    return TreeDistanceReferenceObservation(
        case_id="rf-require-identical-overlap-policy",
        category="overlap-policy",
        method="robinson-foulds-distance",
        reference_tool="Bijux",
        reference_version="runtime",
        reference_source="owned overlap policy",
        rf_mode="rooted",
        taxon_overlap_policy="require-identical",
        input_fixtures=[left_path, right_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else ["require-identical RF comparison no longer rejects mismatched taxa"],
    )


def _validate_identical_taxa_policy_for_branch_score() -> (
    TreeDistanceReferenceObservation
):
    left_path = _trees_root() / "example_tree.nwk"
    right_path = _trees_root() / "example_tree_overlap.nwk"
    expected_metrics = {
        "raised": True,
        "error_substring": "identical taxon sets",
    }
    try:
        compare_branch_score_distance(
            left_path,
            right_path,
            taxon_overlap_policy="require-identical",
        )
    except ValueError as error:
        observed_metrics = {"raised": True, "error": str(error)}
        passed = expected_metrics["error_substring"] in str(error)
    else:  # pragma: no cover - defensive assertion
        observed_metrics = {"raised": False, "error": ""}
        passed = False
    return TreeDistanceReferenceObservation(
        case_id="branch-score-require-identical-overlap-policy",
        category="overlap-policy",
        method="branch-score-distance",
        reference_tool="Bijux",
        reference_version="runtime",
        reference_source="owned overlap policy",
        rf_mode=None,
        taxon_overlap_policy="require-identical",
        input_fixtures=[left_path, right_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "require-identical branch-score comparison no longer rejects mismatched taxa"
        ],
    )


def _expected_metrics_match(
    expected: dict[str, object],
    observed: dict[str, object],
    *,
    tolerance: float,
) -> bool:
    for key, expected_value in expected.items():
        if key not in observed:
            return False
        observed_value = observed[key]
        if isinstance(expected_value, float):
            if not isinstance(observed_value, (int, float)):
                return False
            if abs(float(observed_value) - expected_value) > tolerance:
                return False
            continue
        if observed_value != expected_value:
            return False
    return True


def _metadata_root() -> Path:
    return _package_root() / "tests" / "fixtures" / "metadata"


def _trees_root() -> Path:
    return _package_root() / "tests" / "fixtures" / "trees"


def _package_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _load_robinson_foulds_reference_rows() -> list[dict[str, str]]:
    with (_metadata_root() / "robinson_foulds_reference.tsv").open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_branch_score_reference_rows() -> list[dict[str, str]]:
    with (_metadata_root() / "branch_score_reference.tsv").open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
