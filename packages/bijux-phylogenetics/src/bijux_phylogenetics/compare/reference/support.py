from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology.support import compare_support_values
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class SupportReferenceObservation:
    """One governed support-mapping validation observation."""

    case_id: str
    category: str
    support_kind: str
    reference_source: str
    input_fixtures: list[Path]
    expected_metrics: dict[str, object]
    observed_metrics: dict[str, object]
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class SupportReferenceValidationReport:
    """Integrated support-mapping validation across governed reference cases."""

    case_count: int
    reference_case_count: int
    policy_case_count: int
    all_passed: bool
    observations: list[SupportReferenceObservation]


def validate_support_reference_examples(
    *,
    tolerance: float = 1e-12,
) -> SupportReferenceValidationReport:
    """Validate clade-mapped support parsing across governed support surfaces."""
    observations = [
        _validate_bootstrap_support_reference(tolerance=tolerance),
        _validate_sh_alrt_support_reference(tolerance=tolerance),
        _validate_fasttree_support_reference(tolerance=tolerance),
        _validate_posterior_frequency_reference(tolerance=tolerance),
        _validate_rotated_clade_mapping_case(tolerance=tolerance),
        _validate_bootstrap_posterior_topology_mismatch_case(tolerance=tolerance),
    ]
    reference_case_count = sum(
        observation.category == "reference-pair" for observation in observations
    )
    policy_case_count = sum(
        observation.category != "reference-pair" for observation in observations
    )
    return SupportReferenceValidationReport(
        case_count=len(observations),
        reference_case_count=reference_case_count,
        policy_case_count=policy_case_count,
        all_passed=all(observation.passed for observation in observations),
        observations=observations,
    )


def _validate_bootstrap_support_reference(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    from bijux_phylogenetics.engines.validation import (
        summarize_bootstrap_support_distribution,
    )

    tree_path = _trees_root() / "example_tree_support_left.nwk"
    tree = load_tree(tree_path)
    summary = summarize_bootstrap_support_distribution(tree_path)
    support_by_clade = _bootstrap_support_by_clade(
        summary,
        total_tip_count=tree.tip_count,
    )
    expected_metrics = {
        "support_by_clade": _load_numeric_reference_rows(
            "bootstrap_support_reference.tsv",
            metric_field="support",
        )
    }
    observed_metrics = {
        "support_by_clade": support_by_clade,
        "supported_branch_count": len(support_by_clade),
        "weakly_supported_clade_count": summary.weakly_supported_clade_count,
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    )
    return SupportReferenceObservation(
        case_id="iqtree-ufboot-reference",
        category="reference-pair",
        support_kind="ufboot",
        reference_source="checked iqtree bootstrap support tree fixture",
        input_fixtures=[tree_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else ["ufboot branch support no longer maps to the governed clades"],
    )


def _validate_sh_alrt_support_reference(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    from bijux_phylogenetics.engines.validation import (
        summarize_sh_alrt_support_distribution,
    )

    tree_path = _trees_root() / "example_tree_support_iqtree_composite.nwk"
    tree = load_tree(tree_path)
    summary = summarize_sh_alrt_support_distribution(tree_path)
    support_by_clade = _sh_alrt_support_by_clade(
        summary, total_tip_count=tree.tip_count
    )
    expected_metrics = {
        "support_by_clade": _load_sh_alrt_reference_rows(
            "sh_alrt_support_reference.tsv"
        )
    }
    observed_metrics = {
        "support_by_clade": support_by_clade,
        "annotated_branch_count": len(support_by_clade),
        "conflicting_support_signal_count": summary.conflicting_support_signal_count,
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    )
    return SupportReferenceObservation(
        case_id="iqtree-sh-alrt-reference",
        category="reference-pair",
        support_kind="sh-alrt-ufboot",
        reference_source="checked iqtree sh-alrt and ufboot support tree fixture",
        input_fixtures=[tree_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "compound sh-alrt or ufboot support no longer maps to the governed clades"
        ],
    )


def _validate_fasttree_support_reference(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    from bijux_phylogenetics.engines.validation import (
        summarize_fasttree_support_distribution,
    )

    tree_path = _trees_root() / "example_tree_support_fasttree_local.nwk"
    tree = load_tree(tree_path)
    summary = summarize_fasttree_support_distribution(tree_path)
    support_by_clade = _fasttree_support_by_clade(
        summary, total_tip_count=tree.tip_count
    )
    expected_metrics = {
        "support_by_clade": _load_numeric_reference_rows(
            "fasttree_support_reference.tsv",
            metric_field="local_support",
        )
    }
    observed_metrics = {
        "support_by_clade": support_by_clade,
        "annotated_branch_count": len(support_by_clade),
        "weakly_supported_clade_count": summary.weakly_supported_clade_count,
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    )
    return SupportReferenceObservation(
        case_id="fasttree-local-support-reference",
        category="reference-pair",
        support_kind="fasttree-local-support",
        reference_source="checked fasttree local support tree fixture",
        input_fixtures=[tree_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else ["fasttree local support no longer maps to the governed clades"],
    )


def _validate_posterior_frequency_reference(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    from bijux_phylogenetics.trees import compute_clade_frequency_table

    tree_set_path = _trees_root() / "example_tree_set_left.nwk"
    report = compute_clade_frequency_table(tree_set_path)
    expected_metrics = {
        "frequency_by_clade": _load_posterior_frequency_reference_rows(
            "posterior_clade_frequency_reference.tsv"
        )
    }
    observed_metrics = {
        "frequency_by_clade": {
            row.clade: {
                "tree_count": row.tree_count,
                "frequency": row.frequency,
            }
            for row in report.clade_frequencies
        },
        "tree_count": report.tree_count,
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    )
    return SupportReferenceObservation(
        case_id="posterior-clade-frequency-reference",
        category="reference-pair",
        support_kind="posterior-clade-frequency",
        reference_source="checked posterior tree-set clade frequency fixture",
        input_fixtures=[tree_set_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else ["posterior clade frequencies drifted from the governed tree-set fixture"],
    )


def _validate_rotated_clade_mapping_case(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    left_path = _trees_root() / "example_tree_support_iqtree_composite.nwk"
    right_path = _trees_root() / "example_tree_support_iqtree_composite_rotated.nwk"
    report = compare_support_values(left_path, right_path)
    expected_metrics = {
        "shared_clade_count": 2,
        "conflicting_clade_count": 0,
        "shared_support_pairs": {
            "A|B": {"left_support": 97.0, "right_support": 97.0},
            "C|D": {"left_support": 96.0, "right_support": 96.0},
        },
    }
    observed_metrics = {
        "shared_clade_count": len(report.shared_clades),
        "conflicting_clade_count": len(report.conflicting_clades),
        "shared_support_pairs": {
            row.split_id: {
                "left_support": row.left_support,
                "right_support": row.right_support,
                "support_fraction_delta": row.support_fraction_delta,
            }
            for row in report.shared_clades
        },
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    ) and all(row.support_fraction_delta == 0.0 for row in report.shared_clades)
    return SupportReferenceObservation(
        case_id="clade-mapped-support-rotation",
        category="clade-mapping-policy",
        support_kind="sh-alrt-ufboot",
        reference_source="owned clade-mapping policy",
        input_fixtures=[left_path, right_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "support comparison drifted toward node-order matching instead of clades"
        ],
    )


def _validate_bootstrap_posterior_topology_mismatch_case(
    *,
    tolerance: float,
) -> SupportReferenceObservation:
    from bijux_phylogenetics.trees import compare_bootstrap_and_posterior_uncertainty

    bootstrap_tree_path = _trees_root() / "example_tree_support_iqtree_composite.nwk"
    posterior_tree_set_path = _trees_root() / "example_tree_set_right.nwk"
    report = compare_bootstrap_and_posterior_uncertainty(
        bootstrap_tree_path,
        posterior_tree_set_path,
    )
    expected_metrics = {
        "topology_mismatch_detected": True,
        "topology_mismatch_clade_count": 2,
        "high_conflict_clade_count": 2,
        "row_agreements": {
            "A|B": "strong_conflict",
            "A|D": "method_specific",
            "B|C": "method_specific",
            "C|D": "strong_conflict",
        },
        "bootstrap_support_by_clade": {
            "A|B": 0.97,
            "C|D": 0.96,
        },
    }
    observed_metrics = {
        "topology_mismatch_detected": report.topology_mismatch_detected,
        "topology_mismatch_clade_count": report.topology_mismatch_clade_count,
        "high_conflict_clade_count": report.high_conflict_clade_count,
        "row_agreements": {row.clade: row.agreement for row in report.rows},
        "bootstrap_support_by_clade": {
            row.clade: row.bootstrap_support
            for row in report.rows
            if row.bootstrap_support is not None
        },
    }
    passed = _expected_metrics_match(
        expected_metrics,
        observed_metrics,
        tolerance=tolerance,
    )
    return SupportReferenceObservation(
        case_id="bootstrap-posterior-topology-mismatch",
        category="topology-mismatch-policy",
        support_kind="bootstrap-versus-posterior",
        reference_source="owned topology mismatch policy",
        input_fixtures=[bootstrap_tree_path, posterior_tree_set_path],
        expected_metrics=expected_metrics,
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "bootstrap versus posterior comparison no longer flags topology mismatches explicitly"
        ],
    )


def _bootstrap_support_by_clade(
    report: object,
    *,
    total_tip_count: int,
) -> dict[str, float]:
    return {
        "|".join(node.descendant_taxa): node.support
        for node in report.nodes
        if len(node.descendant_taxa) < total_tip_count
    }


def _sh_alrt_support_by_clade(
    report: object,
    *,
    total_tip_count: int,
) -> dict[str, dict[str, float | None]]:
    return {
        "|".join(node.descendant_taxa): {
            "sh_alrt_support": node.sh_alrt_support,
            "ufboot_support": node.ufboot_support,
        }
        for node in report.nodes
        if len(node.descendant_taxa) < total_tip_count
    }


def _fasttree_support_by_clade(
    report: object,
    *,
    total_tip_count: int,
) -> dict[str, float]:
    return {
        "|".join(node.descendant_taxa): node.local_support
        for node in report.nodes
        if len(node.descendant_taxa) < total_tip_count
    }


def _load_numeric_reference_rows(
    filename: str,
    *,
    metric_field: str,
) -> dict[str, float]:
    with (_metadata_root() / filename).open(encoding="utf-8", newline="") as handle:
        return {
            str(row["clade"]): float(row[metric_field])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _load_sh_alrt_reference_rows(
    filename: str,
) -> dict[str, dict[str, float]]:
    with (_metadata_root() / filename).open(encoding="utf-8", newline="") as handle:
        return {
            str(row["clade"]): {
                "sh_alrt_support": float(row["sh_alrt_support"]),
                "ufboot_support": float(row["ufboot_support"]),
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _load_posterior_frequency_reference_rows(
    filename: str,
) -> dict[str, dict[str, float | int]]:
    with (_metadata_root() / filename).open(encoding="utf-8", newline="") as handle:
        return {
            str(row["clade"]): {
                "tree_count": int(row["tree_count"]),
                "frequency": float(row["frequency"]),
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _expected_metrics_match(
    expected: dict[str, object],
    observed: dict[str, object],
    *,
    tolerance: float,
) -> bool:
    return _values_match(expected, observed, tolerance=tolerance)


def _values_match(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, float):
        return (
            isinstance(observed, (int, float))
            and abs(float(observed) - expected) <= tolerance
        )
    if isinstance(expected, dict):
        if not isinstance(observed, dict):
            return False
        if set(expected) - set(observed):
            return False
        return all(
            _values_match(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list):
        if not isinstance(observed, list) or len(expected) != len(observed):
            return False
        return all(
            _values_match(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return expected == observed


def _metadata_root() -> Path:
    return _package_root() / "tests" / "fixtures" / "metadata"


def _trees_root() -> Path:
    return _package_root() / "tests" / "fixtures" / "trees"


def _package_root() -> Path:
    return Path(__file__).resolve().parents[4]
