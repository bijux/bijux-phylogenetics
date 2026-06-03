from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_tree_sets_structurally,
    compare_tree_structurally,
)
from bijux_phylogenetics.distance import build_tree_from_imported_distance_matrix
from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.phylo.topology import (
    extract_tree_clade_by_node_id,
    root_tree_on_outgroup,
    unroot_tree,
)
from bijux_phylogenetics.trees import (
    compute_consensus_tree,
    compute_strict_consensus_tree,
)

from ..registry import ApeParityCase
from .normalization import (
    _clear_branch_lengths,
    _compare_json,
    _load_json,
    _normalize_reference_summary,
    _normalize_tree_labels,
    _summary_rooted_flag,
)


def _tree_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    expected_tree = load_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _tree_set_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree_set = load_newick_tree_set(
        execution_root / "normalized-tree-set.nwk"
    )
    expected_tree_set = load_newick_tree_set(case.input_fixture)
    expected_tip_labels = {
        tip_name for tree in expected_tree_set for tip_name in tree.tip_names
    }
    for tree in reference_tree_set:
        _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    report = compare_tree_sets_structurally(
        expected_tree_set,
        reference_tree_set,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _root_tree_outgroup_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = True
    expected_tree, _report = root_tree_on_outgroup(
        case.input_fixture,
        outgroup_taxa=list(case.outgroup_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _unroot_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = unroot_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _drop_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = drop_tree_taxa(
        case.input_fixture, list(case.excluded_taxa)
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _keep_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = prune_tree_to_requested_taxa(
        case.input_fixture,
        list(case.requested_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _extract_clade_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    if case.node_id is None:
        raise ValueError(
            f"ape parity case '{case.case_id}' is missing an extraction node id"
        )
    expected_tree, _report = extract_tree_clade_by_node_id(
        case.input_fixture,
        node_id=case.node_id,
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _consensus_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    if case.consensus_method == "strict":
        expected_tree, _report = compute_strict_consensus_tree(case.input_fixture)
    elif case.consensus_method == "majority-rule":
        expected_tree, _report = compute_consensus_tree(case.input_fixture)
    else:
        raise ValueError(
            f"ape parity case '{case.case_id}' has unsupported consensus method "
            f"{case.consensus_method!r}"
        )
    expected_tree.rooted = False
    _clear_branch_lengths(expected_tree.root)
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    _clear_branch_lengths(reference_tree.root)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _neighbor_joining_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = build_tree_from_imported_distance_matrix(
        case.input_fixture,
        method="neighbor-joining",
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _supports_ard_rate_multiset_equivalence(
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> bool:
    if reference_summary is None or bijux_summary is None:
        return False
    return (
        reference_summary.get("model") == "all-rates-different"
        and bijux_summary.get("model") == "all-rates-different"
        and reference_summary.get("overparameterized") is True
        and bijux_summary.get("overparameterized") is True
    )


def _group_transition_rate_rows(
    rows: list[dict[str, object]],
) -> dict[tuple[bool, int], list[float]]:
    grouped: dict[tuple[bool, int], list[float]] = {}
    for row in rows:
        grouped.setdefault(
            (
                bool(row["transition_allowed"]),
                int(row["step_distance"]),
            ),
            [],
        ).append(float(row["rate"]))
    for values in grouped.values():
        values.sort()
    return grouped


def _transition_rate_rows_match(
    *,
    reference_rows: list[dict[str, object]],
    bijux_rows: list[dict[str, object]],
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    tolerance: float,
) -> bool:
    if _compare_json(reference_rows, bijux_rows, tolerance=tolerance):
        return True
    if not _supports_ard_rate_multiset_equivalence(
        reference_summary=reference_summary,
        bijux_summary=bijux_summary,
    ):
        return False
    return _compare_json(
        _group_transition_rate_rows(reference_rows),
        _group_transition_rate_rows(bijux_rows),
        tolerance=tolerance,
    )
