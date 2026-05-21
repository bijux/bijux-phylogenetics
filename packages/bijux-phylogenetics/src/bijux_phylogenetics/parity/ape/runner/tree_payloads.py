from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.compare.topology import (
    compare_topology_distance_trees,
)
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    load_newick_tree_set,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.phylo.topology import (
    assess_tree_monophyly,
    extract_tree_clade_by_node_id,
    find_tree_mrca,
    root_tree_on_outgroup,
    unroot_tree,
)
from bijux_phylogenetics.phylo.topology.clades import informative_unrooted_splits
from bijux_phylogenetics.phylo.topology.tip_distances import (
    compute_tree_tip_distance_matrix,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_reference_tree_clade_support,
    compute_strict_consensus_tree,
    extract_tree_clades,
    extract_tree_set_clades,
)

from ..registry import ApeParityCase


def _node_kind_order(node_kind: str) -> int:
    return {"root": 0, "internal": 1, "tip": 2}.get(node_kind, 9)


def _sort_parity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            0 if row["tree_index"] == "" else int(row["tree_index"]),
            _node_kind_order(str(row["node_kind"])),
            str(row["clade_id"]),
            str(row["node_label"]),
        ),
    )


def _clade_rows_to_parity_rows(rows) -> list[dict[str, object]]:
    parity_rows = [
        {
            "tree_index": "" if row.tree_index is None else row.tree_index,
            "node_kind": row.node_kind,
            "clade_id": row.clade_id,
            "node_label": "" if row.node_label is None else row.node_label,
            "taxon_count": row.taxon_count,
            "taxa": "|".join(row.taxa),
            "support": "" if row.support is None else row.support,
            "branch_length": "" if row.branch_length is None else row.branch_length,
        }
        for row in rows
    ]
    return _sort_parity_rows(parity_rows)


def _build_bijux_tree_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    tree = load_tree(input_fixture)
    inspection = inspect_tree_path(input_fixture)
    clades = extract_tree_clades(input_fixture)
    return _tree_structure_payload(tree, inspection.rooted, clades.rows)


def _tree_structure_payload(
    tree: TreeNode | PhyloTree,
    rooted: bool | None,
    clade_rows,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    phylo_tree = (
        tree if isinstance(tree, PhyloTree) else PhyloTree(root=tree, rooted=rooted)
    )
    summary = {
        "tree_count": 1,
        "tip_count": phylo_tree.tip_count,
        "internal_node_count": phylo_tree.internal_node_count,
        "edge_count": phylo_tree.tip_count + phylo_tree.internal_node_count - 1,
        "rooted": rooted,
        "tip_labels": phylo_tree.tip_names,
        "branch_length_count": sum(
            1
            for branch_length in phylo_tree.branch_lengths()
            if branch_length is not None
        ),
    }
    return summary, _clade_rows_to_parity_rows(clade_rows), dumps_newick(phylo_tree)


def _build_bijux_tree_set_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], None]:
    clades = extract_tree_set_clades(input_fixture)
    parity_rows = _clade_rows_to_parity_rows(clades.rows)
    tree_indices = sorted(
        {
            tree_index
            for row in parity_rows
            for tree_index in [row.get("tree_index")]
            if isinstance(tree_index, int)
        }
    )
    first_tree_tip_labels = [
        node_label
        for row in parity_rows
        for node_label in [row.get("node_label")]
        if row.get("tree_index") == 1
        and row.get("node_kind") == "tip"
        and isinstance(node_label, str)
    ]
    summary = {
        "tree_count": clades.tree_count,
        "source_format": clades.source_format,
        "tree_indices": tree_indices,
        "shared_tip_labels": sorted(first_tree_tip_labels),
        "unique_tip_label_count": len(first_tree_tip_labels),
    }
    return summary, parity_rows, None


def _build_bijux_root_outgroup_structure(
    input_fixture: Path,
    *,
    outgroup_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    rooted_tree, _report = root_tree_on_outgroup(
        input_fixture,
        outgroup_taxa=list(outgroup_taxa),
    )
    with tempfile.TemporaryDirectory(prefix="bijux-ape-root-") as tmpdir:
        rooted_path = Path(tmpdir) / "rooted.nwk"
        write_newick(rooted_path, rooted_tree)
        clades = extract_tree_clades(rooted_path)
    return _tree_structure_payload(rooted_tree, rooted_tree.rooted, clades.rows)


def _build_bijux_unroot_structure(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    unrooted_tree, _report = unroot_tree(input_fixture)
    with tempfile.TemporaryDirectory(prefix="bijux-ape-unroot-") as tmpdir:
        unrooted_path = Path(tmpdir) / "unrooted.nwk"
        write_newick(unrooted_path, unrooted_tree)
        clades = extract_tree_clades(unrooted_path)
    return _tree_structure_payload(unrooted_tree, unrooted_tree.rooted, clades.rows)


def _build_bijux_drop_tip_structure(
    input_fixture: Path,
    *,
    excluded_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    pruned_tree, report = drop_tree_taxa(input_fixture, list(excluded_taxa))
    with tempfile.TemporaryDirectory(prefix="bijux-ape-drop-tip-") as tmpdir:
        pruned_path = Path(tmpdir) / "drop-tip.nwk"
        write_newick(pruned_path, pruned_tree)
        clades = extract_tree_clades(pruned_path)
    summary, rows, normalized_text = _tree_structure_payload(
        pruned_tree,
        pruned_tree.rooted,
        clades.rows,
    )
    summary["dropped_taxa"] = report.removed_taxa
    summary["absent_requested_taxa"] = report.absent_requested_taxa
    return summary, rows, normalized_text


def _build_bijux_keep_tip_structure(
    input_fixture: Path,
    *,
    requested_taxa: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    pruned_tree, report = prune_tree_to_requested_taxa(
        input_fixture, list(requested_taxa)
    )
    with tempfile.TemporaryDirectory(prefix="bijux-ape-keep-tip-") as tmpdir:
        pruned_path = Path(tmpdir) / "keep-tip.nwk"
        write_newick(pruned_path, pruned_tree)
        clades = extract_tree_clades(pruned_path)
    summary, rows, normalized_text = _tree_structure_payload(
        pruned_tree,
        pruned_tree.rooted,
        clades.rows,
    )
    summary["requested_taxa"] = report.requested_taxa
    summary["dropped_taxa"] = report.removed_taxa
    return summary, rows, normalized_text


def _build_bijux_extract_clade_structure(
    input_fixture: Path,
    *,
    node_id: int,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    subtree, report = extract_tree_clade_by_node_id(input_fixture, node_id=node_id)
    with tempfile.TemporaryDirectory(prefix="bijux-ape-extract-clade-") as tmpdir:
        subtree_path = Path(tmpdir) / "extract-clade.nwk"
        write_newick(subtree_path, subtree)
        clades = extract_tree_clades(subtree_path)
    summary, rows, normalized_text = _tree_structure_payload(
        subtree,
        subtree.rooted,
        clades.rows,
    )
    summary["requested_node_id"] = report.requested_node_id
    summary["matched_node_id"] = report.matched_node_id
    summary["matched_node_name"] = report.matched_node_name
    return summary, rows, normalized_text


def _build_bijux_mrca_summary(
    input_fixture: Path,
    *,
    mrca_taxa: tuple[str, ...],
) -> dict[str, object]:
    report = find_tree_mrca(input_fixture, taxa=list(mrca_taxa))
    return {
        "requested_taxa": report.requested_taxa,
        "unique_requested_taxa": report.unique_requested_taxa,
        "duplicate_requested_taxa": report.duplicate_requested_taxa,
        "matched_node_id": report.matched_node_id,
        "matched_node_name": report.matched_node_name or "",
        "matched_taxa": report.matched_taxa,
        "matched_extra_taxa": report.matched_extra_taxa,
        "matched_tip_count": report.matched_tip_count,
        "is_root": report.is_root,
    }


def _build_bijux_monophyly_summary(
    input_fixture: Path,
    *,
    requested_taxa: tuple[str, ...],
    reroot: bool,
) -> dict[str, object]:
    report = assess_tree_monophyly(
        input_fixture,
        taxa=list(requested_taxa),
        reroot=reroot,
    )
    return {
        "requested_taxa": report.requested_taxa,
        "unique_requested_taxa": report.unique_requested_taxa,
        "duplicate_requested_taxa": report.duplicate_requested_taxa,
        "missing_requested_taxa": report.missing_requested_taxa,
        "present_requested_taxa": report.present_requested_taxa,
        "reroot": report.reroot,
        "rooted": report.rooted,
        "monophyletic": report.monophyletic,
        "complementary_clade_used": report.complementary_clade_used,
        "matched_node_id": report.matched_node_id,
        "matched_node_name": report.matched_node_name or "",
        "matched_taxa": report.matched_taxa,
        "matched_extra_taxa": report.matched_extra_taxa,
        "matched_tip_count": report.matched_tip_count,
        "is_root": report.is_root,
    }


def _materialize_reference_input(case: ApeParityCase, working_root: Path) -> Path:
    reference_input_path = working_root / "bijux-reference-input.nwk"
    if case.operation == "write-tree-structure":
        tree = load_tree(case.input_fixture)
        write_newick(reference_input_path, tree)
        return reference_input_path
    if case.operation == "write-tree-set-structure":
        trees = load_newick_tree_set(case.input_fixture)
        write_newick_tree_set(reference_input_path, trees)
        return reference_input_path
    return case.input_fixture


def _build_bijux_tree_tip_distance_rows(
    input_fixture: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_tree_tip_distance_matrix(input_fixture)
    return {
        "tip_count": len(report.identifiers),
        "rooted": report.rooted,
        "tip_labels": report.identifiers,
        "pair_count": report.pair_count,
        "diagonal_zero": report.diagonal_zero,
        "symmetric": report.symmetric,
        "complete_branch_lengths": report.complete_branch_lengths,
        "missing_branch_length_policy": report.missing_branch_length_policy,
    }, [
        {
            "left_identifier": row.left_identifier,
            "right_identifier": row.right_identifier,
            "distance": row.distance,
        }
        for row in report.pairs
    ]


def _build_bijux_consensus_rows(
    input_fixture: Path,
    *,
    consensus_method: str,
) -> tuple[dict[str, object], list[dict[str, object]], str]:
    if consensus_method == "strict":
        tree, report = compute_strict_consensus_tree(input_fixture)
    elif consensus_method == "majority-rule":
        tree, report = compute_consensus_tree(input_fixture)
    else:
        raise ValueError(f"unsupported consensus method '{consensus_method}'")
    tree.rooted = False
    frequency_report = compute_clade_frequency_table(input_fixture)
    included_split_count = len(informative_unrooted_splits(tree, set(tree.tip_names)))
    return (
        {
            "tree_count": report.tree_count,
            "shared_taxa": report.shared_taxa,
            "shared_taxon_count": len(report.shared_taxa),
            "tip_count": len(tree.tip_names),
            "rooted": False,
            "consensus_method": report.consensus_method,
            "consensus_threshold": report.consensus_threshold,
            "included_clade_count": included_split_count,
            "clade_frequency_count": len(frequency_report.clade_frequencies),
        },
        [
            {
                "clade": row.clade,
                "tree_count": row.tree_count,
                "frequency": row.frequency,
            }
            for row in frequency_report.clade_frequencies
        ],
        dumps_newick(tree),
    )


def _build_bijux_prop_clades_rows(
    reference_tree_path: Path,
    comparison_tree_set_path: Path,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    report = compute_reference_tree_clade_support(
        reference_tree_path,
        comparison_tree_set_path,
    )
    return {
        "tree_count": report.tree_count,
        "shared_taxa": report.shared_taxa,
        "shared_taxon_count": len(report.shared_taxa),
        "internal_node_count": len(report.rows),
        "supported_clade_count": report.supported_clade_count,
        "absent_clade_count": report.absent_clade_count,
        "unscored_clade_count": report.unscored_clade_count,
    }, [
        {
            "node_id": row.node_id,
            "node_kind": row.node_kind,
            "node_label": "" if row.node_label is None else row.node_label,
            "descendant_taxa": "|".join(row.descendant_taxa),
            "supporting_tree_count": (
                "" if row.supporting_tree_count is None else row.supporting_tree_count
            ),
            "clade_frequency": ""
            if row.clade_frequency is None
            else row.clade_frequency,
            "support_percent": ""
            if row.support_percent is None
            else row.support_percent,
            "support_status": row.support_status,
            "explanation": row.explanation,
            "reference_branch_length": (
                ""
                if row.reference_branch_length is None
                else row.reference_branch_length
            ),
            "reference_root_depth": (
                "" if row.reference_root_depth is None else row.reference_root_depth
            ),
        }
        for row in report.rows
    ]


def _inspect_tree_set_rooted_flags(input_fixture: Path) -> tuple[bool, bool]:
    lines = [
        line.strip()
        for line in input_fixture.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) != 2:
        raise ValueError(
            "ape topology-distance parity fixtures must contain exactly two trees"
        )
    rooted_flags: list[bool] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-topology-distance-rootedness-"
    ) as tmpdir:
        temporary_root = Path(tmpdir)
        for index, line in enumerate(lines, start=1):
            temporary_path = temporary_root / f"tree-{index}.nwk"
            temporary_path.write_text(f"{line}\n", encoding="utf-8")
            rooted_flags.append(inspect_tree_path(temporary_path).rooted)
    return rooted_flags[0], rooted_flags[1]


def _build_bijux_topology_distance_rows(
    input_fixture: Path,
    *,
    rf_mode: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_set = load_newick_tree_set(input_fixture)
    if len(tree_set) != 2:
        raise ValueError(
            "ape topology-distance parity fixtures must contain exactly two trees"
        )
    report = compare_topology_distance_trees(
        tree_set[0],
        tree_set[1],
        left_path=input_fixture,
        right_path=input_fixture,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
    )
    rooted_left, rooted_right = _inspect_tree_set_rooted_flags(input_fixture)
    return {
        "tip_count": len(report.shared_taxa),
        "shared_taxa": report.shared_taxa,
        "left_only_taxa": report.left_only_taxa,
        "right_only_taxa": report.right_only_taxa,
        "taxon_overlap_policy": report.taxon_overlap_policy,
        "rf_mode": report.rf_mode,
        "rooted_left": rooted_left,
        "rooted_right": rooted_right,
        "polytomy_present_left": report.polytomy_present_left,
        "polytomy_present_right": report.polytomy_present_right,
        "left_split_count": report.left_split_count,
        "right_split_count": report.right_split_count,
        "shared_split_count": report.shared_split_count,
        "left_only_split_count": report.left_only_split_count,
        "right_only_split_count": report.right_only_split_count,
        "robinson_foulds_distance": report.robinson_foulds_distance,
        "normalized_robinson_foulds": report.normalized_robinson_foulds,
        "topology_equal": report.topology_equal,
    }, [
        {
            "split_id": row.split_id,
            "split_kind": row.split_kind,
            "comparison_status": row.comparison_status,
            "taxon_count": row.taxon_count,
            "descendant_taxa": "|".join(row.descendant_taxa),
            "left_present": row.left_present,
            "right_present": row.right_present,
        }
        for row in report.split_rows
    ]
