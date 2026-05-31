from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    rooted_topology_fingerprint,
    rooted_topology_signature_ids,
)
from bijux_phylogenetics.phylo.topology.models import (
    RootedSprEnumerationBudget,
    RootedSprMoveApplicationReport,
    RootedSprNeighborRow,
    RootedSprNeighborhoodReport,
)

from .affected_subtrees import summarize_affected_subtrees
from .tree import PhyloTree, TreeNode, descendant_taxa

_ROOT_REGRAFT_BRANCH_ID = "root"


@dataclass(frozen=True, slots=True)
class RootedSprMoveCandidate:
    """One deterministic rooted SPR move over a binary rooted tree."""

    pruned_node_id: str
    pruned_clade_id: str
    pruned_descendant_taxa: tuple[str, ...]
    regraft_target_branch_id: str
    regraft_target_descendant_taxa: tuple[str, ...] | None


def rooted_spr_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort rooted SPR candidates deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def rooted_spr_clade_id(node: TreeNode) -> str:
    """Render one rooted SPR branch endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_rooted_spr_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a rooted SPR tree."""
    if node.node_id is None:
        raise AssertionError("rooted SPR search requires refreshed node identities")
    return node.node_id


def validate_rooted_spr_enumeration_budget(
    budget: RootedSprEnumerationBudget | None,
) -> RootedSprEnumerationBudget:
    """Normalize one explicit rooted SPR enumeration budget."""
    if budget is None:
        return RootedSprEnumerationBudget()
    if (
        budget.max_pruned_clade_count is not None
        and budget.max_pruned_clade_count <= 0
    ):
        raise ValueError("rooted SPR prune-node budget must be positive when provided")
    if (
        budget.max_regraft_target_count_per_pruned_clade is not None
        and budget.max_regraft_target_count_per_pruned_clade <= 0
    ):
        raise ValueError(
            "rooted SPR regraft-target budget must be positive when provided"
        )
    return budget


def iter_rooted_spr_move_candidates(
    tree: PhyloTree,
    *,
    budget: RootedSprEnumerationBudget | None = None,
):
    """Yield deterministic rooted SPR candidates for one rooted tree."""
    candidates, _skipped_pruned, _skipped_regraft_targets = (
        _collect_rooted_spr_move_candidates(tree, budget=budget)
    )
    yield from candidates


def _collect_rooted_spr_move_candidates(
    tree: PhyloTree,
    *,
    budget: RootedSprEnumerationBudget | None,
) -> tuple[list[RootedSprMoveCandidate], int, int]:
    normalized_budget = validate_rooted_spr_enumeration_budget(budget)
    seen_candidates: set[tuple[str, str]] = set()
    sorted_prune_nodes = sorted(
        (node for node in tree.iter_nodes(order="preorder") if node is not tree.root),
        key=rooted_spr_node_sort_key,
    )
    if normalized_budget.max_pruned_clade_count is None:
        limited_prune_nodes = sorted_prune_nodes
    else:
        limited_prune_nodes = sorted_prune_nodes[
            : normalized_budget.max_pruned_clade_count
        ]
    skipped_pruned_clade_count = len(sorted_prune_nodes) - len(limited_prune_nodes)
    skipped_regraft_target_count = 0
    candidates: list[RootedSprMoveCandidate] = []
    for prune_node in limited_prune_nodes:
        pruned_clade_id = rooted_spr_clade_id(prune_node)
        remainder_tree, _pruned_subtree = prune_rooted_spr_subtree(
            tree,
            require_rooted_spr_node_id(prune_node),
        )
        regraft_targets = [
            (_ROOT_REGRAFT_BRANCH_ID, None),
            *[
                (
                    rooted_spr_clade_id(target_node),
                    tuple(sorted(target_node.descendant_taxa)),
                )
                for target_node in sorted(
                    remainder_tree.iter_nodes(order="preorder"),
                    key=rooted_spr_node_sort_key,
                )
            ],
        ]
        if normalized_budget.max_regraft_target_count_per_pruned_clade is not None:
            regraft_target_limit = (
                normalized_budget.max_regraft_target_count_per_pruned_clade
            )
            skipped_regraft_target_count += max(
                0,
                len(regraft_targets) - regraft_target_limit,
            )
            regraft_targets = regraft_targets[:regraft_target_limit]
        for regraft_target_branch_id, regraft_target_descendant_taxa in regraft_targets:
            signature = (pruned_clade_id, regraft_target_branch_id)
            if signature in seen_candidates:
                continue
            seen_candidates.add(signature)
            candidates.append(
                RootedSprMoveCandidate(
                    pruned_node_id=require_rooted_spr_node_id(prune_node),
                    pruned_clade_id=pruned_clade_id,
                    pruned_descendant_taxa=tuple(sorted(prune_node.descendant_taxa)),
                    regraft_target_branch_id=regraft_target_branch_id,
                    regraft_target_descendant_taxa=regraft_target_descendant_taxa,
                )
            )
    return candidates, skipped_pruned_clade_count, skipped_regraft_target_count


def apply_rooted_spr_move(
    tree: PhyloTree,
    candidate: RootedSprMoveCandidate,
) -> PhyloTree:
    """Return one copied rooted tree with the selected SPR move applied."""
    tree_has_explicit_branch_lengths = any(
        child.branch_length is not None for _parent, child in tree.iter_edges()
    )
    remainder_tree, pruned_subtree = prune_rooted_spr_subtree(
        tree,
        candidate.pruned_node_id,
    )
    if candidate.regraft_target_branch_id == _ROOT_REGRAFT_BRANCH_ID:
        if tree_has_explicit_branch_lengths:
            if remainder_tree.root.branch_length is None:
                remainder_tree.root.branch_length = _seed_missing_rooted_spr_branch_length(
                    pruned_subtree
                )
            if pruned_subtree.branch_length is None:
                pruned_subtree.branch_length = _seed_missing_rooted_spr_branch_length(
                    remainder_tree.root
                )
        return PhyloTree(
            root=TreeNode(children=[remainder_tree.root, pruned_subtree]),
            source_format=remainder_tree.source_format,
            rooted=True,
        ).refresh()
    regraft_target = _find_node_by_clade_id(
        remainder_tree,
        candidate.regraft_target_branch_id,
    )
    target_parent = regraft_target.parent
    original_target_branch_length = regraft_target.branch_length
    inserted = TreeNode(children=[regraft_target, pruned_subtree])
    if tree_has_explicit_branch_lengths:
        if original_target_branch_length is not None:
            inserted.branch_length = original_target_branch_length / 2.0
            regraft_target.branch_length = original_target_branch_length / 2.0
        else:
            seeded_length = _seed_missing_rooted_spr_branch_length(pruned_subtree)
            inserted.branch_length = seeded_length
            regraft_target.branch_length = seeded_length
    if target_parent is None:
        return PhyloTree(
            root=inserted,
            source_format=remainder_tree.source_format,
            rooted=True,
        ).refresh()
    target_parent.replace_children(
        [
            child if child is not regraft_target else inserted
            for child in target_parent.children
        ]
    )
    return remainder_tree.refresh()


def resolve_rooted_spr_move_candidate(
    tree: PhyloTree | Path,
    move_index: int,
    *,
    budget: RootedSprEnumerationBudget | None = None,
) -> tuple[RootedSprMoveCandidate, int]:
    """Resolve one 1-indexed rooted SPR move candidate from a validated tree."""
    resolved_tree, _input_tree_path = _resolve_rooted_spr_tree(tree)
    validate_rooted_spr_tree(resolved_tree)
    candidates = list(iter_rooted_spr_move_candidates(resolved_tree, budget=budget))
    if move_index < 1 or move_index > len(candidates):
        raise ValueError(
            f"rooted SPR move index must be between 1 and {len(candidates)}"
        )
    return candidates[move_index - 1], len(candidates)


def summarize_rooted_spr_move_application(
    tree: PhyloTree | Path,
    move_index: int,
    *,
    budget: RootedSprEnumerationBudget | None = None,
) -> RootedSprMoveApplicationReport:
    """Apply one rooted SPR move and report the resulting valid transformed tree."""
    resolved_tree, input_tree_path = _resolve_rooted_spr_tree(tree)
    validate_rooted_spr_tree(resolved_tree)
    normalized_budget = validate_rooted_spr_enumeration_budget(budget)
    selected_candidate, available_move_count = resolve_rooted_spr_move_candidate(
        resolved_tree,
        move_index,
        budget=normalized_budget,
    )
    moved_tree = apply_rooted_spr_move(resolved_tree, selected_candidate)
    input_tip_taxa = set(resolved_tree.tip_names)
    moved_tip_taxa = set(moved_tree.tip_names)
    input_topology_fingerprint = rooted_topology_fingerprint(resolved_tree)
    moved_topology_fingerprint = rooted_topology_fingerprint(moved_tree)
    affected_clade_ids = _rooted_spr_affected_clade_ids(
        resolved_tree,
        moved_tree,
        selected_candidate,
    )
    report = RootedSprMoveApplicationReport(
        algorithm="rooted-spr-move-application",
        input_tree_path=input_tree_path,
        input_tree_newick=resolved_tree.to_newick(),
        input_topology_fingerprint=input_topology_fingerprint,
        selected_move_index=move_index,
        available_move_count=available_move_count,
        max_pruned_clade_count=normalized_budget.max_pruned_clade_count,
        max_regraft_target_count_per_pruned_clade=(
            normalized_budget.max_regraft_target_count_per_pruned_clade
        ),
        selected_pruned_node_id=selected_candidate.pruned_node_id,
        selected_pruned_clade_id=selected_candidate.pruned_clade_id,
        selected_pruned_descendant_taxa=list(selected_candidate.pruned_descendant_taxa),
        selected_regraft_target_branch_id=selected_candidate.regraft_target_branch_id,
        selected_regraft_target_descendant_taxa=(
            None
            if selected_candidate.regraft_target_descendant_taxa is None
            else list(selected_candidate.regraft_target_descendant_taxa)
        ),
        moved_tree_newick=moved_tree.to_newick(),
        moved_topology_fingerprint=moved_topology_fingerprint,
        moved_topology_changed=moved_topology_fingerprint != input_topology_fingerprint,
        tip_count=resolved_tree.tip_count,
        internal_node_count=resolved_tree.internal_node_count,
        rooted=resolved_tree.rooted,
        strictly_bifurcating=True,
        missing_tip_taxa=sorted(input_tip_taxa - moved_tip_taxa),
        unexpected_tip_taxa=sorted(moved_tip_taxa - input_tip_taxa),
        moved_validation_errors=moved_tree.validation_errors(),
        affected_subtree_report=summarize_affected_subtrees(
            resolved_tree,
            moved_tree,
        ),
        affected_clade_ids=affected_clade_ids,
        pruned_edge_id=selected_candidate.pruned_clade_id,
        regraft_edge_id=selected_candidate.regraft_target_branch_id,
    )
    _validate_rooted_spr_move_application_report(report)
    return report


def prune_rooted_spr_subtree(
    tree: PhyloTree,
    prune_node_id: str,
) -> tuple[PhyloTree, TreeNode]:
    """Detach one non-root subtree and return the reduced tree plus detached clade."""
    working_tree = tree.copy().refresh()
    prune_node = working_tree.node_by_id(prune_node_id)
    parent = prune_node.parent
    if parent is None:
        raise AssertionError("rooted SPR prune node must not be the root")
    sibling = next(child for child in parent.children if child is not prune_node)
    grandparent = parent.parent
    detached_subtree = prune_node.copy()
    if grandparent is None:
        sibling.parent = None
        remainder_tree = PhyloTree(
            root=sibling,
            source_format=working_tree.source_format,
            rooted=True,
        )
        return remainder_tree.refresh(), detached_subtree
    grandparent.replace_children(
        [child if child is not parent else sibling for child in grandparent.children]
    )
    return working_tree.refresh(), detached_subtree


def _find_node_by_clade_id(tree: PhyloTree, clade_id: str) -> TreeNode:
    for node in tree.iter_nodes(order="preorder"):
        if rooted_spr_clade_id(node) == clade_id:
            return node
    raise KeyError(f"tree does not contain clade_id '{clade_id}'")


def _seed_missing_rooted_spr_branch_length(node: TreeNode) -> float:
    branch_length = node.branch_length
    if branch_length is not None:
        return branch_length
    return 0.0


def validate_rooted_spr_tree(tree: PhyloTree) -> None:
    """Require one structurally valid strictly bifurcating binary-root tree."""
    validation_errors = tree.validation_errors()
    if validation_errors:
        raise ValueError(
            "rooted SPR enumeration requires a structurally valid tree: "
            + "; ".join(validation_errors)
        )
    if len(tree.root.children) != 2:
        raise ValueError("rooted SPR enumeration requires a binary root")
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ValueError("rooted SPR enumeration requires a strictly bifurcating tree")


def enumerate_rooted_spr_neighbors(
    tree: PhyloTree | Path,
    *,
    budget: RootedSprEnumerationBudget | None = None,
) -> RootedSprNeighborhoodReport:
    """Enumerate rooted SPR neighbors exactly once by unique non-identity topology."""
    resolved_tree, input_tree_path = _resolve_rooted_spr_tree(tree)
    validate_rooted_spr_tree(resolved_tree)
    normalized_budget = validate_rooted_spr_enumeration_budget(budget)
    input_tip_taxa = set(resolved_tree.tip_names)
    input_topology_fingerprint = rooted_topology_fingerprint(resolved_tree)
    generated_move_candidate_count = 0
    identity_move_candidate_count = 0
    self_regraft_candidate_count = 0
    duplicate_move_topology_counts: dict[str, int] = {}
    neighbor_row_by_fingerprint: dict[str, RootedSprNeighborRow] = {}
    (
        move_candidates,
        skipped_pruned_clade_count,
        skipped_regraft_target_count,
    ) = _collect_rooted_spr_move_candidates(resolved_tree, budget=normalized_budget)
    for candidate in move_candidates:
        generated_move_candidate_count += 1
        if _candidate_is_self_regraft(candidate):
            self_regraft_candidate_count += 1
            continue
        neighbor_tree = apply_rooted_spr_move(resolved_tree, candidate)
        topology_fingerprint = rooted_topology_fingerprint(neighbor_tree)
        if topology_fingerprint == input_topology_fingerprint:
            identity_move_candidate_count += 1
            continue
        duplicate_move_topology_counts[topology_fingerprint] = (
            duplicate_move_topology_counts.get(topology_fingerprint, 0) + 1
        )
        if topology_fingerprint in neighbor_row_by_fingerprint:
            representative_row = neighbor_row_by_fingerprint[topology_fingerprint]
            neighbor_row_by_fingerprint[topology_fingerprint] = RootedSprNeighborRow(
                neighbor_index=representative_row.neighbor_index,
                representative_pruned_node_id=representative_row.representative_pruned_node_id,
                representative_pruned_clade_id=representative_row.representative_pruned_clade_id,
                representative_pruned_descendant_taxa=representative_row.representative_pruned_descendant_taxa,
                representative_regraft_target_branch_id=representative_row.representative_regraft_target_branch_id,
                representative_regraft_target_descendant_taxa=representative_row.representative_regraft_target_descendant_taxa,
                supporting_move_count=representative_row.supporting_move_count + 1,
                neighbor_tree_newick=representative_row.neighbor_tree_newick,
                neighbor_topology_fingerprint=representative_row.neighbor_topology_fingerprint,
                tip_order=representative_row.tip_order,
                validation_errors=representative_row.validation_errors,
            )
            continue
        neighbor_row_by_fingerprint[topology_fingerprint] = RootedSprNeighborRow(
            neighbor_index=len(neighbor_row_by_fingerprint) + 1,
            representative_pruned_node_id=candidate.pruned_node_id,
            representative_pruned_clade_id=candidate.pruned_clade_id,
            representative_pruned_descendant_taxa=list(candidate.pruned_descendant_taxa),
            representative_regraft_target_branch_id=candidate.regraft_target_branch_id,
            representative_regraft_target_descendant_taxa=(
                None
                if candidate.regraft_target_descendant_taxa is None
                else list(candidate.regraft_target_descendant_taxa)
            ),
            supporting_move_count=1,
            neighbor_tree_newick=neighbor_tree.to_newick(),
            neighbor_topology_fingerprint=topology_fingerprint,
            tip_order=neighbor_tree.tip_names,
            validation_errors=neighbor_tree.validation_errors(),
        )
    neighbor_rows = list(neighbor_row_by_fingerprint.values())
    report = RootedSprNeighborhoodReport(
        algorithm="rooted-spr-neighbor-enumeration",
        input_tree_path=input_tree_path,
        input_tree_newick=resolved_tree.to_newick(),
        tip_count=resolved_tree.tip_count,
        internal_node_count=resolved_tree.internal_node_count,
        rooted=resolved_tree.rooted,
        strictly_bifurcating=True,
        max_pruned_clade_count=normalized_budget.max_pruned_clade_count,
        max_regraft_target_count_per_pruned_clade=(
            normalized_budget.max_regraft_target_count_per_pruned_clade
        ),
        skipped_pruned_clade_count=skipped_pruned_clade_count,
        skipped_regraft_target_count=skipped_regraft_target_count,
        generated_move_candidate_count=generated_move_candidate_count,
        identity_move_candidate_count=identity_move_candidate_count,
        self_regraft_candidate_count=self_regraft_candidate_count,
        generated_neighbor_count=len(neighbor_rows),
        unique_neighbor_topology_count=len(neighbor_row_by_fingerprint),
        duplicate_move_neighbor_topologies=sorted(
            fingerprint
            for fingerprint, count in duplicate_move_topology_counts.items()
            if count > 1
        ),
        missing_tip_taxa=sorted(
            taxon
            for taxon in input_tip_taxa
            if any(taxon not in row.tip_order for row in neighbor_rows)
        ),
        unexpected_tip_taxa=sorted(
            {
                taxon
                for row in neighbor_rows
                for taxon in row.tip_order
                if taxon not in input_tip_taxa
            }
        ),
        input_validation_errors=[],
        neighbor_rows=neighbor_rows,
    )
    _validate_rooted_spr_neighbor_report(report)
    return report


def _resolve_rooted_spr_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def _candidate_is_self_regraft(candidate: RootedSprMoveCandidate) -> bool:
    if candidate.regraft_target_branch_id == candidate.pruned_clade_id:
        return True
    if candidate.regraft_target_descendant_taxa is None:
        return False
    return set(candidate.regraft_target_descendant_taxa).issubset(
        set(candidate.pruned_descendant_taxa)
    )


def _rooted_spr_affected_clade_ids(
    original_tree: PhyloTree,
    moved_tree: PhyloTree,
    candidate: RootedSprMoveCandidate,
) -> list[str]:
    changed_clade_ids = sorted(
        set(rooted_topology_signature_ids(original_tree))
        ^ set(rooted_topology_signature_ids(moved_tree))
    )
    affected = set(changed_clade_ids)
    affected.add(candidate.pruned_clade_id)
    if candidate.regraft_target_branch_id != _ROOT_REGRAFT_BRANCH_ID:
        affected.add(candidate.regraft_target_branch_id)
    return sorted(affected)


def _validate_rooted_spr_neighbor_report(
    report: RootedSprNeighborhoodReport,
) -> None:
    if report.self_regraft_candidate_count:
        raise ValueError("rooted SPR enumeration generated one or more self-regraft moves")
    if report.generated_neighbor_count != report.unique_neighbor_topology_count:
        raise ValueError(
            "rooted SPR enumeration did not collapse neighbors to unique topologies"
        )
    if report.missing_tip_taxa or report.unexpected_tip_taxa:
        raise ValueError("rooted SPR enumeration changed the input taxon set")
    invalid_neighbor_rows = [
        row.neighbor_index for row in report.neighbor_rows if row.validation_errors
    ]
    if invalid_neighbor_rows:
        raise ValueError("rooted SPR enumeration generated invalid neighbor trees")


def _validate_rooted_spr_move_application_report(
    report: RootedSprMoveApplicationReport,
) -> None:
    if not report.moved_topology_changed:
        raise ValueError("rooted SPR move application did not change the topology")
    if report.missing_tip_taxa or report.unexpected_tip_taxa:
        raise ValueError("rooted SPR move application changed the input taxon set")
    if report.moved_validation_errors:
        raise ValueError("rooted SPR move application generated an invalid tree")
    if not report.affected_clade_ids:
        raise ValueError(
            "rooted SPR move application did not report any affected clades"
        )


def write_rooted_spr_neighbor_table(
    path: Path,
    report: RootedSprNeighborhoodReport,
) -> Path:
    """Write one row per unique rooted SPR neighbor topology."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "neighbor_index",
                "representative_pruned_node_id",
                "representative_pruned_clade_id",
                "representative_pruned_descendant_taxa",
                "representative_regraft_target_branch_id",
                "representative_regraft_target_descendant_taxa",
                "supporting_move_count",
                "neighbor_topology_fingerprint",
                "tip_order",
                "validation_errors",
                "neighbor_tree_newick",
            ]
        )
    ]
    for row in report.neighbor_rows:
        lines.append(
            "\t".join(
                [
                    str(row.neighbor_index),
                    row.representative_pruned_node_id,
                    row.representative_pruned_clade_id,
                    ",".join(row.representative_pruned_descendant_taxa),
                    row.representative_regraft_target_branch_id,
                    ""
                    if row.representative_regraft_target_descendant_taxa is None
                    else ",".join(row.representative_regraft_target_descendant_taxa),
                    str(row.supporting_move_count),
                    row.neighbor_topology_fingerprint,
                    ",".join(row.tip_order),
                    ",".join(row.validation_errors),
                    row.neighbor_tree_newick,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_rooted_spr_run_json(
    path: Path,
    report: RootedSprNeighborhoodReport,
) -> Path:
    """Write one machine-readable rooted SPR neighborhood payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "input_tree_path": (
            None if report.input_tree_path is None else str(report.input_tree_path)
        ),
        "input_tree_newick": report.input_tree_newick,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "max_pruned_clade_count": report.max_pruned_clade_count,
        "max_regraft_target_count_per_pruned_clade": (
            report.max_regraft_target_count_per_pruned_clade
        ),
        "skipped_pruned_clade_count": report.skipped_pruned_clade_count,
        "skipped_regraft_target_count": report.skipped_regraft_target_count,
        "generated_move_candidate_count": report.generated_move_candidate_count,
        "identity_move_candidate_count": report.identity_move_candidate_count,
        "self_regraft_candidate_count": report.self_regraft_candidate_count,
        "generated_neighbor_count": report.generated_neighbor_count,
        "unique_neighbor_topology_count": report.unique_neighbor_topology_count,
        "duplicate_move_neighbor_topologies": report.duplicate_move_neighbor_topologies,
        "missing_tip_taxa": report.missing_tip_taxa,
        "unexpected_tip_taxa": report.unexpected_tip_taxa,
        "input_validation_errors": report.input_validation_errors,
        "neighbor_rows": [
            {
                "neighbor_index": row.neighbor_index,
                "representative_pruned_node_id": row.representative_pruned_node_id,
                "representative_pruned_clade_id": row.representative_pruned_clade_id,
                "representative_pruned_descendant_taxa": (
                    row.representative_pruned_descendant_taxa
                ),
                "representative_regraft_target_branch_id": (
                    row.representative_regraft_target_branch_id
                ),
                "representative_regraft_target_descendant_taxa": (
                    row.representative_regraft_target_descendant_taxa
                ),
                "supporting_move_count": row.supporting_move_count,
                "neighbor_tree_newick": row.neighbor_tree_newick,
                "neighbor_topology_fingerprint": row.neighbor_topology_fingerprint,
                "tip_order": row.tip_order,
                "validation_errors": row.validation_errors,
            }
            for row in report.neighbor_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_rooted_spr_artifacts(
    out_dir: Path,
    report: RootedSprNeighborhoodReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted SPR enumeration run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    neighbors_path = write_rooted_spr_neighbor_table(out_dir / "neighbors.tsv", report)
    run_json_path = write_rooted_spr_run_json(out_dir / "run.json", report)
    return {
        "input_tree_path": input_tree_path,
        "neighbors_path": neighbors_path,
        "run_json_path": run_json_path,
    }


def write_rooted_spr_move_run_json(
    path: Path,
    report: RootedSprMoveApplicationReport,
) -> Path:
    """Write one machine-readable rooted SPR move-application payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "input_tree_path": (
            None if report.input_tree_path is None else str(report.input_tree_path)
        ),
        "input_tree_newick": report.input_tree_newick,
        "input_topology_fingerprint": report.input_topology_fingerprint,
        "selected_move_index": report.selected_move_index,
        "available_move_count": report.available_move_count,
        "max_pruned_clade_count": report.max_pruned_clade_count,
        "max_regraft_target_count_per_pruned_clade": (
            report.max_regraft_target_count_per_pruned_clade
        ),
        "selected_pruned_node_id": report.selected_pruned_node_id,
        "selected_pruned_clade_id": report.selected_pruned_clade_id,
        "selected_pruned_descendant_taxa": report.selected_pruned_descendant_taxa,
        "selected_regraft_target_branch_id": report.selected_regraft_target_branch_id,
        "selected_regraft_target_descendant_taxa": (
            report.selected_regraft_target_descendant_taxa
        ),
        "moved_tree_newick": report.moved_tree_newick,
        "moved_topology_fingerprint": report.moved_topology_fingerprint,
        "moved_topology_changed": report.moved_topology_changed,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "strictly_bifurcating": report.strictly_bifurcating,
        "missing_tip_taxa": report.missing_tip_taxa,
        "unexpected_tip_taxa": report.unexpected_tip_taxa,
        "moved_validation_errors": report.moved_validation_errors,
        "affected_subtrees": {
            "original_branch_clade_ids": (
                report.affected_subtree_report.original_branch_clade_ids
            ),
            "moved_branch_clade_ids": (
                report.affected_subtree_report.moved_branch_clade_ids
            ),
            "retired_branch_clade_ids": (
                report.affected_subtree_report.retired_branch_clade_ids
            ),
            "introduced_branch_clade_ids": (
                report.affected_subtree_report.introduced_branch_clade_ids
            ),
            "affected_branch_clade_ids": (
                report.affected_subtree_report.affected_branch_clade_ids
            ),
            "unaffected_branch_clade_ids": (
                report.affected_subtree_report.unaffected_branch_clade_ids
            ),
        },
        "affected_clade_ids": report.affected_clade_ids,
        "pruned_edge_id": report.pruned_edge_id,
        "regraft_edge_id": report.regraft_edge_id,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rooted_spr_move_artifacts(
    out_dir: Path,
    report: RootedSprMoveApplicationReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted SPR move application."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    moved_tree_path = write_newick(
        out_dir / "moved_tree.nwk",
        loads_newick(report.moved_tree_newick),
    )
    run_json_path = write_rooted_spr_move_run_json(out_dir / "run.json", report)
    return {
        "input_tree_path": input_tree_path,
        "moved_tree_path": moved_tree_path,
        "run_json_path": run_json_path,
    }
