from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import TreeRootingError

from .models import TreeRootingReport
from .subtree import _common_prefix_length, _mrca_node_from_taxa, _root_to_tip_paths
from .transformation import (
    _clone_node,
    _combine_branch_lengths,
    _descendant_taxa,
    _format_optional_float,
    _join_taxa,
    _summarize_transformation,
)


def _adjacent_nodes(node: TreeNode) -> list[TreeNode]:
    neighbors = list(node.children)
    if node.parent is not None:
        neighbors.append(node.parent)
    return neighbors


def _edge_length_between(left: TreeNode, right: TreeNode) -> float | None:
    if right.parent is left:
        return right.branch_length
    if left.parent is right:
        return left.branch_length
    raise ValueError("requested nodes do not share one edge")


def _copy_node_payload(node: TreeNode, *, branch_length: float | None) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=branch_length,
        children=[],
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )


def _clone_subtree_away_from(
    node: TreeNode, *, from_neighbor: TreeNode | None
) -> TreeNode:
    clone = _copy_node_payload(node, branch_length=node.branch_length)
    for neighbor in _adjacent_nodes(node):
        if neighbor is from_neighbor:
            continue
        child = _clone_subtree_component(
            neighbor,
            from_neighbor=node,
            incoming_length=_edge_length_between(node, neighbor),
        )
        clone.append_child(child)
    return clone


def _clone_subtree_component(
    node: TreeNode,
    *,
    from_neighbor: TreeNode,
    incoming_length: float | None,
) -> TreeNode:
    neighbors = [
        neighbor for neighbor in _adjacent_nodes(node) if neighbor is not from_neighbor
    ]
    if neighbors and len(neighbors) == 1 and not node.is_leaf():
        next_neighbor = neighbors[0]
        return _clone_subtree_component(
            next_neighbor,
            from_neighbor=node,
            incoming_length=_combine_branch_lengths(
                incoming_length,
                _edge_length_between(node, next_neighbor),
            ),
        )
    if node.parent is not None and node.parent in neighbors:
        neighbors = [
            node.parent,
            *[neighbor for neighbor in neighbors if neighbor is not node.parent],
        ]
    clone = _copy_node_payload(node, branch_length=incoming_length)
    for neighbor in neighbors:
        child = _clone_subtree_component(
            neighbor,
            from_neighbor=node,
            incoming_length=_edge_length_between(node, neighbor),
        )
        clone.append_child(child)
    return clone


def _find_monophyletic_outgroup_node(
    tree: PhyloTree,
    *,
    requested_taxa: list[str],
) -> tuple[TreeNode, list[str], list[str]]:
    matched_node = _mrca_node_from_taxa(tree, requested_taxa)
    matched_taxa = _descendant_taxa(matched_node)
    extra_taxa = sorted(set(matched_taxa) - set(requested_taxa))
    if extra_taxa:
        raise TreeRootingError(
            "requested outgroup taxa are not monophyletic in the input tree",
            code="outgroup_not_monophyletic",
            details={
                "matched_taxa": sorted(requested_taxa),
                "outgroup_mrca_taxa": matched_taxa,
                "outgroup_mrca_extra_taxa": extra_taxa,
            },
        )
    return matched_node, matched_taxa, extra_taxa


def _root_tree_by_outgroup_node(
    tree: PhyloTree,
    *,
    outgroup_node: TreeNode,
) -> PhyloTree:
    if outgroup_node is tree.root:
        rooted_tree = tree.copy()
        rooted_tree.rooted = True
        if rooted_tree.root.branch_length is not None:
            rooted_tree.root.branch_length = None
        return rooted_tree

    parent = outgroup_node.parent
    if parent is None:
        raise TreeRootingError("requested outgroup cannot be resolved for rooting")

    if outgroup_node.is_leaf():
        outgroup_branch_length = outgroup_node.branch_length
        ingroup_initial_length = 0.0
    else:
        outgroup_branch_length = 0.0
        ingroup_initial_length = outgroup_node.branch_length

    rooted_tree = PhyloTree(
        root=TreeNode(
            children=[
                _clone_subtree_component(
                    parent,
                    from_neighbor=outgroup_node,
                    incoming_length=ingroup_initial_length,
                ),
                _copy_node_payload(outgroup_node, branch_length=outgroup_branch_length),
            ]
        ),
        source_format=tree.source_format,
        rooted=True,
    )
    rooted_tree.root.children[1].replace_children(
        _clone_subtree_away_from(outgroup_node, from_neighbor=parent).children
    )
    return _normalize_outgroup_rooting_to_ape(
        rooted_tree,
        requested_taxa_set=set(_descendant_taxa(outgroup_node)),
    )


def _analyze_midpoint_path(tree: PhyloTree) -> tuple[list[str], float, float, str]:
    paths = _root_to_tip_paths(tree)
    tip_names = sorted(paths)
    if len(tip_names) < 2:
        raise ValueError("midpoint rerooting requires at least two tips")

    tolerance = 1e-12
    best_anchor_taxa: list[str] = []
    best_path_length = -1.0
    best_distance_from_anchor = 0.0
    best_position_kind = "branch"

    for index, left_taxon in enumerate(tip_names[:-1]):
        for right_taxon in tip_names[index + 1 :]:
            left_path = paths[left_taxon]
            right_path = paths[right_taxon]
            prefix_length = _common_prefix_length(left_path, right_path)
            left_suffix = left_path[prefix_length:]
            right_suffix = right_path[prefix_length:]
            path_nodes = [*reversed(left_suffix), *right_suffix]
            path_length = sum((node.branch_length or 0.0) for node in path_nodes)
            if path_length < best_path_length - tolerance:
                continue

            midpoint_distance = path_length / 2.0
            traversed = 0.0
            position_kind = "branch"
            for node in path_nodes:
                traversed += node.branch_length or 0.0
                if abs(traversed - midpoint_distance) <= tolerance:
                    position_kind = "node"
                    break
                if traversed > midpoint_distance + tolerance:
                    position_kind = "branch"
                    break

            candidate_anchor_taxa = [left_taxon, right_taxon]
            candidate = (
                round(path_length, 15),
                candidate_anchor_taxa,
                position_kind,
            )
            best = (
                round(best_path_length, 15),
                best_anchor_taxa,
                best_position_kind,
            )
            if path_length > best_path_length + tolerance or candidate < best:
                best_anchor_taxa = candidate_anchor_taxa
                best_path_length = path_length
                best_distance_from_anchor = midpoint_distance
                best_position_kind = position_kind

    return (
        best_anchor_taxa,
        round(best_path_length, 15),
        round(best_distance_from_anchor, 15),
        best_position_kind,
    )


def _biophylo_clade_taxa(clade: object) -> list[str]:
    get_terminals = clade.get_terminals
    return sorted(
        terminal.name
        for terminal in get_terminals()
        if getattr(terminal, "name", None) is not None
    )


def _normalize_outgroup_rooting_to_ape(
    rooted_tree: PhyloTree,
    *,
    requested_taxa_set: set[str],
) -> PhyloTree:
    root_children = list(rooted_tree.root.children)
    outgroup_children = [
        child
        for child in root_children
        if set(_descendant_taxa(child)).issubset(requested_taxa_set)
    ]
    ingroup_children = [
        child for child in root_children if child not in outgroup_children
    ]
    if not outgroup_children:
        raise TreeRootingError(
            "rooted tree does not isolate the requested outgroup after rerooting",
            code="outgroup_root_not_isolated",
            details={"requested_taxa": sorted(requested_taxa_set)},
        )

    if (
        len(outgroup_children) == 1
        and outgroup_children[0].is_leaf()
        and len(root_children) == 2
    ):
        outgroup_child = outgroup_children[0]
        ingroup_child = ingroup_children[0]
        if (
            outgroup_child.branch_length or 0.0
        ) == 0.0 and ingroup_child.branch_length is not None:
            outgroup_child.branch_length = ingroup_child.branch_length
            ingroup_child.branch_length = 0.0
        return rooted_tree

    if len(outgroup_children) > 1:
        rooted_tree.root.replace_children(
            [
                *ingroup_children,
                TreeNode(branch_length=0.0, children=outgroup_children),
            ]
        )
    return rooted_tree


def write_tree_rooting_report(path: Path, report: TreeRootingReport) -> Path:
    """Write one durable TSV summary for a rooting transform."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "tree_path",
                "strategy",
                "requested_taxa",
                "matched_taxa",
                "absent_taxa",
                "ingroup_taxa",
                "outgroup_monophyletic",
                "outgroup_mrca_taxa",
                "outgroup_mrca_extra_taxa",
                "rooted_outgroup_taxa",
                "rooted_ingroup_taxa",
                "tip_order",
                "warnings",
                "midpoint_anchor_taxa",
                "midpoint_path_length",
                "midpoint_distance_from_anchor",
                "midpoint_position_kind",
                "midpoint_anchor_side_taxa",
                "midpoint_opposite_side_taxa",
                "midpoint_suitable",
            ]
        ),
        "\t".join(
            [
                str(report.tree_path),
                report.strategy,
                _join_taxa(report.requested_taxa),
                _join_taxa(report.matched_taxa),
                _join_taxa(report.absent_taxa),
                _join_taxa(report.ingroup_taxa),
                ""
                if report.outgroup_monophyletic is None
                else ("true" if report.outgroup_monophyletic else "false"),
                _join_taxa(report.outgroup_mrca_taxa),
                _join_taxa(report.outgroup_mrca_extra_taxa),
                _join_taxa(report.rooted_outgroup_taxa),
                _join_taxa(report.rooted_ingroup_taxa),
                _join_taxa(report.tip_order),
                " | ".join(report.warnings),
                _join_taxa(report.midpoint_anchor_taxa),
                _format_optional_float(report.midpoint_path_length),
                _format_optional_float(report.midpoint_distance_from_anchor),
                report.midpoint_position_kind or "",
                _join_taxa(report.midpoint_anchor_side_taxa),
                _join_taxa(report.midpoint_opposite_side_taxa),
                ""
                if report.midpoint_suitable is None
                else ("true" if report.midpoint_suitable else "false"),
            ]
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes())


def root_tree_on_outgroup(
    tree_path: Path,
    *,
    outgroup_taxa: list[str],
) -> tuple[PhyloTree, TreeRootingReport]:
    """Root a tree on one or more named outgroup taxa."""
    if not outgroup_taxa:
        raise TreeRootingError("at least one outgroup taxon is required")

    tree = load_tree(tree_path)
    matched_taxa = sorted(
        taxon for taxon in outgroup_taxa if taxon in set(tree.tip_names)
    )
    absent_taxa = sorted(
        taxon for taxon in outgroup_taxa if taxon not in set(tree.tip_names)
    )
    if not matched_taxa:
        raise TreeRootingError(
            f"none of the requested outgroup taxa were found in {tree_path}",
            code="outgroup_taxa_missing",
            details={
                "tree_path": str(tree_path),
                "requested_taxa": sorted(outgroup_taxa),
            },
        )

    requested_taxa_set = set(matched_taxa)
    ingroup_taxa = sorted(set(tree.tip_names) - requested_taxa_set)
    outgroup_mrca_taxa: list[str] = []
    outgroup_mrca_extra_taxa: list[str] = []
    outgroup_monophyletic: bool | None = None
    warnings: list[str] = []
    outgroup_node, outgroup_mrca_taxa, outgroup_mrca_extra_taxa = (
        _find_monophyletic_outgroup_node(
            tree,
            requested_taxa=matched_taxa,
        )
    )
    outgroup_monophyletic = True
    if absent_taxa:
        warnings.append(
            "one or more requested outgroup taxa were absent from the input tree"
        )

    rooted_tree = _root_tree_by_outgroup_node(tree, outgroup_node=outgroup_node)
    rooted_outgroup_taxa_set: set[str] = set()
    for child in rooted_tree.root.children:
        child_taxa = _descendant_taxa(child)
        child_taxa_set = set(child_taxa)
        if child_taxa_set and child_taxa_set.issubset(requested_taxa_set):
            rooted_outgroup_taxa_set.update(child_taxa_set)
    rooted_outgroup_taxa = sorted(rooted_outgroup_taxa_set)
    rooted_ingroup_taxa = sorted(set(rooted_tree.tip_names) - rooted_outgroup_taxa_set)
    if matched_taxa and not rooted_outgroup_taxa:
        raise TreeRootingError(
            "rooted tree does not isolate every matched outgroup taxon on one root child",
            code="outgroup_root_not_isolated",
            details={"matched_taxa": sorted(matched_taxa)},
        )
    summary = _summarize_transformation(
        tree, rooted_tree, transformation="root-outgroup"
    )
    return rooted_tree, TreeRootingReport(
        tree_path=tree_path,
        strategy="outgroup",
        requested_taxa=sorted(outgroup_taxa),
        matched_taxa=sorted(matched_taxa),
        absent_taxa=absent_taxa,
        ingroup_taxa=ingroup_taxa,
        outgroup_monophyletic=outgroup_monophyletic,
        outgroup_mrca_taxa=outgroup_mrca_taxa,
        outgroup_mrca_extra_taxa=outgroup_mrca_extra_taxa,
        rooted_outgroup_taxa=rooted_outgroup_taxa,
        rooted_ingroup_taxa=rooted_ingroup_taxa,
        tip_order=rooted_tree.tip_names,
        warnings=warnings,
        midpoint_anchor_taxa=[],
        midpoint_path_length=None,
        midpoint_distance_from_anchor=None,
        midpoint_position_kind=None,
        midpoint_anchor_side_taxa=[],
        midpoint_opposite_side_taxa=[],
        midpoint_suitable=None,
        summary=summary,
    )


def reroot_tree_by_midpoint(tree_path: Path) -> tuple[PhyloTree, TreeRootingReport]:
    """Reroot a tree by the midpoint of its longest tip-to-tip path."""
    tree = load_tree(tree_path)
    branch_lengths = tree.branch_lengths()
    if not branch_lengths or any(length is None for length in branch_lengths):
        raise ValueError(
            f"midpoint rerooting requires complete branch lengths in {tree_path}"
        )

    biophylo_tree = tree_to_biophylo(tree)
    biophylo_tree.root_at_midpoint()
    rerooted_tree = tree_from_biophylo(biophylo_tree, source_format=tree.source_format)
    (
        midpoint_anchor_taxa,
        midpoint_path_length,
        midpoint_distance_from_anchor,
        midpoint_position_kind,
    ) = _analyze_midpoint_path(tree)
    midpoint_root_partitions = sorted(
        _descendant_taxa(child) for child in rerooted_tree.root.children
    )
    midpoint_anchor_side_taxa: list[str] = []
    midpoint_opposite_side_taxa: list[str] = []
    warnings: list[str] = []
    midpoint_suitable = True
    if midpoint_anchor_taxa and len(midpoint_root_partitions) >= 2:
        anchor_taxon = midpoint_anchor_taxa[0]
        for partition in midpoint_root_partitions:
            if anchor_taxon in partition:
                midpoint_anchor_side_taxa = partition
            else:
                midpoint_opposite_side_taxa.extend(partition)
        midpoint_opposite_side_taxa = sorted(midpoint_opposite_side_taxa)
    if not _is_strictly_bifurcating(tree):
        warnings.append(
            "midpoint rooting is exploratory because the input tree is not strictly bifurcating"
        )
        midpoint_suitable = False
    summary = _summarize_transformation(
        tree, rerooted_tree, transformation="reroot-midpoint"
    )
    return rerooted_tree, TreeRootingReport(
        tree_path=tree_path,
        strategy="midpoint",
        requested_taxa=[],
        matched_taxa=[],
        absent_taxa=[],
        ingroup_taxa=sorted(rerooted_tree.tip_names),
        outgroup_monophyletic=None,
        outgroup_mrca_taxa=[],
        outgroup_mrca_extra_taxa=[],
        rooted_outgroup_taxa=[],
        rooted_ingroup_taxa=sorted(rerooted_tree.tip_names),
        tip_order=rerooted_tree.tip_names,
        warnings=warnings,
        midpoint_anchor_taxa=midpoint_anchor_taxa,
        midpoint_path_length=midpoint_path_length,
        midpoint_distance_from_anchor=midpoint_distance_from_anchor,
        midpoint_position_kind=midpoint_position_kind,
        midpoint_anchor_side_taxa=midpoint_anchor_side_taxa,
        midpoint_opposite_side_taxa=midpoint_opposite_side_taxa,
        midpoint_suitable=midpoint_suitable,
        summary=summary,
    )


def unroot_tree(tree_path: Path) -> tuple[PhyloTree, TreeRootingReport]:
    """Convert a rooted binary tree into an explicit unrooted trifurcation."""
    tree = load_tree(tree_path)
    root_child_count = len(tree.root.children)
    if root_child_count < 2:
        raise ValueError(f"tree is invalid and cannot be unrooted: {tree_path}")
    if root_child_count != 2:
        unchanged_tree = PhyloTree(
            root=_clone_node(tree.root),
            source_format=tree.source_format,
            rooted=False,
        )
        summary = _summarize_transformation(
            tree, unchanged_tree, transformation="unroot-tree"
        )
        return unchanged_tree, TreeRootingReport(
            tree_path=tree_path,
            strategy="unroot",
            requested_taxa=[],
            matched_taxa=[],
            absent_taxa=[],
            ingroup_taxa=sorted(unchanged_tree.tip_names),
            outgroup_monophyletic=None,
            outgroup_mrca_taxa=[],
            outgroup_mrca_extra_taxa=[],
            rooted_outgroup_taxa=[],
            rooted_ingroup_taxa=sorted(unchanged_tree.tip_names),
            tip_order=unchanged_tree.tip_names,
            warnings=[
                "input tree already behaves as an unrooted representation; returned unchanged"
            ],
            midpoint_anchor_taxa=[],
            midpoint_path_length=None,
            midpoint_distance_from_anchor=None,
            midpoint_position_kind=None,
            midpoint_anchor_side_taxa=[],
            midpoint_opposite_side_taxa=[],
            midpoint_suitable=None,
            summary=summary,
        )

    left, right = tree.root.children
    expandable_children = [child for child in (left, right) if not child.is_leaf()]
    if not expandable_children:
        raise ValueError(
            f"cannot unroot a two-tip tree under the repository convention: {tree_path}"
        )

    expanded = sorted(expandable_children, key=_descendant_taxa)[0]
    retained = right if expanded is left else left
    new_children = [_clone_node(child) for child in expanded.children]
    retained_child = _clone_node(retained)
    retained_child.branch_length = _combine_branch_lengths(
        retained.branch_length,
        expanded.branch_length,
    )
    new_children.append(retained_child)
    unrooted_tree = PhyloTree(
        root=TreeNode(name=tree.root.name, branch_length=None, children=new_children),
        source_format=tree.source_format,
        rooted=False,
    )
    warnings: list[str] = []
    if expanded.branch_length not in {None, 0.0}:
        warnings.append(
            "unrooting merged the removed root-edge length into the retained sibling branch to match ape::unroot"
        )
    summary = _summarize_transformation(
        tree, unrooted_tree, transformation="unroot-tree"
    )
    return unrooted_tree, TreeRootingReport(
        tree_path=tree_path,
        strategy="unroot",
        requested_taxa=[],
        matched_taxa=[],
        absent_taxa=[],
        ingroup_taxa=sorted(unrooted_tree.tip_names),
        outgroup_monophyletic=None,
        outgroup_mrca_taxa=[],
        outgroup_mrca_extra_taxa=[],
        rooted_outgroup_taxa=[],
        rooted_ingroup_taxa=sorted(unrooted_tree.tip_names),
        tip_order=unrooted_tree.tip_names,
        warnings=warnings,
        midpoint_anchor_taxa=[],
        midpoint_path_length=None,
        midpoint_distance_from_anchor=None,
        midpoint_position_kind=None,
        midpoint_anchor_side_taxa=[],
        midpoint_opposite_side_taxa=[],
        midpoint_suitable=None,
        summary=summary,
    )
