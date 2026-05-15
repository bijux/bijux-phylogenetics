from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import TreeRootingError
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class CladeExtractionReport:
    """Explicit record of subtree extraction for a named clade."""

    tree_path: Path
    clade_name: str
    tip_count: int
    taxa: list[str]
    retained_all_requested_descendants: bool
    missing_requested_descendants: list[str]
    unexpected_retained_taxa: list[str]
    summary: TreeTransformationSummary


@dataclass(slots=True)
class BranchCollapseReport:
    """Explicit record of internal branches collapsed by a length threshold."""

    tree_path: Path
    threshold: float
    collapsed_clades: list[str]
    topology_preserved: bool
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeOrderingReport:
    """Explicit record of a deterministic child-ordering transform."""

    tree_path: Path
    strategy: str
    tip_order: list[str]
    rooted_topology_preserved: bool
    unrooted_topology_preserved: bool
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeRootingReport:
    """Explicit record of a tree rooting transform."""

    tree_path: Path
    strategy: str
    requested_taxa: list[str]
    matched_taxa: list[str]
    absent_taxa: list[str]
    ingroup_taxa: list[str]
    outgroup_monophyletic: bool | None
    outgroup_mrca_taxa: list[str]
    outgroup_mrca_extra_taxa: list[str]
    rooted_outgroup_taxa: list[str]
    rooted_ingroup_taxa: list[str]
    tip_order: list[str]
    warnings: list[str]
    midpoint_anchor_taxa: list[str]
    midpoint_path_length: float | None
    midpoint_distance_from_anchor: float | None
    midpoint_position_kind: str | None
    midpoint_anchor_side_taxa: list[str]
    midpoint_opposite_side_taxa: list[str]
    midpoint_suitable: bool | None
    summary: TreeTransformationSummary


@dataclass(slots=True)
class TreeTransformationSummary:
    """Before/after summary for a tree transformation."""

    transformation: str
    original_tip_count: int
    transformed_tip_count: int
    retained_taxa: list[str]
    removed_taxa: list[str]
    added_taxa: list[str]
    original_internal_node_count: int
    transformed_internal_node_count: int
    nodes_changed: list[str]
    original_total_branch_length: float
    transformed_total_branch_length: float
    branch_length_delta: float
    branch_lengths_affected: list[str]


@dataclass(slots=True)
class _TopologyComparison:
    topology_equal: bool
    same_unrooted_topology: bool


def _clone_node(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[_clone_node(child) for child in node.children],
    )


def _find_named_nodes(node: TreeNode, *, clade_name: str) -> list[TreeNode]:
    matches: list[TreeNode] = []
    if node.name == clade_name:
        matches.append(node)
    for child in node.children:
        matches.extend(_find_named_nodes(child, clade_name=clade_name))
    return matches


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _node_signature(node: TreeNode) -> str:
    taxa = _descendant_taxa(node)
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def _root_to_tip_paths(tree: PhyloTree) -> dict[str, list[TreeNode]]:
    paths: dict[str, list[TreeNode]] = {}

    def visit(node: TreeNode, path: list[TreeNode]) -> None:
        next_path = [*path, node]
        if node.is_leaf():
            if node.name is not None:
                paths[node.name] = next_path
            return
        for child in node.children:
            visit(child, next_path)

    visit(tree.root, [])
    return paths


def _common_prefix_length(left: list[TreeNode], right: list[TreeNode]) -> int:
    length = 0
    for left_node, right_node in zip(left, right, strict=False):
        if left_node is not right_node:
            break
        length += 1
    return length


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


def _join_taxa(taxa: list[str]) -> str:
    return ",".join(taxa)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(round(value, 15))


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

    if len(outgroup_children) == 1 and outgroup_children[0].is_leaf() and len(root_children) == 2:
        outgroup_child = outgroup_children[0]
        ingroup_child = ingroup_children[0]
        if (
            (outgroup_child.branch_length or 0.0) == 0.0
            and ingroup_child.branch_length is not None
        ):
            outgroup_child.branch_length = ingroup_child.branch_length
            ingroup_child.branch_length = 0.0
        return rooted_tree

    if len(outgroup_children) > 1:
        rooted_tree.root.children = [
            *ingroup_children,
            TreeNode(branch_length=0.0, children=outgroup_children),
        ]
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


def _leaf_count(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_leaf_count(child) for child in node.children)


def _is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes())


def _combine_branch_lengths(base: float | None, extra: float | None) -> float | None:
    if base is None and extra is None:
        return None
    if base is None:
        return extra
    if extra is None:
        return base
    return base + extra


def _branch_length_affecting_nodes(
    original: TreeNode, transformed: TreeNode
) -> list[str]:
    original_nodes = {
        _node_signature(node): node.branch_length
        for node in original.iter_nodes()
        if node is not original
    }
    transformed_nodes = {
        _node_signature(node): node.branch_length
        for node in transformed.iter_nodes()
        if node is not transformed
    }
    affected = {
        signature
        for signature in set(original_nodes) | set(transformed_nodes)
        if original_nodes.get(signature) != transformed_nodes.get(signature)
    }
    return sorted(affected)


def _changed_node_signatures(original: TreeNode, transformed: TreeNode) -> list[str]:
    original_nodes = {_node_signature(node) for node in original.iter_nodes()}
    transformed_nodes = {_node_signature(node) for node in transformed.iter_nodes()}
    return sorted(original_nodes.symmetric_difference(transformed_nodes))


def _summarize_transformation(
    original: PhyloTree,
    transformed: PhyloTree,
    *,
    transformation: str,
    extra_changed_nodes: list[str] | None = None,
) -> TreeTransformationSummary:
    original_taxa = sorted(original.tip_names)
    transformed_taxa = sorted(transformed.tip_names)
    original_taxa_set = set(original_taxa)
    transformed_taxa_set = set(transformed_taxa)
    nodes_changed = _changed_node_signatures(original.root, transformed.root)
    if extra_changed_nodes:
        nodes_changed = sorted(set(nodes_changed) | set(extra_changed_nodes))
    branch_lengths_affected = _branch_length_affecting_nodes(
        original.root, transformed.root
    )
    original_total = round(original.total_branch_length(), 15)
    transformed_total = round(transformed.total_branch_length(), 15)
    return TreeTransformationSummary(
        transformation=transformation,
        original_tip_count=original.tip_count,
        transformed_tip_count=transformed.tip_count,
        retained_taxa=sorted(original_taxa_set & transformed_taxa_set),
        removed_taxa=sorted(original_taxa_set - transformed_taxa_set),
        added_taxa=sorted(transformed_taxa_set - original_taxa_set),
        original_internal_node_count=original.internal_node_count,
        transformed_internal_node_count=transformed.internal_node_count,
        nodes_changed=nodes_changed,
        original_total_branch_length=original_total,
        transformed_total_branch_length=transformed_total,
        branch_length_delta=round(transformed_total - original_total, 15),
        branch_lengths_affected=branch_lengths_affected,
    )


def _informative_clades(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()
        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if 1 < len(taxa) < len(shared_taxa):
            clades.add(frozenset(taxa))
        return taxa

    visit(tree.root)
    return clades


def _canonical_bipartition(taxa: set[str], universe: set[str]) -> frozenset[str]:
    complement = universe - taxa
    left = sorted(taxa)
    right = sorted(complement)
    if (len(left), left) <= (len(right), right):
        return frozenset(taxa)
    return frozenset(complement)


def _unrooted_splits(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    splits: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()
        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if node is not tree.root and 1 < len(taxa) < len(shared_taxa) - 1:
            splits.add(_canonical_bipartition(taxa, shared_taxa))
        return taxa

    visit(tree.root)
    return splits


def _compare_tree_topology(
    original: PhyloTree, transformed: PhyloTree
) -> _TopologyComparison:
    shared_taxa = set(original.tip_names) & set(transformed.tip_names)
    if len(shared_taxa) < 2:
        return _TopologyComparison(
            topology_equal=original.tip_names == transformed.tip_names,
            same_unrooted_topology=True,
        )
    left_clades = _informative_clades(original, shared_taxa)
    right_clades = _informative_clades(transformed, shared_taxa)
    return _TopologyComparison(
        topology_equal=left_clades == right_clades,
        same_unrooted_topology=_unrooted_splits(original, shared_taxa)
        == _unrooted_splits(transformed, shared_taxa),
    )


def _collapse_short_internal_branches(
    node: TreeNode,
    *,
    threshold: float,
    collapsed_clades: list[str],
) -> TreeNode:
    if node.is_leaf():
        return TreeNode(name=node.name, branch_length=node.branch_length, children=[])

    rewritten_children: list[TreeNode] = []
    for child in node.children:
        rewritten_child = _collapse_short_internal_branches(
            child,
            threshold=threshold,
            collapsed_clades=collapsed_clades,
        )
        should_collapse = (
            not rewritten_child.is_leaf()
            and rewritten_child.branch_length is not None
            and rewritten_child.branch_length < threshold
        )
        if should_collapse:
            collapsed_clades.append(_node_signature(rewritten_child))
            rewritten_children.extend(
                TreeNode(
                    name=grandchild.name,
                    branch_length=grandchild.branch_length,
                    children=[
                        _clone_node(child_node) for child_node in grandchild.children
                    ],
                )
                for grandchild in rewritten_child.children
            )
            continue
        rewritten_children.append(rewritten_child)

    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=rewritten_children,
    )


def extract_named_clade(
    tree_path: Path,
    *,
    clade_name: str,
) -> tuple[PhyloTree, CladeExtractionReport]:
    """Extract a named clade as a standalone subtree."""
    tree = load_tree(tree_path)
    matches = _find_named_nodes(tree.root, clade_name=clade_name)
    if not matches:
        raise ValueError(f"clade '{clade_name}' was not found in {tree_path}")
    if len(matches) > 1:
        raise ValueError(f"clade '{clade_name}' is ambiguous in {tree_path}")

    subtree_root = _clone_node(matches[0])
    subtree_root.branch_length = None
    subtree = PhyloTree(
        root=subtree_root, source_format=tree.source_format, rooted=tree.rooted
    )
    expected_taxa = sorted(_descendant_taxa(matches[0]))
    observed_taxa = sorted(subtree.tip_names)
    summary = _summarize_transformation(
        tree, subtree, transformation="extract-named-clade"
    )
    return subtree, CladeExtractionReport(
        tree_path=tree_path,
        clade_name=clade_name,
        tip_count=subtree.tip_count,
        taxa=observed_taxa,
        retained_all_requested_descendants=observed_taxa == expected_taxa,
        missing_requested_descendants=sorted(set(expected_taxa) - set(observed_taxa)),
        unexpected_retained_taxa=sorted(set(observed_taxa) - set(expected_taxa)),
        summary=summary,
    )


def collapse_branches_below_length(
    tree_path: Path,
    *,
    threshold: float,
) -> tuple[PhyloTree, BranchCollapseReport]:
    """Collapse short internal branches into parent-level polytomies."""
    if threshold < 0:
        raise ValueError(f"threshold must be non-negative, got {threshold}")

    tree = load_tree(tree_path)
    collapsed_clades: list[str] = []
    collapsed_root = _collapse_short_internal_branches(
        tree.root,
        threshold=threshold,
        collapsed_clades=collapsed_clades,
    )
    collapsed_root.branch_length = None
    collapsed_tree = PhyloTree(
        root=collapsed_root, source_format=tree.source_format, rooted=tree.rooted
    )
    comparison = _compare_tree_topology(tree, collapsed_tree)
    summary = _summarize_transformation(
        tree,
        collapsed_tree,
        transformation="collapse-short-branches",
        extra_changed_nodes=collapsed_clades,
    )
    return collapsed_tree, BranchCollapseReport(
        tree_path=tree_path,
        threshold=threshold,
        collapsed_clades=sorted(collapsed_clades),
        topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def _order_tree(node: TreeNode, *, strategy: str) -> TreeNode:
    if node.is_leaf():
        return TreeNode(name=node.name, branch_length=node.branch_length, children=[])

    ordered_children = [
        _order_tree(child, strategy=strategy) for child in node.children
    ]
    if strategy == "ladderize":
        ordered_children.sort(
            key=lambda child: (-_leaf_count(child), _descendant_taxa(child))
        )
    elif strategy == "alphabetical":
        ordered_children.sort(key=_descendant_taxa)
    else:
        raise ValueError(f"unsupported ordering strategy: {strategy}")

    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=ordered_children,
    )


def _rotate_named_node(
    node: TreeNode,
    *,
    clade_name: str,
) -> tuple[TreeNode, int]:
    rotated_children = []
    match_count = 1 if node.name == clade_name else 0
    for child in node.children:
        rotated_child, child_matches = _rotate_named_node(child, clade_name=clade_name)
        rotated_children.append(rotated_child)
        match_count += child_matches
    if node.name == clade_name:
        rotated_children = list(reversed(rotated_children))
    return (
        TreeNode(
            name=node.name,
            branch_length=node.branch_length,
            children=rotated_children,
        ),
        match_count,
    )


def _rotate_all_nodes(node: TreeNode) -> TreeNode:
    rotated_children = [_rotate_all_nodes(child) for child in node.children]
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=list(reversed(rotated_children)),
    )


def ladderize_tree(tree_path: Path) -> tuple[PhyloTree, TreeOrderingReport]:
    """Ladderize a tree deterministically by descendant clade size."""
    tree = load_tree(tree_path)
    ladderized_tree = PhyloTree(
        root=_order_tree(tree.root, strategy="ladderize"),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, ladderized_tree)
    summary = _summarize_transformation(
        tree, ladderized_tree, transformation="ladderize-tree"
    )
    return ladderized_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="ladderize",
        tip_order=ladderized_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def rotate_named_node(
    tree_path: Path, *, clade_name: str
) -> tuple[PhyloTree, TreeOrderingReport]:
    """Reverse the child order at one named internal node."""
    tree = load_tree(tree_path)
    rotated_root, match_count = _rotate_named_node(tree.root, clade_name=clade_name)
    if match_count == 0:
        raise ValueError(f"clade '{clade_name}' was not found in {tree_path}")
    if match_count > 1:
        raise ValueError(f"clade '{clade_name}' is ambiguous in {tree_path}")
    rotated_tree = PhyloTree(
        root=rotated_root,
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, rotated_tree)
    summary = _summarize_transformation(
        tree,
        rotated_tree,
        transformation="rotate-named-node",
        extra_changed_nodes=[clade_name],
    )
    return rotated_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy=f"rotate:{clade_name}",
        tip_order=rotated_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def rotate_all_internal_nodes(tree_path: Path) -> tuple[PhyloTree, TreeOrderingReport]:
    """Reverse the child order at every internal node."""
    tree = load_tree(tree_path)
    rotated_tree = PhyloTree(
        root=_rotate_all_nodes(tree.root),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, rotated_tree)
    summary = _summarize_transformation(
        tree,
        rotated_tree,
        transformation="rotate-all-internal-nodes",
    )
    return rotated_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="rotate-all",
        tip_order=rotated_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def sort_tree_tips_alphabetically(
    tree_path: Path,
) -> tuple[PhyloTree, TreeOrderingReport]:
    """Sort tree children recursively by alphabetical descendant tip order."""
    tree = load_tree(tree_path)
    sorted_tree = PhyloTree(
        root=_order_tree(tree.root, strategy="alphabetical"),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, sorted_tree)
    summary = _summarize_transformation(
        tree, sorted_tree, transformation="sort-tree-tips"
    )
    return sorted_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="alphabetical",
        tip_order=sorted_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def root_tree_on_outgroup(
    tree_path: Path,
    *,
    outgroup_taxa: list[str],
) -> tuple[PhyloTree, TreeRootingReport]:
    """Root a tree on one or more named outgroup taxa."""
    if not outgroup_taxa:
        raise TreeRootingError("at least one outgroup taxon is required")

    tree = load_tree(tree_path)
    biophylo_tree = tree_to_biophylo(tree)
    matched_clades = [
        next(biophylo_tree.find_clades(name=taxon), None) for taxon in outgroup_taxa
    ]
    matched_taxa = [
        taxon
        for taxon, clade in zip(outgroup_taxa, matched_clades, strict=True)
        if clade is not None
    ]
    absent_taxa = sorted(
        taxon
        for taxon, clade in zip(outgroup_taxa, matched_clades, strict=True)
        if clade is None
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
    if matched_taxa:
        outgroup_mrca = biophylo_tree.common_ancestor(
            *[clade for clade in matched_clades if clade is not None]
        )
        outgroup_mrca_taxa = _biophylo_clade_taxa(outgroup_mrca)
        outgroup_mrca_extra_taxa = sorted(set(outgroup_mrca_taxa) - requested_taxa_set)
        outgroup_monophyletic = outgroup_mrca_extra_taxa == []
        if not outgroup_monophyletic:
            raise TreeRootingError(
                "requested outgroup taxa are not monophyletic in the input tree",
                code="outgroup_not_monophyletic",
                details={
                    "matched_taxa": sorted(matched_taxa),
                    "outgroup_mrca_taxa": outgroup_mrca_taxa,
                    "outgroup_mrca_extra_taxa": outgroup_mrca_extra_taxa,
                },
            )
        if absent_taxa:
            warnings.append(
                "one or more requested outgroup taxa were absent from the input tree"
            )

    biophylo_tree.root_with_outgroup(
        biophylo_tree.common_ancestor(
            *[clade for clade in matched_clades if clade is not None]
        )
    )
    rooted_tree = tree_from_biophylo(biophylo_tree, source_format=tree.source_format)
    rooted_tree = _normalize_outgroup_rooting_to_ape(
        rooted_tree,
        requested_taxa_set=requested_taxa_set,
    )
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
    new_children = [
        _clone_node(child) for child in expanded.children
    ]
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
