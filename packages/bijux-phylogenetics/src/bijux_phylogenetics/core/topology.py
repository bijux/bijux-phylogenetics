from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
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
    tip_order: list[str]
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


def _leaf_count(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_leaf_count(child) for child in node.children)


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
        raise ValueError("at least one outgroup taxon is required")

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
        raise ValueError(
            f"none of the requested outgroup taxa were found in {tree_path}"
        )

    biophylo_tree.root_with_outgroup(
        *[clade for clade in matched_clades if clade is not None]
    )
    rooted_tree = tree_from_biophylo(biophylo_tree, source_format=tree.source_format)
    summary = _summarize_transformation(
        tree, rooted_tree, transformation="root-outgroup"
    )
    return rooted_tree, TreeRootingReport(
        tree_path=tree_path,
        strategy="outgroup",
        requested_taxa=sorted(outgroup_taxa),
        matched_taxa=sorted(matched_taxa),
        absent_taxa=absent_taxa,
        tip_order=rooted_tree.tip_names,
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
    summary = _summarize_transformation(
        tree, rerooted_tree, transformation="reroot-midpoint"
    )
    return rerooted_tree, TreeRootingReport(
        tree_path=tree_path,
        strategy="midpoint",
        requested_taxa=[],
        matched_taxa=[],
        absent_taxa=[],
        tip_order=rerooted_tree.tip_names,
        summary=summary,
    )


def unroot_tree(tree_path: Path) -> tuple[PhyloTree, TreeRootingReport]:
    """Convert a rooted binary tree into an explicit unrooted trifurcation."""
    tree = load_tree(tree_path)
    if len(tree.root.children) != 2:
        raise ValueError(
            f"tree is not rooted under the repository convention: {tree_path}"
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
        TreeNode(
            name=child.name,
            branch_length=_combine_branch_lengths(
                child.branch_length, expanded.branch_length
            ),
            children=[_clone_node(grandchild) for grandchild in child.children],
        )
        for child in expanded.children
    ]
    new_children.append(_clone_node(retained))
    unrooted_tree = PhyloTree(
        root=TreeNode(name=tree.root.name, branch_length=None, children=new_children),
        source_format=tree.source_format,
        rooted=False,
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
        tip_order=unrooted_tree.tip_names,
        summary=summary,
    )
