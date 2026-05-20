from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.taxa import inspect_tree_taxon_identity
from bijux_phylogenetics.trees import summarize_tree_shape_from_tree

from .branch_review import (
    _branch_length_health,
    _branch_length_status,
    _branch_length_summary,
    _branch_outliers,
    _internal_node_child_counts,
    _long_branch_taxa,
    _missing_internal_branch_nodes,
    _missing_terminal_branch_taxa,
    _singleton_internal_nodes,
    _tree_diameter,
    _tree_quality_score,
    _tree_quality_warnings,
)
from .label_review import (
    _internal_label_conflicts,
    _internal_label_diagnostics,
    _unsafe_external_labels,
)
from .models import TreeInspectionReport
from .structure import (
    _edge_count,
    _load_tree,
    _node_count,
    _polytomy_nodes,
    _root_state_confidence,
    _stable_node_identities,
    _ultrametric,
)


def inspect_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeInspectionReport:
    """Inspect a tree file and return lightweight summary metrics."""
    tree = _load_tree(path, source_format=source_format)
    shape = summarize_tree_shape_from_tree(tree, source_path=path)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    branch_lengths = [
        node.branch_length for node in tree.iter_nodes() if node is not tree.root
    ]
    polytomy_nodes = _polytomy_nodes(tree)
    branch_length_status = _branch_length_status(tree)
    internal_child_counts = _internal_node_child_counts(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    zero_length_branch_count = sum(1 for length in branch_lengths if length == 0)
    ultrametric = _ultrametric(tree)
    branch_length_summary = _branch_length_summary(tree)
    long_branch_taxa = _long_branch_taxa(tree)
    long_branch_outliers, short_branch_outliers = _branch_outliers(tree)
    (
        likely_support_labels,
        likely_named_internal_labels,
        suspicious_support_value_ranges,
        mixed_support_scales,
    ) = _internal_label_diagnostics(tree)
    internal_label_conflicts = _internal_label_conflicts(
        likely_support_labels,
        likely_named_internal_labels,
        suspicious_support_value_ranges,
        mixed_support_scales,
    )
    root_state_confidence = _root_state_confidence(tree)
    stable_node_identities = _stable_node_identities(tree)
    unsafe_external_labels = _unsafe_external_labels(tree)
    taxon_identity_audit = inspect_tree_taxon_identity(tree)
    _, _, negative_branch_count = _branch_length_health(tree)
    tree_quality_warnings = _tree_quality_warnings(
        tree,
        branch_length_status=branch_length_status,
        zero_length_branch_count=zero_length_branch_count,
        negative_branch_count=negative_branch_count,
        polytomy_nodes=polytomy_nodes,
        unusually_imbalanced=shape.unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        star_like=shape.star_like,
        comb_like=shape.comb_like,
    )
    warnings: list[str] = []
    if zero_length_branch_count:
        warnings.append("tree contains zero-length branches")
    if missing_internal_branch_nodes:
        warnings.append("tree contains internal branches without lengths")
    if missing_terminal_branch_taxa:
        warnings.append("tree contains terminal branches without lengths")
    if singleton_internal_nodes:
        warnings.append("tree contains singleton internal nodes")
    if polytomy_nodes:
        warnings.append("tree contains one or more polytomies")
    if branch_length_status == "partial":
        warnings.append("tree contains partial branch lengths")
    if branch_length_status == "absent":
        warnings.append("tree contains no branch lengths")
    if (
        taxon_identity_audit.whitespace_variants
        or taxon_identity_audit.underscore_space_collisions
    ):
        warnings.append(
            "tree contains taxon labels with whitespace or underscore identity collisions"
        )
    if (
        taxon_identity_audit.case_collisions
        or taxon_identity_audit.suspicious_near_duplicates
    ):
        warnings.append(
            "tree contains potentially ambiguous near-duplicate taxon labels"
        )
    return TreeInspectionReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        node_count=_node_count(tree),
        internal_node_count=tree.internal_node_count,
        edge_count=_edge_count(tree),
        clade_count=tree.internal_node_count,
        rooted=root_state_confidence.classification
        in {"explicitly_rooted", "apparently_rooted"},
        root_state_confidence=root_state_confidence,
        is_binary=all(
            node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes()
        ),
        internal_child_counts=internal_child_counts,
        singleton_internal_nodes=singleton_internal_nodes,
        polytomy_count=len(polytomy_nodes),
        polytomy_nodes=polytomy_nodes,
        has_branch_lengths=any(length is not None for length in branch_lengths),
        branch_length_status=branch_length_status,
        missing_internal_branch_nodes=missing_internal_branch_nodes,
        missing_terminal_branch_taxa=missing_terminal_branch_taxa,
        is_ultrametric=ultrametric,
        total_branch_length=tree.total_branch_length(),
        branch_length_summary=branch_length_summary,
        tree_diameter=_tree_diameter(tree),
        zero_length_branch_count=zero_length_branch_count,
        min_root_to_tip=min(lengths) if lengths else None,
        max_root_to_tip=max(lengths) if lengths else None,
        max_depth=shape.tree_height_edges,
        mean_depth=shape.mean_tip_depth_edges,
        colless_imbalance_index=shape.colless_imbalance_index,
        normalized_colless_imbalance=shape.normalized_colless_imbalance,
        sackin_imbalance_index=shape.sackin_imbalance_index,
        unusually_imbalanced=shape.unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        long_branch_outliers=long_branch_outliers,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        likely_support_labels=likely_support_labels,
        likely_named_internal_labels=likely_named_internal_labels,
        internal_label_conflicts=internal_label_conflicts,
        stable_node_identities=stable_node_identities,
        unsafe_external_labels=unsafe_external_labels,
        taxon_identity_audit=taxon_identity_audit,
        star_like=shape.star_like,
        comb_like=shape.comb_like,
        tree_quality_score=_tree_quality_score(tree_quality_warnings),
        tree_quality_warnings=tree_quality_warnings,
        imbalance_summary=shape.imbalance_summary,
        cherry_count=shape.cherry_count,
        taxa=sorted(tree.tip_names),
        warnings=warnings,
    )
