from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import stable_covariance
from bijux_phylogenetics.datasets.study_inputs import (
    align_tree_and_trait_table,
    load_taxon_table,
    validate_traits_table,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class NumericTraitSummary:
    """Numeric trait values retained after explicit phylogenetic pruning."""

    trait: str
    taxon_count: int
    taxa: list[str]
    mean: float
    variance: float
    minimum: float
    maximum: float


@dataclass(slots=True)
class ComparativeReadinessReport:
    """Readiness summary for a rooted tree plus one numeric trait."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    rooted: bool
    binary: bool
    complete_branch_lengths: bool
    negative_branch_lengths: bool
    minimum_branch_length: float | None
    tree_taxa: int
    analysis_taxa: list[str]
    missing_from_traits: list[str]
    extra_trait_taxa: list[str]
    pruned_missing_value_taxa: list[str]
    pruned_non_numeric_taxa: list[str]
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativeDataset:
    """Tree-plus-trait dataset normalized for comparative methods."""

    tree_path: Path
    traits_path: Path
    tree: PhyloTree
    taxon_column: str
    trait: str
    taxa: list[str]
    trait_values: list[float]
    covariance_matrix: list[list[float]]
    readiness: ComparativeReadinessReport


def summarize_numeric_trait_readiness(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> ComparativeReadinessReport:
    """Report whether a rooted tree and one numeric trait are ready for comparative analysis."""
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise ComparativeMethodError(f"trait table does not contain column '{trait}'")
    alignment = align_tree_and_trait_table(
        tree_path,
        traits_path,
        taxon_column=taxon_column,
        required_trait_columns=(trait,),
    )
    trait_report = validate_traits_table(traits_path, taxon_column=taxon_column)
    trait_column = next(
        (column for column in trait_report.trait_columns if column.name == trait), None
    )
    if trait_column is None:
        raise ComparativeMethodError(
            f"trait column '{trait}' is not available for validation"
        )

    overlapping_rows = {row[table.taxon_column]: row for row in alignment.rows}
    analysis_taxa: list[str] = []
    pruned_missing_value_taxa: list[str] = []
    pruned_non_numeric_taxa: list[str] = []
    for taxon in alignment.report.aligned_taxa:
        row = overlapping_rows.get(taxon)
        if row is None:
            continue
        raw = row[trait]
        if not raw:
            pruned_missing_value_taxa.append(taxon)
            continue
        try:
            float(raw)
        except ValueError:
            pruned_non_numeric_taxa.append(taxon)
            continue
        analysis_taxa.append(taxon)
    missing_from_traits = list(alignment.report.dropped_tree_taxa)
    extra_trait_taxa = list(alignment.report.dropped_trait_taxa)

    blockers: list[str] = []
    warnings: list[str] = []
    rooted = len(tree.root.children) == 2
    binary = all(
        node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes()
    )
    complete_branch_lengths = _has_complete_branch_lengths(tree)
    minimum_branch_length = _minimum_branch_length(tree)
    negative_branch_lengths = (
        minimum_branch_length is not None and minimum_branch_length < 0.0
    )
    if trait_column.kind != "numeric":
        blockers.append(
            f"trait column '{trait}' must be numeric for comparative analysis"
        )
    if not rooted:
        blockers.append("tree must be rooted for comparative analysis")
    if not complete_branch_lengths:
        blockers.append(
            "tree requires complete branch lengths for comparative analysis"
        )
    if negative_branch_lengths:
        blockers.append(
            "tree contains negative branch lengths that invalidate comparative analysis"
        )
    if len(analysis_taxa) < 3:
        blockers.append(
            "comparative analysis requires at least three taxa with numeric trait values"
        )
    if missing_from_traits:
        warnings.append(
            "trait table is missing one or more tree taxa and those taxa will be pruned"
        )
    if extra_trait_taxa:
        warnings.append("trait table contains taxa absent from the tree")
    if pruned_missing_value_taxa:
        warnings.append(
            "one or more overlapping taxa have missing trait values and will be pruned"
        )
    if pruned_non_numeric_taxa:
        warnings.append(
            "one or more overlapping taxa have non-numeric trait values and will be pruned"
        )

    return ComparativeReadinessReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        rooted=rooted,
        binary=binary,
        complete_branch_lengths=complete_branch_lengths,
        negative_branch_lengths=negative_branch_lengths,
        minimum_branch_length=minimum_branch_length,
        tree_taxa=alignment.report.original_tree_taxa,
        analysis_taxa=analysis_taxa,
        missing_from_traits=missing_from_traits,
        extra_trait_taxa=extra_trait_taxa,
        pruned_missing_value_taxa=pruned_missing_value_taxa,
        pruned_non_numeric_taxa=pruned_non_numeric_taxa,
        ready=not blockers,
        blockers=blockers,
        warnings=warnings,
    )


def summarize_numeric_trait(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> NumericTraitSummary:
    """Compute mean and sample variance for a numeric trait after phylogenetic pruning."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if len(readiness.analysis_taxa) < 2:
        raise ComparativeMethodError(
            "at least two taxa are required to summarize a numeric trait"
        )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    values_by_taxon = {
        row[table.taxon_column]: float(row[trait])
        for row in table.rows
        if row[table.taxon_column] in readiness.analysis_taxa and row[trait]
    }
    values = [values_by_taxon[taxon] for taxon in readiness.analysis_taxa]
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
    return NumericTraitSummary(
        trait=trait,
        taxon_count=len(readiness.analysis_taxa),
        taxa=readiness.analysis_taxa,
        mean=mean_value,
        variance=variance,
        minimum=min(values),
        maximum=max(values),
    )


def load_comparative_dataset(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    minimum_taxa: int = 3,
    require_rooted: bool = True,
    require_binary: bool = False,
) -> ComparativeDataset:
    """Load a comparative dataset over the tree/trait taxon intersection."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if require_rooted and not readiness.rooted:
        raise ComparativeMethodError(
            "tree must be rooted for this comparative method",
            details={
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "failure_reason": "comparative_tree_rooting_invalid",
                "scientific_explanation": (
                    "This comparative method needs a rooted phylogeny so ancestor-descendant covariance is defined."
                ),
                "likely_causes": [
                    "the supplied tree is unrooted or does not expose one biological root",
                ],
                "actionable_fixes": [
                    "root the tree on an explicit outgroup or midpoint policy before rerunning the comparative method",
                ],
                "evidence": {
                    "analysis_taxa": readiness.analysis_taxa,
                    "missing_from_traits": readiness.missing_from_traits,
                    "extra_trait_taxa": readiness.extra_trait_taxa,
                },
            },
        )
    if not readiness.complete_branch_lengths:
        raise ComparativeMethodError(
            "tree must contain complete branch lengths for this comparative method",
            details={
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "failure_reason": "comparative_branch_lengths_incomplete",
                "scientific_explanation": (
                    "Comparative covariance cannot be computed because one or more branches are missing lengths."
                ),
                "likely_causes": [
                    "the tree was exported without branch lengths",
                    "one or more internal or terminal branches have blank lengths",
                ],
                "actionable_fixes": [
                    "rerun the upstream inference or tree preparation step with branch lengths preserved",
                    "inspect the tree file for branches that lack numeric lengths",
                ],
                "evidence": {
                    "analysis_taxa": readiness.analysis_taxa,
                    "missing_from_traits": readiness.missing_from_traits,
                    "extra_trait_taxa": readiness.extra_trait_taxa,
                },
            },
        )
    if readiness.negative_branch_lengths:
        raise ComparativeMethodError(
            "tree contains negative branch lengths that invalidate this comparative method",
            details={
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "failure_reason": "comparative_negative_branch_lengths",
                "scientific_explanation": (
                    "This comparative method treats branch length as evolutionary variance or shared path length, so negative values make the analysis scientifically invalid."
                ),
                "likely_causes": [
                    "the tree file contains one or more negative branch lengths",
                    "branch scaling or export introduced a negative value on a non-root edge",
                ],
                "actionable_fixes": [
                    "repair or re-estimate the tree so every non-root branch length is non-negative",
                    "inspect the tree for scaling or export errors before rerunning the comparative method",
                ],
                "evidence": {
                    "analysis_taxa": readiness.analysis_taxa,
                    "minimum_branch_length": readiness.minimum_branch_length,
                },
            },
        )
    if require_binary and not readiness.binary:
        raise ComparativeMethodError(
            "tree must be strictly binary for this comparative method"
        )
    if len(readiness.analysis_taxa) < minimum_taxa:
        raise ComparativeMethodError(
            f"this comparative method requires at least {minimum_taxa} taxa",
            details={
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "failure_reason": "comparative_taxon_overlap_insufficient",
                "scientific_explanation": (
                    "Too few taxa remain after matching the tree and trait table, so the comparative fit would not be scientifically identifiable."
                ),
                "likely_causes": [
                    "too many tree taxa are missing from the trait table",
                    "too many overlapping taxa have missing or non-numeric trait values",
                ],
                "actionable_fixes": [
                    "add trait values for the missing tree taxa or intentionally prune the tree",
                    "repair missing or non-numeric trait values before rerunning the comparative method",
                ],
                "evidence": {
                    "minimum_taxa": minimum_taxa,
                    "analysis_taxa": readiness.analysis_taxa,
                    "missing_from_traits": readiness.missing_from_traits,
                    "extra_trait_taxa": readiness.extra_trait_taxa,
                    "pruned_missing_value_taxa": readiness.pruned_missing_value_taxa,
                    "pruned_non_numeric_taxa": readiness.pruned_non_numeric_taxa,
                },
            },
        )

    tree, _pruning_report = prune_tree_to_requested_taxa(
        tree_path,
        readiness.analysis_taxa,
    )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    values_by_taxon = {
        row[table.taxon_column]: float(row[trait])
        for row in table.rows
        if row[table.taxon_column] in readiness.analysis_taxa and row[trait]
    }
    taxa = list(readiness.analysis_taxa)
    covariance_matrix = stable_covariance(build_brownian_covariance_matrix(tree, taxa))
    return ComparativeDataset(
        tree_path=tree_path,
        traits_path=traits_path,
        tree=tree,
        taxon_column=table.taxon_column,
        trait=trait,
        taxa=taxa,
        trait_values=[values_by_taxon[taxon] for taxon in taxa],
        covariance_matrix=covariance_matrix,
        readiness=readiness,
    )


def build_brownian_covariance_matrix(
    tree: PhyloTree, taxa: list[str]
) -> list[list[float]]:
    """Build the Brownian shared-path covariance matrix for an ordered tip list."""
    leaf_paths = _leaf_ancestor_depths(tree)
    matrix: list[list[float]] = []
    for left_taxon in taxa:
        left_path = leaf_paths[left_taxon]
        row: list[float] = []
        for right_taxon in taxa:
            right_path = leaf_paths[right_taxon]
            shared_ancestor_ids = set(left_path) & set(right_path)
            shared_depth = max(left_path[node_id] for node_id in shared_ancestor_ids)
            row.append(shared_depth)
        matrix.append(row)
    return matrix


def tip_root_depths(tree: PhyloTree, taxa: list[str]) -> dict[str, float]:
    """Return root-to-tip path lengths for an explicit ordered taxon set."""
    leaf_paths = _leaf_ancestor_depths(tree)
    return {taxon: max(leaf_paths[taxon].values()) for taxon in taxa}


def build_ou_covariance_matrix(
    tree: PhyloTree, taxa: list[str], *, alpha: float
) -> list[list[float]]:
    """Build a stationary-root OU covariance matrix for an ordered tip list."""
    if alpha <= 0.0:
        raise ComparativeMethodError("OU alpha must be positive")
    leaf_paths = _leaf_ancestor_depths(tree)
    root_depths = tip_root_depths(tree, taxa)
    matrix: list[list[float]] = []
    for left_taxon in taxa:
        left_path = leaf_paths[left_taxon]
        left_depth = root_depths[left_taxon]
        row: list[float] = []
        for right_taxon in taxa:
            right_path = leaf_paths[right_taxon]
            right_depth = root_depths[right_taxon]
            shared_ancestor_ids = set(left_path) & set(right_path)
            shared_depth = max(left_path[node_id] for node_id in shared_ancestor_ids)
            if left_taxon == right_taxon:
                covariance = (1.0 - math.exp(-2.0 * alpha * left_depth)) / (2.0 * alpha)
            else:
                covariance = (
                    math.exp(
                        -alpha
                        * ((left_depth - shared_depth) + (right_depth - shared_depth))
                    )
                    * (1.0 - math.exp(-2.0 * alpha * shared_depth))
                    / (2.0 * alpha)
                )
            row.append(covariance)
        matrix.append(row)
    return matrix


def descendant_taxa(node: TreeNode) -> list[str]:
    """Return the sorted descendant-tip names for a node."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(descendant_taxa(child))
    return sorted(taxa)


def node_signature(node: TreeNode) -> str:
    """Return a deterministic signature for one node."""
    taxa = descendant_taxa(node)
    return "|".join(taxa) if taxa else node.name or "<unnamed>"


def tip_index_lookup(tree: PhyloTree) -> dict[str, TreeNode]:
    """Index tree leaves by exact taxon label."""
    lookup: dict[str, TreeNode] = {}
    for leaf in tree.iter_leaves():
        if leaf.name is not None:
            lookup[leaf.name] = leaf
    return lookup


def lambda_transform_covariance(
    covariance_matrix: list[list[float]],
    lambda_value: float,
) -> list[list[float]]:
    """Apply Pagel's lambda transformation to a covariance matrix."""
    transformed: list[list[float]] = []
    for row_index, row in enumerate(covariance_matrix):
        transformed_row: list[float] = []
        for column_index, value in enumerate(row):
            if row_index == column_index:
                transformed_row.append(value)
            else:
                transformed_row.append(value * lambda_value)
        transformed.append(transformed_row)
    return stable_covariance(transformed)


def _has_complete_branch_lengths(tree: PhyloTree) -> bool:
    return all(
        node.branch_length is not None
        for node in tree.iter_nodes()
        if node is not tree.root
    )


def _minimum_branch_length(tree: PhyloTree) -> float | None:
    branch_lengths = [
        node.branch_length
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    if not branch_lengths:
        return None
    return min(branch_lengths)


def _leaf_ancestor_depths(tree: PhyloTree) -> dict[str, dict[str, float]]:
    depths_by_leaf: dict[str, dict[str, float]] = {}

    def visit(node: TreeNode, ancestors: dict[str, float], depth: float) -> None:
        if node is not tree.root:
            if node.branch_length is None:
                raise ComparativeMethodError(
                    "tree must contain complete branch lengths for comparative analysis"
                )
            depth += node.branch_length
        current_ancestors = dict(ancestors)
        current_ancestors[node.node_id or node_signature(node)] = depth
        if node.is_leaf():
            if node.name is None:
                raise ComparativeMethodError("tree contains an unnamed terminal taxon")
            depths_by_leaf[node.name] = current_ancestors
            return
        for child in node.children:
            visit(child, current_ancestors, depth)

    visit(tree.root, {}, 0.0)
    return depths_by_leaf
