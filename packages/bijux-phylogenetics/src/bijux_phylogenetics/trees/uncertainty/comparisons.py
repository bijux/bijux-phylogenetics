from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.iqtree_support import (
    parse_iqtree_branch_support_label,
)
from bijux_phylogenetics.io.iqtree_support import (
    support_fraction as normalize_support_fraction,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clade_nodes
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from ..tree_sets.clade_support import compute_clade_frequency_table
from ..tree_sets.inventory import (
    _require_tree_set,
    _validate_same_taxa,
)
from ..tree_sets.topology import (
    _clade_counts,
    _format_clade,
    _rooted_topology_id,
    _tree_distance,
)
from .models import (
    BootstrapPosteriorCladeComparison,
    BootstrapPosteriorSupportComparisonReport,
    CladeFrequencyDelta,
    PosteriorTreeSetComparisonReport,
)


def compare_posterior_tree_sets(
    left_path: Path,
    right_path: Path,
) -> PosteriorTreeSetComparisonReport:
    """Compare two tree sets over shared taxa, clade support, and cross-set topology distance."""
    _, left_trees = _require_tree_set(left_path)
    _, right_trees = _require_tree_set(right_path)
    left_taxa = _validate_same_taxa(left_trees)
    right_taxa = _validate_same_taxa(right_trees)
    if left_taxa != right_taxa:
        raise InvalidAlignmentError(
            "posterior tree-set comparison requires identical taxon sets across both inputs"
        )

    shared_taxa = set(left_taxa)
    left_counts = _clade_counts(left_trees, shared_taxa)
    right_counts = _clade_counts(right_trees, shared_taxa)
    all_clades = left_counts.keys() | right_counts.keys()
    deltas = [
        CladeFrequencyDelta(
            clade=_format_clade(clade),
            left_frequency=round(left_counts.get(clade, 0) / len(left_trees), 15),
            right_frequency=round(right_counts.get(clade, 0) / len(right_trees), 15),
            delta=round(
                (right_counts.get(clade, 0) / len(right_trees))
                - (left_counts.get(clade, 0) / len(left_trees)),
                15,
            ),
        )
        for clade in sorted(all_clades, key=_format_clade)
    ]
    comparisons = [
        _tree_distance(left, right, shared_taxa)
        for left in left_trees
        for right in right_trees
    ]
    left_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in left_trees}
    right_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in right_trees}
    return PosteriorTreeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=left_taxa,
        left_tree_count=len(left_trees),
        right_tree_count=len(right_trees),
        left_rooted_topology_count=len(left_topologies),
        right_rooted_topology_count=len(right_topologies),
        shared_rooted_topology_count=len(left_topologies & right_topologies),
        mean_between_set_robinson_foulds=round(
            sum(distance for distance, _ in comparisons) / len(comparisons),
            15,
        ),
        mean_between_set_normalized_robinson_foulds=round(
            sum(normalized for _, normalized in comparisons) / len(comparisons),
            15,
        ),
        clade_frequency_deltas=deltas,
    )


def _require_tree(path: Path) -> PhyloTree:
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    return load_tree(path)


def _parse_support_label(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    parsed_iqtree_label = parse_iqtree_branch_support_label(text)
    if parsed_iqtree_label is not None:
        support_value = (
            parsed_iqtree_label.ufboot_support
            if parsed_iqtree_label.ufboot_support is not None
            else parsed_iqtree_label.sh_alrt_support
        )
        if support_value is None:
            return None
        normalized = normalize_support_fraction(support_value)
        return None if normalized is None else round(normalized, 15)
    try:
        parsed = float(text)
    except ValueError:
        return None
    return round(parsed / 100.0, 15) if parsed > 1.0 else round(parsed, 15)


def _support_agreement_label(
    bootstrap_support: float | None,
    posterior_frequency: float | None,
    absolute_delta: float | None,
) -> str:
    if bootstrap_support is None and posterior_frequency is None:
        return "not_observed"
    if bootstrap_support is None or posterior_frequency is None:
        return "method_specific"
    if absolute_delta is None:
        return "not_comparable"
    if absolute_delta >= 0.35:
        return "strong_conflict"
    if absolute_delta >= 0.15:
        return "moderate_difference"
    return "broad_agreement"


def compare_bootstrap_and_posterior_uncertainty(
    bootstrap_tree_path: Path,
    posterior_tree_set_path: Path,
) -> BootstrapPosteriorSupportComparisonReport:
    """Compare bootstrap support on one summary tree against posterior clade frequencies from a tree set."""
    bootstrap_tree = _require_tree(bootstrap_tree_path)
    posterior_report = compute_clade_frequency_table(posterior_tree_set_path)
    shared_taxa = set(bootstrap_tree.tip_names)
    posterior_taxa = set(posterior_report.shared_taxa)
    if shared_taxa != posterior_taxa:
        raise InvalidAlignmentError(
            "bootstrap versus posterior comparison requires identical taxon sets"
        )
    bootstrap_nodes = informative_rooted_clade_nodes(bootstrap_tree, shared_taxa)
    bootstrap_support_by_clade = {
        _format_clade(clade): _parse_support_label(node.name)
        for clade, node in bootstrap_nodes.items()
    }
    posterior_frequency_by_clade = {
        row.clade: row.frequency for row in posterior_report.clade_frequencies
    }
    all_clades = sorted(
        set(bootstrap_support_by_clade) | set(posterior_frequency_by_clade)
    )
    rows: list[BootstrapPosteriorCladeComparison] = []
    for clade in all_clades:
        bootstrap_support = bootstrap_support_by_clade.get(clade)
        posterior_frequency = posterior_frequency_by_clade.get(clade)
        absolute_delta = None
        if bootstrap_support is not None and posterior_frequency is not None:
            absolute_delta = abs(bootstrap_support - posterior_frequency)
        agreement = _support_agreement_label(
            bootstrap_support,
            posterior_frequency,
            absolute_delta,
        )
        rows.append(
            BootstrapPosteriorCladeComparison(
                clade=clade,
                bootstrap_support=bootstrap_support,
                posterior_frequency=posterior_frequency,
                absolute_delta=absolute_delta,
                agreement=agreement,
            )
        )
    topology_mismatch_clade_count = sum(
        1 for row in rows if row.agreement == "method_specific"
    )
    return BootstrapPosteriorSupportComparisonReport(
        bootstrap_tree_path=bootstrap_tree_path,
        posterior_tree_set_path=posterior_tree_set_path,
        posterior_tree_count=posterior_report.tree_count,
        shared_taxa=posterior_report.shared_taxa,
        high_conflict_clade_count=sum(
            1 for row in rows if row.agreement == "strong_conflict"
        ),
        topology_mismatch_detected=topology_mismatch_clade_count > 0,
        topology_mismatch_clade_count=topology_mismatch_clade_count,
        rows=rows,
    )
