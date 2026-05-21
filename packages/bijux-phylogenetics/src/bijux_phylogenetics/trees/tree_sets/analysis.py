from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .budgets import (
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    enforce_tree_set_tree_budget,
)
from .consensus import (
    _build_consensus_tree_with_threshold,
    compute_consensus_tree,
    compute_consensus_tree_with_threshold,
    compute_strict_consensus_tree,
    write_consensus_tree,
)
from .clade_support import (
    _build_clade_frequency_report,
    _support_classification,
    compute_clade_frequency_table,
    compute_reference_tree_clade_support,
    write_clade_frequency_table,
    write_reference_tree_clade_support_table,
)
from .contracts import ConsensusTreeReport, TreeDistanceMatrixReport, TreeDistancePair
from .contracts import TreeSetReport as TreeSetReport
from .inventory import (
    _TreeSetAnalysis,
    _analyze_tree_set,
    _require_exact_taxa,
    _require_tree_set,
    _validate_same_taxa,
    load_tree_set,
)
from .topology import (
    _clade_counts,
    _clade_signature,
    _clades_conflict,
    _format_clade,
    _rooted_topology_id,
    _tree_distance,
)


def _build_tree_distance_matrix_report(
    analysis: _TreeSetAnalysis,
) -> TreeDistanceMatrixReport:
    shared_taxa = set(_require_exact_taxa(analysis))
    pairs: list[TreeDistancePair] = []
    for left_index, left in enumerate(analysis.trees, start=1):
        for right_index, right in enumerate(
            analysis.trees[left_index - 1 :], start=left_index
        ):
            distance, normalized = _tree_distance(left, right, shared_taxa)
            pairs.append(
                TreeDistancePair(
                    left_index=left_index,
                    right_index=right_index,
                    robinson_foulds_distance=distance,
                    normalized_robinson_foulds=normalized,
                )
            )
    return TreeDistanceMatrixReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=sorted(shared_taxa),
        pairs=pairs,
    )




def compute_tree_distance_matrix(path: Path) -> TreeDistanceMatrixReport:
    """Compute a pairwise RF-distance matrix across a tree set."""
    return _build_tree_distance_matrix_report(_analyze_tree_set(path))


def write_tree_distance_matrix(path: Path, report: TreeDistanceMatrixReport) -> Path:
    """Write a pairwise tree-distance matrix as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_index",
                "right_index",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.pairs:
            writer.writerow(
                {
                    "left_index": row.left_index,
                    "right_index": row.right_index,
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": format(
                        row.normalized_robinson_foulds, ".15g"
                    ),
                }
            )
    return path
