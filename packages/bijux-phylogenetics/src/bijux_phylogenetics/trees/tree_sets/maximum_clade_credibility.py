from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .contracts import (
    TreeSetMaximumCladeCredibilityCandidateRow,
    TreeSetMaximumCladeCredibilityReport,
)
from .inventory import _analyze_tree_set, _require_exact_taxa, _TreeSetAnalysis


def _score_candidate_tree(
    tree: PhyloTree,
    *,
    shared_taxa: set[str],
    clade_frequencies: dict[frozenset[str], float],
) -> float:
    # Sum log clade frequencies for numerical stability instead of multiplying.
    return round(
        sum(
            math.log(clade_frequencies[clade])
            for clade in informative_rooted_clades(tree, shared_taxa)
        ),
        15,
    )


def _compute_maximum_clade_credibility_from_analysis(
    analysis: _TreeSetAnalysis,
) -> tuple[PhyloTree, TreeSetMaximumCladeCredibilityReport]:
    exact_taxa = _require_exact_taxa(analysis)
    clade_counts = analysis.clade_counts or {}
    tree_count = len(analysis.trees)
    shared_taxa = set(exact_taxa)
    clade_frequencies = {
        clade: count / tree_count for clade, count in clade_counts.items()
    }

    raw_tree_counts: dict[str, int] = {}
    for tree in analysis.trees:
        newick = dumps_newick(tree)
        raw_tree_counts[newick] = raw_tree_counts.get(newick, 0) + 1

    scored_rows: list[
        tuple[float, int, TreeSetMaximumCladeCredibilityCandidateRow]
    ] = []
    for record, tree in zip(analysis.records, analysis.trees, strict=True):
        candidate_newick = dumps_newick(tree)
        score = _score_candidate_tree(
            tree,
            shared_taxa=shared_taxa,
            clade_frequencies=clade_frequencies,
        )
        scored_rows.append(
            (
                score,
                record.index,
                TreeSetMaximumCladeCredibilityCandidateRow(
                    score_rank=0,
                    source_tree_index=record.index,
                    rooted_topology_id=record.rooted_topology_id,
                    raw_tree_count=raw_tree_counts[candidate_newick],
                    raw_tree_frequency=round(
                        raw_tree_counts[candidate_newick] / tree_count, 15
                    ),
                    clade_credibility_score=score,
                    candidate_newick=candidate_newick,
                ),
            )
        )

    ranked_rows = sorted(scored_rows, key=lambda item: (-item[0], item[1]))
    selected_score, selected_index, selected_row = ranked_rows[0]
    selected_tree = analysis.trees[selected_index - 1]
    selected_rooted_topology_id = analysis.records[
        selected_index - 1
    ].rooted_topology_id
    rows = [
        TreeSetMaximumCladeCredibilityCandidateRow(
            score_rank=rank,
            source_tree_index=row.source_tree_index,
            rooted_topology_id=row.rooted_topology_id,
            raw_tree_count=row.raw_tree_count,
            raw_tree_frequency=row.raw_tree_frequency,
            clade_credibility_score=row.clade_credibility_score,
            candidate_newick=row.candidate_newick,
        )
        for rank, (_score, _index, row) in enumerate(ranked_rows, start=1)
    ]
    report = TreeSetMaximumCladeCredibilityReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        rooted_topology_count=len(analysis.rooted_topology_counts),
        selected_tree_index=selected_index,
        selected_rooted_topology_id=selected_rooted_topology_id,
        clade_credibility_score=selected_score,
        maximum_clade_credibility_newick=selected_row.candidate_newick,
        rows=rows,
    )
    return selected_tree, report


def compute_maximum_clade_credibility_tree(
    path: Path,
) -> tuple[PhyloTree, TreeSetMaximumCladeCredibilityReport]:
    """Select the sampled tree with the highest summed posterior clade credibility."""
    return _compute_maximum_clade_credibility_from_analysis(_analyze_tree_set(path))


def write_maximum_clade_credibility_score_table(
    path: Path,
    report: TreeSetMaximumCladeCredibilityReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "score_rank",
                "source_tree_index",
                "rooted_topology_id",
                "raw_tree_count",
                "raw_tree_frequency",
                "clade_credibility_score",
                "candidate_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "score_rank": row.score_rank,
                    "source_tree_index": row.source_tree_index,
                    "rooted_topology_id": row.rooted_topology_id,
                    "raw_tree_count": row.raw_tree_count,
                    "raw_tree_frequency": format(row.raw_tree_frequency, ".15g"),
                    "clade_credibility_score": format(
                        row.clade_credibility_score, ".15g"
                    ),
                    "candidate_newick": row.candidate_newick,
                }
            )
    return path


def write_maximum_clade_credibility_artifacts(
    out_dir: Path,
    tree: PhyloTree,
    report: TreeSetMaximumCladeCredibilityReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_path = write_newick(out_dir / "maximum-clade-credibility-tree.nwk", tree)
    score_table_path = write_maximum_clade_credibility_score_table(
        out_dir / "candidate-score-table.tsv",
        report,
    )
    return {
        "maximum_clade_credibility_tree_path": tree_path,
        "candidate_score_table_path": score_table_path,
    }
