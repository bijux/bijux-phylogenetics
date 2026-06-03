from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.compare.topology.branch_lengths import (
    _compare_branch_lengths_for_trees,
)
from bijux_phylogenetics.compare.topology.distance import (
    compare_topology_distance_trees,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .consensus import _build_consensus_tree_with_threshold
from .contracts import (
    PosteriorTreeDistanceDiagnosticRow,
    PosteriorTreeDistanceDiagnosticsReport,
    PosteriorTreeDistanceDistributionRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa
from .maximum_clade_credibility import _compute_maximum_clade_credibility_from_analysis

_MCC_REFERENCE_KIND = "maximum-clade-credibility"
_CONSENSUS_REFERENCE_KIND = "consensus"


def _distance_distribution_rows(
    rows: list[PosteriorTreeDistanceDiagnosticRow],
) -> list[PosteriorTreeDistanceDistributionRow]:
    counts: dict[tuple[str, str, float], int] = {}
    total_tree_count = len(rows)
    for row in rows:
        observations = [
            (
                _MCC_REFERENCE_KIND,
                "robinson-foulds",
                float(row.mcc_robinson_foulds_distance),
            ),
            (
                _MCC_REFERENCE_KIND,
                "normalized-robinson-foulds",
                row.mcc_normalized_robinson_foulds,
            ),
            (
                _CONSENSUS_REFERENCE_KIND,
                "robinson-foulds",
                float(row.consensus_robinson_foulds_distance),
            ),
            (
                _CONSENSUS_REFERENCE_KIND,
                "normalized-robinson-foulds",
                row.consensus_normalized_robinson_foulds,
            ),
        ]
        if row.mcc_branch_score_distance is not None:
            observations.append(
                (_MCC_REFERENCE_KIND, "branch-score", row.mcc_branch_score_distance)
            )
        if row.consensus_branch_score_distance is not None:
            observations.append(
                (
                    _CONSENSUS_REFERENCE_KIND,
                    "branch-score",
                    row.consensus_branch_score_distance,
                )
            )
        for reference_tree_kind, distance_metric, observed_value in observations:
            key = (reference_tree_kind, distance_metric, observed_value)
            counts[key] = counts.get(key, 0) + 1
    return [
        PosteriorTreeDistanceDistributionRow(
            reference_tree_kind=reference_tree_kind,
            distance_metric=distance_metric,
            observed_value=observed_value,
            tree_count=count,
            frequency=round(count / total_tree_count, 15),
        )
        for (reference_tree_kind, distance_metric, observed_value), count in sorted(
            counts.items(),
            key=lambda item: (item[0][0], item[0][1], item[0][2]),
        )
    ]


def _rank_rows(
    rows: list[dict[str, float | int | str | None]],
    *,
    normalized_key: str,
    branch_score_key: str,
    rf_key: str,
) -> dict[int, int]:
    ranked = sorted(
        rows,
        key=lambda row: (
            -float(row[normalized_key]),
            -(-1.0 if row[branch_score_key] is None else float(row[branch_score_key])),
            -int(row[rf_key]),
            int(row["source_tree_index"]),
        ),
    )
    return {
        int(row["source_tree_index"]): rank for rank, row in enumerate(ranked, start=1)
    }


def _diagnostic_row(
    source_tree_index: int,
    rooted_topology_id: str,
    tree: PhyloTree,
    mcc_tree: PhyloTree,
    consensus_tree: PhyloTree,
    *,
    left_path: Path,
    right_path: Path,
) -> dict[str, float | int | str | None]:
    mcc_topology = compare_topology_distance_trees(
        tree,
        mcc_tree,
        left_path=left_path,
        right_path=right_path,
    )
    mcc_branch_lengths = _compare_branch_lengths_for_trees(
        left_path,
        right_path,
        tree,
        mcc_tree,
        taxon_overlap_policy="require-identical",
    )
    consensus_topology = compare_topology_distance_trees(
        tree,
        consensus_tree,
        left_path=left_path,
        right_path=right_path,
    )
    consensus_branch_lengths = _compare_branch_lengths_for_trees(
        left_path,
        right_path,
        tree,
        consensus_tree,
        taxon_overlap_policy="require-identical",
    )
    return {
        "source_tree_index": source_tree_index,
        "rooted_topology_id": rooted_topology_id,
        "mcc_robinson_foulds_distance": mcc_topology.robinson_foulds_distance,
        "mcc_normalized_robinson_foulds": mcc_topology.normalized_robinson_foulds,
        "mcc_branch_score_distance": mcc_branch_lengths.branch_score_distance,
        "consensus_robinson_foulds_distance": (
            consensus_topology.robinson_foulds_distance
        ),
        "consensus_normalized_robinson_foulds": (
            consensus_topology.normalized_robinson_foulds
        ),
        "consensus_branch_score_distance": (
            consensus_branch_lengths.branch_score_distance
        ),
    }


def compute_posterior_tree_distance_diagnostics(
    path: Path,
) -> PosteriorTreeDistanceDiagnosticsReport:
    """Compare every posterior tree against MCC and consensus references."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    mcc_tree, mcc_report = _compute_maximum_clade_credibility_from_analysis(analysis)
    consensus_tree, consensus_report = _build_consensus_tree_with_threshold(
        analysis,
        threshold=0.5,
    )

    raw_rows = [
        _diagnostic_row(
            record.index,
            record.rooted_topology_id,
            tree,
            mcc_tree,
            consensus_tree,
            left_path=analysis.path,
            right_path=analysis.path,
        )
        for record, tree in zip(analysis.records, analysis.trees, strict=True)
    ]
    mcc_ranks = _rank_rows(
        raw_rows,
        normalized_key="mcc_normalized_robinson_foulds",
        branch_score_key="mcc_branch_score_distance",
        rf_key="mcc_robinson_foulds_distance",
    )
    consensus_ranks = _rank_rows(
        raw_rows,
        normalized_key="consensus_normalized_robinson_foulds",
        branch_score_key="consensus_branch_score_distance",
        rf_key="consensus_robinson_foulds_distance",
    )
    rows = [
        PosteriorTreeDistanceDiagnosticRow(
            source_tree_index=int(row["source_tree_index"]),
            rooted_topology_id=str(row["rooted_topology_id"]),
            mcc_outlier_rank=mcc_ranks[int(row["source_tree_index"])],
            consensus_outlier_rank=consensus_ranks[int(row["source_tree_index"])],
            mcc_robinson_foulds_distance=int(row["mcc_robinson_foulds_distance"]),
            mcc_normalized_robinson_foulds=float(row["mcc_normalized_robinson_foulds"]),
            mcc_branch_score_distance=(
                None
                if row["mcc_branch_score_distance"] is None
                else float(row["mcc_branch_score_distance"])
            ),
            consensus_robinson_foulds_distance=int(
                row["consensus_robinson_foulds_distance"]
            ),
            consensus_normalized_robinson_foulds=float(
                row["consensus_normalized_robinson_foulds"]
            ),
            consensus_branch_score_distance=(
                None
                if row["consensus_branch_score_distance"] is None
                else float(row["consensus_branch_score_distance"])
            ),
        )
        for row in sorted(raw_rows, key=lambda row: int(row["source_tree_index"]))
    ]
    distributions = _distance_distribution_rows(rows)
    return PosteriorTreeDistanceDiagnosticsReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        maximum_clade_credibility_tree_index=mcc_report.selected_tree_index,
        maximum_clade_credibility_rooted_topology_id=(
            mcc_report.selected_rooted_topology_id
        ),
        maximum_clade_credibility_newick=mcc_report.maximum_clade_credibility_newick,
        consensus_method=consensus_report.consensus_method,
        consensus_newick=consensus_report.consensus_newick,
        row_count=len(rows),
        distribution_row_count=len(distributions),
        rows=rows,
        distributions=distributions,
    )


def write_posterior_tree_distance_diagnostic_table(
    path: Path,
    report: PosteriorTreeDistanceDiagnosticsReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_tree_index",
                "rooted_topology_id",
                "mcc_outlier_rank",
                "consensus_outlier_rank",
                "mcc_robinson_foulds_distance",
                "mcc_normalized_robinson_foulds",
                "mcc_branch_score_distance",
                "consensus_robinson_foulds_distance",
                "consensus_normalized_robinson_foulds",
                "consensus_branch_score_distance",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "source_tree_index": row.source_tree_index,
                    "rooted_topology_id": row.rooted_topology_id,
                    "mcc_outlier_rank": row.mcc_outlier_rank,
                    "consensus_outlier_rank": row.consensus_outlier_rank,
                    "mcc_robinson_foulds_distance": row.mcc_robinson_foulds_distance,
                    "mcc_normalized_robinson_foulds": format(
                        row.mcc_normalized_robinson_foulds,
                        ".15g",
                    ),
                    "mcc_branch_score_distance": (
                        ""
                        if row.mcc_branch_score_distance is None
                        else format(row.mcc_branch_score_distance, ".15g")
                    ),
                    "consensus_robinson_foulds_distance": (
                        row.consensus_robinson_foulds_distance
                    ),
                    "consensus_normalized_robinson_foulds": format(
                        row.consensus_normalized_robinson_foulds,
                        ".15g",
                    ),
                    "consensus_branch_score_distance": (
                        ""
                        if row.consensus_branch_score_distance is None
                        else format(row.consensus_branch_score_distance, ".15g")
                    ),
                }
            )
    return path


def write_posterior_tree_distance_distribution_table(
    path: Path,
    report: PosteriorTreeDistanceDiagnosticsReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reference_tree_kind",
                "distance_metric",
                "observed_value",
                "tree_count",
                "frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.distributions:
            writer.writerow(
                {
                    "reference_tree_kind": row.reference_tree_kind,
                    "distance_metric": row.distance_metric,
                    "observed_value": format(row.observed_value, ".15g"),
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                }
            )
    return path


def write_posterior_tree_distance_artifacts(
    out_dir: Path,
    mcc_tree: PhyloTree,
    consensus_tree: PhyloTree,
    report: PosteriorTreeDistanceDiagnosticsReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    mcc_tree_path = write_newick(
        out_dir / "maximum-clade-credibility-tree.nwk", mcc_tree
    )
    consensus_tree_path = write_newick(out_dir / "consensus-tree.nwk", consensus_tree)
    diagnostic_table_path = write_posterior_tree_distance_diagnostic_table(
        out_dir / "posterior-tree-distance-diagnostics.tsv",
        report,
    )
    distribution_table_path = write_posterior_tree_distance_distribution_table(
        out_dir / "posterior-tree-distance-distribution.tsv",
        report,
    )
    return {
        "maximum_clade_credibility_tree_path": mcc_tree_path,
        "consensus_tree_path": consensus_tree_path,
        "posterior_tree_distance_diagnostic_table_path": diagnostic_table_path,
        "posterior_tree_distance_distribution_table_path": distribution_table_path,
    }
