from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

from .models import (
    RelaxedRateBranchOutlier,
    RelaxedRateBranchSummaryReport,
    RelaxedRateBranchSummaryRow,
)


def summarize_relaxed_rate_branches(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    *,
    substitution_tree_path: Path | None = None,
    dated_tree_path: Path | None = None,
    outlier_threshold: float = 2.0,
) -> RelaxedRateBranchSummaryReport:
    """Summarize branch-specific rates from one substitution tree and one dated tree."""
    if substitution_tree.rooted is not True:
        raise UnrootedTreeError(
            "relaxed-rate branch summary requires one rooted substitution tree",
            code="relaxed_rate_branch_summary_requires_rooted_tree",
        )
    if dated_tree.rooted is not True:
        raise UnrootedTreeError(
            "relaxed-rate branch summary requires one rooted dated tree",
            code="relaxed_rate_branch_summary_requires_rooted_tree",
        )
    if outlier_threshold <= 0.0:
        raise PhylogeneticsError(
            "relaxed-rate branch summary requires a strictly positive outlier threshold",
            code="relaxed_rate_branch_summary_error",
        )

    _require_complete_branch_lengths(
        substitution_tree,
        message="relaxed-rate branch summary requires complete substitution branch lengths",
    )
    _require_complete_branch_lengths(
        dated_tree,
        message="relaxed-rate branch summary requires complete dated branch durations",
    )
    _require_matching_topology(substitution_tree, dated_tree)

    dated_branch_by_signature = _branch_lookup_by_descendant_taxa(dated_tree)
    raw_rows: list[tuple[str, str | None, list[str], float, float, float]] = []
    for _parent, child in substitution_tree.iter_edges():
        branch_signature = tuple(child.descendant_taxa)
        matched_dated_child = dated_branch_by_signature[branch_signature]
        substitution_branch_length = float(child.branch_length or 0.0)
        dated_time_duration = float(matched_dated_child.branch_length or 0.0)
        if dated_time_duration <= 0.0:
            raise InvalidBranchLengthError(
                "relaxed-rate branch summary requires strictly positive dated branch durations"
            )
        if child.node_id is None:
            raise PhylogeneticsError(
                "relaxed-rate branch summary requires stable branch ids",
                code="relaxed_rate_branch_summary_error",
            )
        raw_rows.append(
            (
                child.node_id,
                child.name,
                child.descendant_taxa,
                substitution_branch_length,
                dated_time_duration,
                substitution_branch_length / dated_time_duration,
            )
        )

    branch_rates = [row[5] for row in raw_rows]
    mean_branch_rate = sum(branch_rates) / len(branch_rates)
    rate_sum_squares = sum(
        (branch_rate - mean_branch_rate) ** 2 for branch_rate in branch_rates
    )
    standard_deviation_branch_rate = math.sqrt(rate_sum_squares / len(branch_rates))
    branch_rows: list[RelaxedRateBranchSummaryRow] = []
    for (
        branch_id,
        child_name,
        descendant_taxa,
        substitution_branch_length,
        dated_time_duration,
        branch_rate,
    ) in raw_rows:
        rate_z_score = (
            (branch_rate - mean_branch_rate) / standard_deviation_branch_rate
            if standard_deviation_branch_rate > 0.0
            else 0.0
        )
        branch_rows.append(
            RelaxedRateBranchSummaryRow(
                branch_id=branch_id,
                child_name=child_name,
                descendant_taxa=descendant_taxa,
                substitution_branch_length=substitution_branch_length,
                dated_time_duration=dated_time_duration,
                branch_rate=branch_rate,
                rate_z_score=rate_z_score,
                outlier=abs(rate_z_score) >= outlier_threshold,
            )
        )

    ranked_outliers = sorted(
        (row for row in branch_rows if row.outlier),
        key=lambda row: (-abs(row.rate_z_score), row.branch_id),
    )
    outlier_rows = [
        RelaxedRateBranchOutlier(
            rank=index,
            branch_id=row.branch_id,
            child_name=row.child_name,
            descendant_taxa=row.descendant_taxa,
            substitution_branch_length=row.substitution_branch_length,
            dated_time_duration=row.dated_time_duration,
            branch_rate=row.branch_rate,
            rate_z_score=row.rate_z_score,
        )
        for index, row in enumerate(ranked_outliers, start=1)
    ]
    return RelaxedRateBranchSummaryReport(
        substitution_tree_newick=dumps_newick(substitution_tree),
        dated_tree_newick=dumps_newick(dated_tree),
        taxa=sorted(substitution_tree.tip_names),
        tip_count=substitution_tree.tip_count,
        internal_node_count=substitution_tree.internal_node_count,
        branch_count=len(branch_rows),
        substitution_tree_path=(
            None if substitution_tree_path is None else str(substitution_tree_path)
        ),
        dated_tree_path=None if dated_tree_path is None else str(dated_tree_path),
        outlier_threshold=outlier_threshold,
        mean_branch_rate=mean_branch_rate,
        standard_deviation_branch_rate=standard_deviation_branch_rate,
        minimum_branch_rate=min(branch_rates),
        maximum_branch_rate=max(branch_rates),
        outlier_count=len(outlier_rows),
        branch_rows=branch_rows,
        outlier_rows=outlier_rows,
    )


def summarize_relaxed_rate_branches_from_paths(
    substitution_tree_path: Path,
    dated_tree_path: Path,
    *,
    outlier_threshold: float = 2.0,
) -> RelaxedRateBranchSummaryReport:
    """Summarize branch-specific rates from one rooted tree path pair."""
    validate_tree_path(substitution_tree_path, require_rooted=True)
    validate_tree_path(dated_tree_path, require_rooted=True)
    substitution_tree = load_tree(substitution_tree_path)
    dated_tree = load_tree(dated_tree_path)
    substitution_tree.rooted = True
    dated_tree.rooted = True
    return summarize_relaxed_rate_branches(
        substitution_tree,
        dated_tree,
        substitution_tree_path=substitution_tree_path,
        dated_tree_path=dated_tree_path,
        outlier_threshold=outlier_threshold,
    )


def _require_complete_branch_lengths(tree: PhyloTree, *, message: str) -> None:
    if any(branch_length is None for branch_length in tree.branch_lengths()):
        raise InvalidBranchLengthError(message)


def _require_matching_topology(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
) -> None:
    substitution_taxa = sorted(substitution_tree.tip_names)
    dated_taxa = sorted(dated_tree.tip_names)
    if substitution_taxa != dated_taxa:
        raise PhylogeneticsError(
            "relaxed-rate branch summary requires matching tree taxa",
            code="relaxed_rate_branch_summary_error",
        )
    substitution_signatures = set(_branch_lookup_by_descendant_taxa(substitution_tree))
    dated_signatures = set(_branch_lookup_by_descendant_taxa(dated_tree))
    if substitution_signatures != dated_signatures:
        raise PhylogeneticsError(
            "relaxed-rate branch summary requires identical rooted topology between substitution and dated trees",
            code="relaxed_rate_branch_summary_error",
        )


def _branch_lookup_by_descendant_taxa(
    tree: PhyloTree,
) -> Mapping[tuple[str, ...], TreeNode]:
    return {tuple(child.descendant_taxa): child for _parent, child in tree.iter_edges()}


def write_relaxed_rate_branch_summary_tsv(
    path: Path,
    report: RelaxedRateBranchSummaryReport,
) -> Path:
    """Write one summary row for one relaxed-rate branch summary run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "substitution_tree_path",
        "dated_tree_path",
        "tip_count",
        "internal_node_count",
        "branch_count",
        "outlier_threshold",
        "mean_branch_rate",
        "standard_deviation_branch_rate",
        "minimum_branch_rate",
        "maximum_branch_rate",
        "outlier_count",
        "outlier_branch_ids",
    ]
    values = [
        report.substitution_tree_path or "",
        report.dated_tree_path or "",
        str(report.tip_count),
        str(report.internal_node_count),
        str(report.branch_count),
        format(report.outlier_threshold, ".15g"),
        format(report.mean_branch_rate, ".15g"),
        format(report.standard_deviation_branch_rate, ".15g"),
        format(report.minimum_branch_rate, ".15g"),
        format(report.maximum_branch_rate, ".15g"),
        str(report.outlier_count),
        "|".join(row.branch_id for row in report.outlier_rows),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_relaxed_rate_branch_table(
    path: Path,
    report: RelaxedRateBranchSummaryReport,
) -> Path:
    """Write one branch-level relaxed-rate row per edge in tree preorder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "child_name",
        "descendant_taxa",
        "substitution_branch_length",
        "dated_time_duration",
        "branch_rate",
        "rate_z_score",
        "outlier",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.branch_id,
                row.child_name or "",
                "|".join(row.descendant_taxa),
                format(row.substitution_branch_length, ".15g"),
                format(row.dated_time_duration, ".15g"),
                format(row.branch_rate, ".15g"),
                format(row.rate_z_score, ".15g"),
                str(row.outlier).lower(),
            ]
        )
        for row in report.branch_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_relaxed_rate_outliers_tsv(
    path: Path,
    report: RelaxedRateBranchSummaryReport,
) -> Path:
    """Write one ranked relaxed-rate outlier row per flagged branch."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "rank",
        "branch_id",
        "child_name",
        "descendant_taxa",
        "substitution_branch_length",
        "dated_time_duration",
        "branch_rate",
        "rate_z_score",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                str(row.rank),
                row.branch_id,
                row.child_name or "",
                "|".join(row.descendant_taxa),
                format(row.substitution_branch_length, ".15g"),
                format(row.dated_time_duration, ".15g"),
                format(row.branch_rate, ".15g"),
                format(row.rate_z_score, ".15g"),
            ]
        )
        for row in report.outlier_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_relaxed_rate_run_json(
    path: Path,
    report: RelaxedRateBranchSummaryReport,
) -> Path:
    """Write the full relaxed-rate branch summary report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_relaxed_rate_branch_summary_artifacts(
    out_dir: Path,
    report: RelaxedRateBranchSummaryReport,
) -> dict[str, Path]:
    """Write governed outputs for one relaxed-rate branch summary run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = write_relaxed_rate_branch_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    branch_rates_path = write_relaxed_rate_branch_table(
        out_dir / "branch_rates.tsv",
        report,
    )
    outliers_path = write_relaxed_rate_outliers_tsv(out_dir / "outliers.tsv", report)
    run_json_path = write_relaxed_rate_run_json(out_dir / "run.json", report)
    return {
        "summary_path": summary_path,
        "branch_rates_path": branch_rates_path,
        "outliers_path": outliers_path,
        "run_json_path": run_json_path,
    }
