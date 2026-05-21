from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.influence import TaxonInfluenceReport

from .contracts import (
    InferenceComparisonConclusionRow,
    InferenceComparisonConclusionSummary,
    InferenceComparisonConflictRow,
    InferenceComparisonSharedCladeRow,
    InferenceComparisonWeightedConflictRow,
)


def _render_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def write_inference_comparison_taxon_influence_table(
    path: Path,
    report: TaxonInfluenceReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "influence_rank",
                "taxon",
                "retained_taxa",
                "baseline_support_disagreements",
                "leave_one_out_support_disagreements",
                "support_disagreement_delta",
                "baseline_conflicting_clades",
                "leave_one_out_conflicting_clades",
                "conflicting_clade_delta",
                "baseline_high_support_conflicts",
                "leave_one_out_high_support_conflicts",
                "high_support_conflict_delta",
                "topology_changed",
                "support_changed",
                "influence_score",
            ]
        )
    ]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.influence_rank),
                    row.taxon,
                    "|".join(row.retained_taxa),
                    str(row.baseline_support_disagreements),
                    str(row.leave_one_out_support_disagreements),
                    str(row.support_disagreement_delta),
                    str(row.baseline_conflicting_clades),
                    str(row.leave_one_out_conflicting_clades),
                    str(row.conflicting_clade_delta),
                    str(row.baseline_high_support_conflicts),
                    str(row.leave_one_out_high_support_conflicts),
                    str(row.high_support_conflict_delta),
                    "true" if row.topology_changed else "false",
                    "true" if row.support_changed else "false",
                    _render_float(row.influence_score),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_weighted_conflict_table(
    path: Path,
    rows: list[InferenceComparisonWeightedConflictRow],
) -> Path:
    """Write one ranked support-weighted conflict ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "split_id",
                "comparison_status",
                "conflict_kind",
                "severity_class",
                "fasttree_support_fraction",
                "iqtree_support_fraction",
                "support_fraction_delta",
                "strongest_support_fraction",
                "support_weight",
                "serious_conflict",
                "detail",
            ]
        )
    ]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row.split_id,
                    row.comparison_status,
                    row.conflict_kind,
                    row.severity_class,
                    _render_float(row.fasttree_support_fraction),
                    _render_float(row.iqtree_support_fraction),
                    _render_float(row.support_fraction_delta),
                    _render_float(row.strongest_support_fraction),
                    _render_float(row.support_weight),
                    "true" if row.serious_conflict else "false",
                    row.detail,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_conclusion_table(
    path: Path,
    rows: list[InferenceComparisonConclusionRow],
) -> Path:
    """Write one reviewer-facing clade stability ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "split_id",
                "conclusion_class",
                "evidence_class",
                "comparison_status",
                "fasttree_present",
                "iqtree_present",
                "fasttree_support_fraction",
                "iqtree_support_fraction",
                "support_fraction_delta",
                "serious_conflict",
                "detail",
            ]
        )
    ]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row.split_id,
                    row.conclusion_class,
                    row.evidence_class,
                    row.comparison_status,
                    "true" if row.fasttree_present else "false",
                    "true" if row.iqtree_present else "false",
                    _render_float(row.fasttree_support_fraction),
                    _render_float(row.iqtree_support_fraction),
                    _render_float(row.support_fraction_delta),
                    "true" if row.serious_conflict else "false",
                    row.detail,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_summary_table(
    path: Path,
    summary: InferenceComparisonConclusionSummary,
) -> Path:
    """Write one compact summary row for the compared inference conclusions."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "shared_taxa_count",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
                "branch_score_distance",
                "stable_clade_count",
                "unstable_clade_count",
                "engine_specific_clade_count",
                "support_weighted_conflict_count",
                "low_support_disagreement_count",
                "moderate_support_disagreement_count",
                "high_support_conflict_count",
                "high_support_disagreement_count",
                "serious_conflict_count",
                "top_conflict_driver_taxa",
            ]
        ),
        "\t".join(
            [
                str(summary.shared_taxa_count),
                str(summary.robinson_foulds_distance),
                _render_float(summary.normalized_robinson_foulds),
                _render_float(summary.branch_score_distance),
                str(summary.stable_clade_count),
                str(summary.unstable_clade_count),
                str(summary.engine_specific_clade_count),
                str(summary.support_weighted_conflict_count),
                str(summary.low_support_disagreement_count),
                str(summary.moderate_support_disagreement_count),
                str(summary.high_support_conflict_count),
                str(summary.high_support_disagreement_count),
                str(summary.serious_conflict_count),
                "|".join(summary.top_conflict_driver_taxa),
            ]
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_inference_comparison_clade_table(
    path: Path,
    *,
    shared_rows: list[InferenceComparisonSharedCladeRow] | None = None,
    conflict_rows: list[InferenceComparisonConflictRow] | None = None,
) -> Path:
    """Write one shared-clade or conflicting-clade comparison table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str]
    if shared_rows is not None:
        lines = [
            "\t".join(
                [
                    "split_id",
                    "fasttree_support",
                    "fasttree_support_fraction",
                    "fasttree_support_label_kind",
                    "iqtree_support",
                    "iqtree_support_fraction",
                    "iqtree_support_label_kind",
                    "support_fraction_delta",
                    "support_disagreement",
                ]
            )
        ]
        for row in shared_rows:
            lines.append(
                "\t".join(
                    [
                        row.split_id,
                        _render_float(row.fasttree_support),
                        _render_float(row.fasttree_support_fraction),
                        row.fasttree_support_label_kind,
                        _render_float(row.iqtree_support),
                        _render_float(row.iqtree_support_fraction),
                        row.iqtree_support_label_kind,
                        _render_float(row.support_fraction_delta),
                        "true" if row.support_disagreement else "false",
                    ]
                )
            )
    elif conflict_rows is not None:
        lines = [
            "\t".join(
                [
                    "split_id",
                    "conflict_kind",
                    "fasttree_present",
                    "iqtree_present",
                    "fasttree_support",
                    "fasttree_support_fraction",
                    "iqtree_support",
                    "iqtree_support_fraction",
                    "detail",
                ]
            )
        ]
        for row in conflict_rows:
            lines.append(
                "\t".join(
                    [
                        row.split_id,
                        row.conflict_kind,
                        "true" if row.fasttree_present else "false",
                        "true" if row.iqtree_present else "false",
                        _render_float(row.fasttree_support),
                        _render_float(row.fasttree_support_fraction),
                        _render_float(row.iqtree_support),
                        _render_float(row.iqtree_support_fraction),
                        row.detail,
                    ]
                )
            )
    else:
        raise ValueError("one row set must be provided when writing comparison tables")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
