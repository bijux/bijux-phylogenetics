from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.compare.influence import TaxonInfluenceReport
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    InferenceComparisonConclusionRow,
    InferenceComparisonConclusionSummary,
    InferenceComparisonWeightedConflictRow,
)


def rewrite_inference_comparison_report_html(
    *,
    base_report: ComparisonReportBuildResult,
    summary: InferenceComparisonConclusionSummary,
    conclusion_rows: list[InferenceComparisonConclusionRow],
    weighted_conflict_rows: list[InferenceComparisonWeightedConflictRow],
    taxon_influence_report: TaxonInfluenceReport | None,
) -> Path:
    """Rewrite the comparison HTML with reviewer-facing stability synthesis."""
    sections = [
        (
            "comparison-summary",
            json.dumps(asdict(summary), indent=2, sort_keys=True),
        ),
        (
            "clade-conclusions",
            json.dumps(
                [asdict(row) for row in conclusion_rows],
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "support-weighted-conflicts",
            json.dumps(
                [asdict(row) for row in weighted_conflict_rows],
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "taxon-influence",
            json.dumps(
                None
                if taxon_influence_report is None
                else asdict(taxon_influence_report),
                indent=2,
                sort_keys=True,
                default=str,
            ),
        ),
        (
            "topology-metrics",
            json.dumps(
                asdict(base_report.topology), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "clade-comparison",
            json.dumps(
                asdict(base_report.clades), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "support-comparison",
            json.dumps(
                asdict(base_report.support), indent=2, sort_keys=True, default=str
            ),
        ),
        (
            "branch-length-comparison",
            json.dumps(
                asdict(base_report.branch_lengths),
                indent=2,
                sort_keys=True,
                default=str,
            ),
        ),
    ]
    return write_html_report(
        title="Bijux Tree Inference Comparison Report",
        sections=sections,
        out_path=base_report.output_path,
        embedded_json={
            "summary": asdict(summary),
            "conclusions": [asdict(row) for row in conclusion_rows],
            "weighted_conflicts": [asdict(row) for row in weighted_conflict_rows],
        },
    )
