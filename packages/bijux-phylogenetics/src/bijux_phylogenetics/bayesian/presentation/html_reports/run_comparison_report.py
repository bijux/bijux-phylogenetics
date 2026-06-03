from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.posterior_sets.comparison import (
    compare_independent_bayesian_runs,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report

from .contracts import BayesianRunComparisonReportBuildResult
from .report_policy import (
    method_tier_section,
    method_tier_summary_metrics,
    run_comparison_limitations,
)


def render_bayesian_run_comparison_report(
    *,
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    left_trace_path: Path,
    right_trace_path: Path,
    out_path: Path,
    trace_kind: str = "mrbayes",
    burnin_fraction: float = 0.25,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BayesianRunComparisonReportBuildResult:
    """Render a deterministic HTML report comparing two Bayesian runs."""
    comparison = compare_independent_bayesian_runs(
        left_tree_set_path,
        right_tree_set_path,
        left_trace_path=left_trace_path,
        right_trace_path=right_trace_path,
        trace_kind=trace_kind,
        burnin_fraction=burnin_fraction,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    title = "Bijux Bayesian Run Comparison Report"
    method_tier = bayesian_report_method_tier("bayesian-run-comparison")
    limitations = run_comparison_limitations(comparison.warnings)
    sections = [
        method_tier_section(method_tier),
        (
            "run-comparison",
            json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
        ),
        (
            "tree-comparison",
            json.dumps(
                asdict(comparison.tree_comparison),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "left-convergence",
            json.dumps(
                asdict(comparison.left_convergence),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "right-convergence",
            json.dumps(
                asdict(comparison.right_convergence),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "parameter-differences",
            json.dumps(
                [asdict(row) for row in comparison.parameter_differences],
                indent=2,
                sort_keys=True,
            ),
        ),
        ("limitations", json.dumps(limitations, indent=2)),
    ]
    machine_manifest = {
        "report_kind": "bayesian-run-comparison",
        "title": title,
        "left_tree_set_path": str(left_tree_set_path),
        "right_tree_set_path": str(right_tree_set_path),
        "left_trace_path": str(left_trace_path),
        "right_trace_path": str(right_trace_path),
        "trace_kind": trace_kind,
        "warning_count": len(comparison.warnings),
        "limitations": limitations,
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=method_tier_summary_metrics(method_tier),
    )
    return BayesianRunComparisonReportBuildResult(
        output_path=out_path,
        report_kind="bayesian-run-comparison",
        title=title,
        left_tree_set_path=left_tree_set_path,
        right_tree_set_path=right_tree_set_path,
        left_trace_path=left_trace_path,
        right_trace_path=right_trace_path,
        trace_kind=trace_kind,
        warning_count=len(comparison.warnings),
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )
