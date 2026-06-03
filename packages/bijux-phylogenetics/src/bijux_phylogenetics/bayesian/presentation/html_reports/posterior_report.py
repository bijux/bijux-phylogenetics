from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.mrbayes import (
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_parameter_traces,
    summarize_mrbayes_posterior_trees,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.trees import compute_clade_frequency_table

from .contracts import BayesianPosteriorReportBuildResult
from .report_policy import (
    method_tier_section,
    method_tier_summary_metrics,
    posterior_report_limitations,
)


def render_bayesian_posterior_report(
    *,
    posterior_tree_path: Path,
    trace_path: Path,
    out_path: Path,
    burnin_fraction: float = 0.25,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BayesianPosteriorReportBuildResult:
    """Render a deterministic HTML report for one posterior tree and trace pair."""
    _, posterior_summary = summarize_mrbayes_posterior_trees(
        posterior_tree_path,
        burnin_fraction=burnin_fraction,
    )
    trace_report = parse_mrbayes_parameter_traces(trace_path)
    ess_report = compute_mrbayes_effective_sample_sizes(trace_path)
    convergence = assess_mrbayes_convergence(
        trace_path,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    clade_frequencies = compute_clade_frequency_table(
        posterior_summary.filtered_tree_set_path
    )
    title = "Bijux Bayesian Posterior Report"
    method_tier = bayesian_report_method_tier("bayesian-posterior")
    limitations = posterior_report_limitations(convergence.warnings)
    sections = [
        method_tier_section(method_tier),
        (
            "posterior-summary",
            json.dumps(
                asdict(posterior_summary), default=str, indent=2, sort_keys=True
            ),
        ),
        (
            "trace-summary",
            json.dumps(asdict(trace_report), default=str, indent=2, sort_keys=True),
        ),
        (
            "effective-sample-sizes",
            json.dumps(asdict(ess_report), default=str, indent=2, sort_keys=True),
        ),
        (
            "convergence",
            json.dumps(asdict(convergence), default=str, indent=2, sort_keys=True),
        ),
        (
            "clade-frequencies",
            json.dumps(
                asdict(clade_frequencies), default=str, indent=2, sort_keys=True
            ),
        ),
        ("limitations", json.dumps(limitations, indent=2)),
    ]
    machine_manifest = {
        "report_kind": "bayesian-posterior",
        "title": title,
        "posterior_tree_path": str(posterior_tree_path),
        "trace_path": str(trace_path),
        "filtered_tree_set_path": str(posterior_summary.filtered_tree_set_path),
        "kept_tree_count": posterior_summary.kept_tree_count,
        "warning_count": len(convergence.warnings),
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
    return BayesianPosteriorReportBuildResult(
        output_path=out_path,
        report_kind="bayesian-posterior",
        title=title,
        posterior_tree_path=posterior_tree_path,
        trace_path=trace_path,
        kept_tree_count=posterior_summary.kept_tree_count,
        warning_count=len(convergence.warnings),
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )
