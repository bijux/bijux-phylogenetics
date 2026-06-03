from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.posterior_sets.comparison import (
    compare_ml_tree_to_bayesian_posterior,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report

from .contracts import BayesianMlComparisonReportBuildResult
from .report_policy import (
    method_tier_section,
    method_tier_summary_metrics,
    ml_vs_bayesian_limitations,
)


def render_ml_vs_bayesian_tree_report(
    *,
    ml_tree_path: Path,
    posterior_tree_path: Path,
    out_path: Path,
    burnin_fraction: float = 0.25,
) -> BayesianMlComparisonReportBuildResult:
    """Render a reviewer-facing comparison between one ML tree and one Bayesian MCC summary."""
    comparison = compare_ml_tree_to_bayesian_posterior(
        ml_tree_path,
        posterior_tree_path,
        burnin_fraction=burnin_fraction,
    )
    title = "Bijux ML Versus Bayesian Tree Report"
    method_tier = bayesian_report_method_tier("ml-vs-bayesian-tree")
    limitations = ml_vs_bayesian_limitations(comparison.warnings)
    sections = [
        method_tier_section(method_tier),
        (
            "ml-versus-bayesian-summary",
            json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
        ),
        (
            "topology-comparison",
            json.dumps(
                asdict(comparison.topology), default=str, indent=2, sort_keys=True
            ),
        ),
        (
            "branch-length-comparison",
            json.dumps(
                asdict(comparison.branch_lengths), default=str, indent=2, sort_keys=True
            ),
        ),
        ("limitations", json.dumps(limitations, indent=2)),
    ]
    machine_manifest = {
        "report_kind": "ml-vs-bayesian-tree",
        "title": title,
        "ml_tree_path": str(ml_tree_path),
        "posterior_tree_path": str(posterior_tree_path),
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
    return BayesianMlComparisonReportBuildResult(
        output_path=out_path,
        report_kind="ml-vs-bayesian-tree",
        title=title,
        ml_tree_path=ml_tree_path,
        posterior_tree_path=posterior_tree_path,
        warning_count=len(comparison.warnings),
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )
