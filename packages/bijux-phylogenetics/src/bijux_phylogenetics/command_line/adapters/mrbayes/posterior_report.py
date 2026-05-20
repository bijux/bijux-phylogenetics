from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import render_bayesian_posterior_report
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_posterior_report_commands(adapter_subparsers: Any) -> None:
    adapter_mrbayes_report = adapter_subparsers.add_parser(
        "mrbayes-report",
        help="Render an HTML Bayesian posterior report from posterior trees and traces.",
    )
    adapter_mrbayes_report.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_report.add_argument("--traces", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_report.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_mrbayes_report.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_report)


def run_mrbayes_posterior_report_command(args: Any) -> int | None:
    if args.adapter_command != "mrbayes-report":
        return None

    report = render_bayesian_posterior_report(
        posterior_tree_path=args.posterior_trees,
        trace_path=args.traces,
        out_path=args.out,
        burnin_fraction=args.burnin_fraction,
        ess_threshold=args.ess_threshold,
        mean_shift_threshold=args.mean_shift_threshold,
    )
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=[args.posterior_trees, args.traces],
        outputs=[report.output_path],
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=[args.posterior_trees, args.traces],
            outputs=outputs,
            warnings=method_tier_warnings(report.method_tier),
            metrics={
                "kept_tree_count": report.kept_tree_count,
                "warning_count": report.warning_count
                + len(method_tier_warnings(report.method_tier)),
                **method_tier_metrics(report.method_tier),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
