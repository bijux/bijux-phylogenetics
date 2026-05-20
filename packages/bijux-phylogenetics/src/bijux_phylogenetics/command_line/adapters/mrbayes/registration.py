from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    render_bayesian_posterior_report,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .diagnostics import (
    add_mrbayes_diagnostic_commands,
    run_mrbayes_diagnostic_command,
)
from .execution import (
    add_mrbayes_execution_commands,
    run_mrbayes_execution_command,
)
from .parameter_review import (
    add_mrbayes_parameter_review_commands,
    run_mrbayes_parameter_review_command,
)
from .posterior_trees import (
    add_mrbayes_posterior_tree_commands,
    run_mrbayes_posterior_tree_command,
)


def add_mrbayes_adapter_commands(adapter_subparsers: Any) -> None:
    add_mrbayes_execution_commands(adapter_subparsers)
    add_mrbayes_posterior_tree_commands(adapter_subparsers)
    add_mrbayes_diagnostic_commands(adapter_subparsers)
    add_mrbayes_parameter_review_commands(adapter_subparsers)

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


def run_mrbayes_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("mrbayes-"):
        return None

    execution_result = run_mrbayes_execution_command(args)
    if execution_result is not None:
        return execution_result

    posterior_tree_result = run_mrbayes_posterior_tree_command(args)
    if posterior_tree_result is not None:
        return posterior_tree_result

    diagnostic_result = run_mrbayes_diagnostic_command(args)
    if diagnostic_result is not None:
        return diagnostic_result

    parameter_review_result = run_mrbayes_parameter_review_command(args)
    if parameter_review_result is not None:
        return parameter_review_result

    if args.adapter_command == "mrbayes-report":
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

    return None
