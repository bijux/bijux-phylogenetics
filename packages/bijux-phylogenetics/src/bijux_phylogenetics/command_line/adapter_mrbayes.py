from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    parse_mrbayes_posterior_tree_samples,
    prepare_mrbayes_analysis,
    render_bayesian_posterior_report,
    run_mrbayes_posterior_inference,
    subsample_mrbayes_posterior_tree_set,
    summarize_mrbayes_parameter_diagnostics,
    summarize_mrbayes_posterior_trees,
    write_mrbayes_burnin_sensitivity_slice_table,
    write_mrbayes_parameter_summary_table,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.bayesian.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    write_burnin_clade_shift_table,
    write_burnin_parameter_shift_table,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_adapter_commands(adapter_subparsers: Any) -> None:
    adapter_mrbayes_prepare = adapter_subparsers.add_parser(
        "mrbayes-prepare",
        help="Prepare a MrBayes NEXUS analysis from an aligned FASTA file.",
    )
    adapter_mrbayes_prepare.add_argument("input_path", type=Path)
    adapter_mrbayes_prepare.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_prepare.add_argument("--partitions", type=Path)
    adapter_mrbayes_prepare.add_argument("--model", default="gtr")
    adapter_mrbayes_prepare.add_argument("--rates", default="gamma")
    adapter_mrbayes_prepare.add_argument("--ngen", type=int, default=10000)
    adapter_mrbayes_prepare.add_argument("--nchains", type=int, default=4)
    adapter_mrbayes_prepare.add_argument("--samplefreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--printfreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_prepare)

    adapter_mrbayes_run = adapter_subparsers.add_parser(
        "mrbayes-run",
        help="Run a prepared MrBayes posterior inference workflow.",
    )
    adapter_mrbayes_run.add_argument("input_path", type=Path)
    adapter_mrbayes_run.add_argument("--executable", type=str)
    adapter_mrbayes_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_mrbayes_run)
    _add_manifest_argument(adapter_mrbayes_run)

    adapter_mrbayes_summarize = adapter_subparsers.add_parser(
        "mrbayes-summarize",
        help="Summarize MrBayes posterior trees after burn-in removal.",
    )
    adapter_mrbayes_summarize.add_argument("input_path", type=Path)
    adapter_mrbayes_summarize.add_argument(
        "--burnin-fraction", type=float, default=0.25
    )
    adapter_mrbayes_summarize.add_argument(
        "--json", action="store_true", help="Emit the posterior summary as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_summarize)

    adapter_mrbayes_traces = adapter_subparsers.add_parser(
        "mrbayes-traces",
        help="Parse a MrBayes parameter trace table.",
    )
    adapter_mrbayes_traces.add_argument("input_path", type=Path)
    adapter_mrbayes_traces.add_argument(
        "--json", action="store_true", help="Emit the trace report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_traces)

    adapter_mrbayes_trees = adapter_subparsers.add_parser(
        "mrbayes-trees",
        help="Parse a MrBayes posterior tree set into sampled trees.",
    )
    adapter_mrbayes_trees.add_argument("input_path", type=Path)
    adapter_mrbayes_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_trees)

    adapter_mrbayes_subsample = adapter_subparsers.add_parser(
        "mrbayes-subsample",
        help="Subsample MrBayes posterior trees while preserving generation metadata.",
    )
    adapter_mrbayes_subsample.add_argument("input_path", type=Path)
    adapter_mrbayes_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_mrbayes_subsample.add_argument("--burnin-fraction", type=float, default=0.0)
    adapter_mrbayes_subsample.add_argument("--thinning-interval", type=int)
    adapter_mrbayes_subsample.add_argument("--sample-count", type=int)
    adapter_mrbayes_subsample.add_argument("--seed", type=int)
    adapter_mrbayes_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_subsample)

    adapter_mrbayes_mcmc = adapter_subparsers.add_parser(
        "mrbayes-mcmc",
        help="Parse a MrBayes MCMC diagnostics table.",
    )
    adapter_mrbayes_mcmc.add_argument("input_path", type=Path)
    adapter_mrbayes_mcmc.add_argument(
        "--json", action="store_true", help="Emit the MCMC diagnostics report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_mcmc)

    adapter_mrbayes_consensus = adapter_subparsers.add_parser(
        "mrbayes-consensus",
        help="Parse a MrBayes consensus tree with posterior-probability annotations.",
    )
    adapter_mrbayes_consensus.add_argument("input_path", type=Path)
    adapter_mrbayes_consensus.add_argument(
        "--json", action="store_true", help="Emit the consensus tree report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_consensus)

    adapter_mrbayes_ess = adapter_subparsers.add_parser(
        "mrbayes-ess",
        help="Compute effective sample sizes from a MrBayes trace table.",
    )
    adapter_mrbayes_ess.add_argument("input_path", type=Path)
    adapter_mrbayes_ess.add_argument(
        "--json", action="store_true", help="Emit the ESS report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_ess)

    adapter_mrbayes_parameters = adapter_subparsers.add_parser(
        "mrbayes-parameters",
        help="Summarize burn-in-aware posterior parameter diagnostics from a MrBayes trace table.",
    )
    adapter_mrbayes_parameters.add_argument("input_path", type=Path)
    adapter_mrbayes_parameters.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting posterior summaries.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the retained trace samples.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--json", action="store_true", help="Emit the parameter diagnostics as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_parameters)

    adapter_mrbayes_burnin = adapter_subparsers.add_parser(
        "mrbayes-burnin-sensitivity",
        help="Compare MrBayes posterior summaries across multiple burn-in fractions.",
    )
    adapter_mrbayes_burnin.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_burnin.add_argument("--traces", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--burnin-fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_BURNIN_FRACTIONS),
    )
    adapter_mrbayes_burnin.add_argument("--slice-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--parameter-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--clade-out", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--json",
        action="store_true",
        help="Emit the burn-in sensitivity report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_burnin)

    adapter_mrbayes_convergence = adapter_subparsers.add_parser(
        "mrbayes-convergence",
        help="Assess MrBayes trace convergence from ESS and trace drift.",
    )
    adapter_mrbayes_convergence.add_argument("input_path", type=Path)
    adapter_mrbayes_convergence.add_argument(
        "--ess-threshold", type=float, default=200.0
    )
    adapter_mrbayes_convergence.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_convergence.add_argument(
        "--json", action="store_true", help="Emit the convergence report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_convergence)

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

    if args.adapter_command == "mrbayes-prepare":
        report = prepare_mrbayes_analysis(
            args.input_path,
            args.out,
            partition_path=args.partitions,
            model=args.model,
            rates=args.rates,
            ngen=args.ngen,
            nchains=args.nchains,
            samplefreq=args.samplefreq,
            printfreq=args.printfreq,
            burnin_fraction=args.burnin_fraction,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[args.out],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "character_count": report.character_count,
                    "partitioned": report.partition_path is not None,
                    "partition_count": report.partition_count,
                    "partition_warning_count": len(report.partition_warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-run":
        report = run_mrbayes_posterior_inference(
            args.input_path,
            executable=args.executable or "mb",
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
            incomplete_run_policy=args.incomplete_run_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[*report.output_paths.values(), report.manifest_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=report.run.warning_lines,
                metrics={
                    "warning_count": len(report.run.warning_lines),
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-summarize":
        consensus_tree, report = summarize_mrbayes_posterior_trees(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[report.filtered_tree_set_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_topology_count": report.rooted_topology_count,
                    "tip_count": consensus_tree.tip_count,
                },
                data={
                    "summary": report,
                    "consensus_newick": report.consensus_newick,
                },
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-traces":
        report = parse_mrbayes_parameter_traces(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "column_count": len(report.columns),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-trees":
        report = parse_mrbayes_posterior_tree_samples(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "rooted_tree_count": report.rooted_tree_count,
                    "sampled_generation_count": len(report.sampled_generations),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-subsample":
        report = subsample_mrbayes_posterior_tree_set(
            args.input_path,
            method=args.method,
            thinning_interval=args.thinning_interval,
            sample_count=args.sample_count,
            burnin_fraction=args.burnin_fraction,
            random_seed=args.seed,
        )
        outputs: list[Path | str] = []
        if args.tree_set_out is not None:
            outputs.append(write_posterior_tree_subsample(args.tree_set_out, report))
        if args.sample_table_out is not None:
            outputs.append(
                write_posterior_tree_subsample_table(
                    args.sample_table_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "total_tree_count": report.total_tree_count,
                    "burnin_tree_count": report.burnin_tree_count,
                    "pre_subsampling_tree_count": report.pre_subsampling_tree_count,
                    "retained_tree_count": report.retained_tree_count,
                    "selection_method": report.selection_method,
                    "retained_generation_count": len(
                        [tree for tree in report.trees if tree.generation is not None]
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-mcmc":
        report = parse_mrbayes_mcmc_diagnostics(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "column_count": len(report.columns),
                    "comment_count": len(report.comment_lines),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-consensus":
        tree, report = parse_mrbayes_consensus_tree(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "tip_count": tree.tip_count,
                    "annotated_node_count": report.annotated_node_count,
                    "maximum_posterior_probability": (
                        report.maximum_posterior_probability
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-ess":
        report = compute_mrbayes_effective_sample_sizes(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={"parameter_count": len(report.effective_sample_sizes)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-parameters":
        report = summarize_mrbayes_parameter_diagnostics(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                write_mrbayes_parameter_summary_table(
                    args.summary_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "burnin_fraction": report.burnin_fraction,
                    "kept_row_count": report.kept_row_count,
                    "parameter_count": len(report.parameter_summaries),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-burnin-sensitivity":
        report = assess_mrbayes_burnin_sensitivity(
            args.posterior_trees,
            trace_path=args.traces,
            burnin_fractions=tuple(args.burnin_fractions),
        )
        inputs = [args.posterior_trees, *([args.traces] if args.traces is not None else [])]
        outputs: list[Path | str] = []
        if args.slice_out is not None:
            outputs.append(
                write_mrbayes_burnin_sensitivity_slice_table(
                    args.slice_out,
                    report,
                )
            )
        if args.parameter_out is not None:
            outputs.append(
                write_burnin_parameter_shift_table(
                    args.parameter_out,
                    report.parameter_shifts,
                )
            )
        if args.clade_out is not None:
            outputs.append(
                write_burnin_clade_shift_table(
                    args.clade_out,
                    report.clade_shifts,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=inputs,
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "slice_count": len(report.slices),
                    "parameter_shift_count": len(report.parameter_shifts),
                    "unstable_parameter_count": report.unstable_parameter_count,
                    "clade_shift_count": len(report.clade_shifts),
                    "unstable_clade_count": report.unstable_clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-convergence":
        report = assess_mrbayes_convergence(
            args.input_path,
            ess_threshold=args.ess_threshold,
            mean_shift_threshold=args.mean_shift_threshold,
        )
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=[warning["message"] for warning in report.warnings],
                metrics={
                    "warning_count": len(report.warnings),
                    "converged": report.converged,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

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
