from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    assess_beast_burnin_sensitivity,
    assess_beast_convergence,
    parse_beast_log,
    parse_beast_posterior_tree_samples,
    prepare_beast_time_tree_analysis,
    render_calibration_audit_report,
    run_beast_posterior_inference,
    subsample_beast_posterior_tree_set,
    summarize_beast_analysis_xml,
    summarize_beast_log,
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
    write_beast_burnin_sensitivity_slice_table,
    write_beast_log_summary_table,
    write_beast_posterior_tree_set,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
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
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees.tree_set import (
    compute_clade_frequency_table,
    compute_tree_distance_matrix,
    write_clade_frequency_table,
    write_tree_distance_matrix,
)
from bijux_phylogenetics.trees.uncertainty import (
    cluster_trees_by_topology,
    detect_unstable_clades,
    write_topology_cluster_table,
    write_unstable_clade_table,
)


def add_beast_adapter_commands(adapter_subparsers: Any) -> None:
    adapter_beast_prepare = adapter_subparsers.add_parser(
        "beast-prepare",
        help="Prepare a BEAST-style time-tree XML analysis from aligned sequences and dating inputs.",
    )
    adapter_beast_prepare.add_argument("input_path", type=Path)
    adapter_beast_prepare.add_argument("--out", required=True, type=Path)
    adapter_beast_prepare.add_argument("--tree", type=Path)
    adapter_beast_prepare.add_argument("--calibrations", type=Path)
    adapter_beast_prepare.add_argument("--tip-dates", type=Path)
    adapter_beast_prepare.add_argument("--clock-model", default="strict")
    adapter_beast_prepare.add_argument("--tree-prior", default="yule")
    adapter_beast_prepare.add_argument("--chain-length", type=int, default=1000000)
    adapter_beast_prepare.add_argument("--log-every", type=int, default=1000)
    adapter_beast_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_beast_prepare)

    adapter_beast_xml = adapter_subparsers.add_parser(
        "beast-xml",
        help="Summarize and validate one prepared BEAST analysis XML.",
    )
    adapter_beast_xml.add_argument("input_path", type=Path)
    adapter_beast_xml.add_argument(
        "--json", action="store_true", help="Emit the XML summary report as JSON."
    )
    _add_manifest_argument(adapter_beast_xml)

    adapter_beast_run = adapter_subparsers.add_parser(
        "beast-run",
        help="Run a prepared BEAST posterior inference workflow.",
    )
    adapter_beast_run.add_argument("input_path", type=Path)
    adapter_beast_run.add_argument("--executable", type=str)
    adapter_beast_run.add_argument("--threads", type=int, default=1)
    adapter_beast_run.add_argument("--seed", type=int, default=1)
    adapter_beast_run.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Keep any existing posterior outputs instead of passing the BEAST overwrite flag.",
    )
    adapter_beast_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_beast_run)
    _add_manifest_argument(adapter_beast_run)

    adapter_beast_calibrations = adapter_subparsers.add_parser(
        "beast-calibrations",
        help="Validate a fossil calibration table against a tree.",
    )
    adapter_beast_calibrations.add_argument("tree_path", type=Path)
    adapter_beast_calibrations.add_argument("calibration_path", type=Path)
    adapter_beast_calibrations.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_calibrations)

    adapter_beast_tip_dates = adapter_subparsers.add_parser(
        "beast-tip-dates",
        help="Validate tip-dating metadata against a tree and optional alignment.",
    )
    adapter_beast_tip_dates.add_argument("tree_path", type=Path)
    adapter_beast_tip_dates.add_argument("tip_dates_path", type=Path)
    adapter_beast_tip_dates.add_argument("--alignment", type=Path)
    adapter_beast_tip_dates.add_argument("--date-column", default="date")
    adapter_beast_tip_dates.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_tip_dates)

    adapter_beast_log = adapter_subparsers.add_parser(
        "beast-log",
        help="Parse a BEAST log file into a deterministic numeric trace table.",
    )
    adapter_beast_log.add_argument("input_path", type=Path)
    adapter_beast_log.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting summaries.",
    )
    adapter_beast_log.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the parsed log.",
    )
    adapter_beast_log.add_argument(
        "--json", action="store_true", help="Emit the parsed log report as JSON."
    )
    _add_manifest_argument(adapter_beast_log)

    adapter_beast_burnin = adapter_subparsers.add_parser(
        "beast-burnin-sensitivity",
        help="Compare BEAST posterior summaries across multiple burn-in fractions.",
    )
    adapter_beast_burnin.add_argument("posterior_trees", type=Path)
    adapter_beast_burnin.add_argument("--log", type=Path)
    adapter_beast_burnin.add_argument(
        "--burnin-fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_BURNIN_FRACTIONS),
    )
    adapter_beast_burnin.add_argument("--slice-out", type=Path)
    adapter_beast_burnin.add_argument("--parameter-out", type=Path)
    adapter_beast_burnin.add_argument("--clade-out", type=Path)
    adapter_beast_burnin.add_argument(
        "--json",
        action="store_true",
        help="Emit the burn-in sensitivity report as JSON.",
    )
    _add_manifest_argument(adapter_beast_burnin)

    adapter_beast_parameters = adapter_subparsers.add_parser(
        "beast-parameters",
        help="Summarize burn-in-aware posterior parameter diagnostics from a BEAST log.",
    )
    adapter_beast_parameters.add_argument("input_path", type=Path)
    adapter_beast_parameters.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting posterior summaries.",
    )
    adapter_beast_parameters.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the retained log samples.",
    )
    adapter_beast_parameters.add_argument(
        "--json", action="store_true", help="Emit the parameter diagnostics as JSON."
    )
    _add_manifest_argument(adapter_beast_parameters)

    adapter_beast_trees = adapter_subparsers.add_parser(
        "beast-trees",
        help="Parse a BEAST posterior tree set into state-tagged normalized trees.",
    )
    adapter_beast_trees.add_argument("input_path", type=Path)
    adapter_beast_trees.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before reporting summaries.",
    )
    adapter_beast_trees.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_beast_trees)

    adapter_beast_subsample = adapter_subparsers.add_parser(
        "beast-subsample",
        help="Subsample BEAST posterior trees while preserving state metadata.",
    )
    adapter_beast_subsample.add_argument("input_path", type=Path)
    adapter_beast_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_beast_subsample.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before subsampling.",
    )
    adapter_beast_subsample.add_argument("--thinning-interval", type=int)
    adapter_beast_subsample.add_argument("--sample-count", type=int)
    adapter_beast_subsample.add_argument("--seed", type=int)
    adapter_beast_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_beast_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_beast_subsample)

    adapter_beast_consensus = adapter_subparsers.add_parser(
        "beast-consensus",
        help="Build a majority-rule consensus tree from BEAST posterior tree samples.",
    )
    adapter_beast_consensus.add_argument("input_path", type=Path)
    adapter_beast_consensus.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before consensus building.",
    )
    adapter_beast_consensus.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Write the posterior-probability-annotated consensus tree as Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--clade-table-out",
        type=Path,
        help="Write the retained clade-frequency ledger as TSV.",
    )
    adapter_beast_consensus.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior consensus report as JSON.",
    )
    _add_manifest_argument(adapter_beast_consensus)

    adapter_beast_diversity = adapter_subparsers.add_parser(
        "beast-diversity",
        help="Summarize topology diversity across BEAST posterior tree samples.",
    )
    adapter_beast_diversity.add_argument("input_path", type=Path)
    adapter_beast_diversity.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before topology review.",
    )
    adapter_beast_diversity.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_diversity.add_argument(
        "--distance-out",
        type=Path,
        help="Write the pairwise RF distance ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--topology-out",
        type=Path,
        help="Write the rooted topology cluster ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--unstable-clade-out",
        type=Path,
        help="Write the unstable-clade ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--json",
        action="store_true",
        help="Emit the topology diversity report as JSON.",
    )
    _add_manifest_argument(adapter_beast_diversity)

    adapter_beast_convergence = adapter_subparsers.add_parser(
        "beast-convergence",
        help="Assess BEAST log convergence from ESS and trace drift.",
    )
    adapter_beast_convergence.add_argument("input_path", type=Path)
    adapter_beast_convergence.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before assessing convergence.",
    )
    adapter_beast_convergence.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_beast_convergence.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_beast_convergence.add_argument(
        "--json", action="store_true", help="Emit the convergence report as JSON."
    )
    _add_manifest_argument(adapter_beast_convergence)

    adapter_beast_calibration_report = adapter_subparsers.add_parser(
        "beast-calibration-report",
        help="Render an HTML calibration audit report.",
    )
    adapter_beast_calibration_report.add_argument("tree_path", type=Path)
    adapter_beast_calibration_report.add_argument("calibration_path", type=Path)
    adapter_beast_calibration_report.add_argument("--out", required=True, type=Path)
    adapter_beast_calibration_report.add_argument("--tip-dates", type=Path)
    adapter_beast_calibration_report.add_argument("--alignment", type=Path)
    adapter_beast_calibration_report.add_argument("--date-column", default="date")
    adapter_beast_calibration_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_beast_calibration_report)


def run_beast_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("beast-"):
        return None

    if args.adapter_command == "beast-prepare":
        report = prepare_beast_time_tree_analysis(
            args.input_path,
            args.out,
            tree_path=args.tree,
            calibration_path=args.calibrations,
            tip_dates_path=args.tip_dates,
            clock_model=args.clock_model,
            tree_prior=args.tree_prior,
            chain_length=args.chain_length,
            log_every=args.log_every,
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
                    "calibration_count": report.calibration_count,
                    "tip_date_count": report.tip_date_count,
                    "warning_count": report.warning_count,
                    "starting_tree_source": report.starting_tree_source,
                    "beast_data_type": report.beast_data_type,
                    "substitution_model": report.substitution_model,
                    "clock_model": report.clock_model,
                    "tree_prior": report.tree_prior,
                    "chain_length": report.chain_length,
                    "log_every": report.log_every,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-xml":
        report = summarize_beast_analysis_xml(args.input_path)
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "valid": report.valid,
                    "issue_count": len(report.issues),
                    "taxon_count": report.taxon_count,
                    "character_count": report.character_count,
                    "calibration_count": report.calibration_count,
                    "tip_date_count": report.tip_date_count,
                    "chain_length": report.chain_length,
                    "logger_count": report.logger_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-run":
        report = run_beast_posterior_inference(
            args.input_path,
            executable=args.executable or "beast",
            overwrite=not args.no_overwrite,
            threads=args.threads,
            seed=args.seed,
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
                    "threads": args.threads,
                    "seed": args.seed,
                    "overwrite": not args.no_overwrite,
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-calibrations":
        report = validate_fossil_calibration_table(
            args.tree_path, args.calibration_path
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.tree_path, args.calibration_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.tree_path, args.calibration_path],
                outputs=outputs,
                metrics={
                    "calibration_count": report.calibration_count,
                    "invalid_calibration_count": report.invalid_calibration_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-tip-dates":
        report = validate_tip_dating_metadata(
            args.tree_path,
            args.tip_dates_path,
            alignment_path=args.alignment,
            date_column=args.date_column,
        )
        inputs = [
            args.tree_path,
            args.tip_dates_path,
            *([args.alignment] if args.alignment is not None else []),
        ]
        outputs = _finalize_outputs(args, command="adapter", inputs=inputs)
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "valid_tip_count": report.valid_tip_count,
                    "invalid_tip_count": report.invalid_tip_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-log":
        report = parse_beast_log(args.input_path)
        summary = summarize_beast_log(
            args.input_path, burnin_fraction=args.burnin_fraction
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(write_beast_log_summary_table(args.summary_out, summary))
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
                    "row_count": report.row_count,
                    "column_count": len(report.columns),
                    "burnin_fraction": summary.burnin_fraction,
                    "kept_row_count": summary.kept_row_count,
                    "posterior_parameter_count": len(summary.posterior_parameters),
                    "likelihood_parameter_count": len(summary.likelihood_parameters),
                    "prior_parameter_count": len(summary.prior_parameters),
                    "clock_parameter_count": len(summary.clock_parameters),
                    "tree_parameter_count": len(summary.tree_parameters),
                },
                data={"log": report, "summary": summary},
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-burnin-sensitivity":
        report = assess_beast_burnin_sensitivity(
            args.posterior_trees,
            log_path=args.log,
            burnin_fractions=tuple(args.burnin_fractions),
        )
        inputs = [args.posterior_trees, *([args.log] if args.log is not None else [])]
        outputs: list[Path | str] = []
        if args.slice_out is not None:
            outputs.append(write_beast_burnin_sensitivity_slice_table(args.slice_out, report))
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

    if args.adapter_command == "beast-parameters":
        report = summarize_beast_log(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(write_beast_log_summary_table(args.summary_out, report))
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
                    "posterior_parameter_count": len(report.posterior_parameters),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-trees":
        report = parse_beast_posterior_tree_samples(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.tree_set_out is not None:
            outputs.append(write_beast_posterior_tree_set(args.tree_set_out, report))
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
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_tree_count": report.rooted_tree_count,
                    "burnin_fraction": report.burnin_fraction,
                    "clade_count": len(report.clades),
                    "sampled_state_count": len(report.sampled_states),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-subsample":
        report = subsample_beast_posterior_tree_set(
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
                    "retained_state_count": len(
                        [tree for tree in report.trees if tree.state is not None]
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-consensus":
        consensus_tree, report = summarize_beast_posterior_trees(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = [write_newick(args.out, consensus_tree)]
        if args.tree_set_out is not None:
            args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
            args.tree_set_out.write_text(
                report.retained_tree_set_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            outputs.append(args.tree_set_out)
        if args.clade_table_out is not None:
            outputs.append(
                write_clade_frequency_table(
                    args.clade_table_out,
                    compute_clade_frequency_table(report.retained_tree_set_path),
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
                    "kept_tree_count": report.kept_tree_count,
                    "annotated_node_count": report.annotated_node_count,
                    "clade_frequency_count": report.clade_frequency_count,
                    "burnin_fraction": report.burnin_fraction,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-diversity":
        report = summarize_beast_posterior_topology_diversity(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.tree_set_out is not None:
            args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
            args.tree_set_out.write_text(
                Path(report.retained_tree_set_path).read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            outputs.append(args.tree_set_out)
        if args.distance_out is not None:
            outputs.append(
                write_tree_distance_matrix(
                    args.distance_out,
                    compute_tree_distance_matrix(report.retained_tree_set_path),
                )
            )
        if args.topology_out is not None:
            outputs.append(
                write_topology_cluster_table(
                    args.topology_out,
                    cluster_trees_by_topology(report.retained_tree_set_path),
                )
            )
        if args.unstable_clade_out is not None:
            outputs.append(
                write_unstable_clade_table(
                    args.unstable_clade_out,
                    detect_unstable_clades(report.retained_tree_set_path),
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
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_topology_count": report.rooted_topology_count,
                    "dominant_topology_frequency": report.dominant_topology_frequency,
                    "pair_count": report.pair_count,
                    "unstable_clade_count": report.unstable_clade_count,
                    "burnin_fraction": report.burnin_fraction,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-convergence":
        report = assess_beast_convergence(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
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
                    "burnin_fraction": report.burnin_fraction,
                    "sample_count": report.sample_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-calibration-report":
        report = render_calibration_audit_report(
            tree_path=args.tree_path,
            calibration_path=args.calibration_path,
            out_path=args.out,
            tip_dates_path=args.tip_dates,
            alignment_path=args.alignment,
            date_column=args.date_column,
        )
        inputs = [
            args.tree_path,
            args.calibration_path,
            *([args.tip_dates] if args.tip_dates is not None else []),
            *([args.alignment] if args.alignment is not None else []),
        ]
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=inputs,
            outputs=[report.output_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                warnings=method_tier_warnings(report.method_tier),
                metrics={
                    "invalid_calibration_count": report.invalid_calibration_count,
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
