from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def _command_line_api() -> Any:
    import bijux_phylogenetics.command_line as command_line_api

    return command_line_api


def add_benchmark_commands(subparsers: Any) -> None:
    benchmark = subparsers.add_parser(
        get_command_spec("benchmark").name, help=get_command_spec("benchmark").summary
    )
    benchmark_subparsers = benchmark.add_subparsers(
        dest="benchmark_command", required=True
    )
    benchmark_validate = benchmark_subparsers.add_parser(
        "tree-validation",
        help="Benchmark tree validation across size classes.",
    )
    benchmark_validate.add_argument("--replicates", type=int, default=3)
    benchmark_validate.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_validate)

    benchmark_compare = benchmark_subparsers.add_parser(
        "tree-comparison",
        help="Benchmark tree comparison across increasing taxon counts.",
    )
    benchmark_compare.add_argument("--replicates", type=int, default=3)
    benchmark_compare.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_compare)

    benchmark_large_tree = benchmark_subparsers.add_parser(
        "large-tree-scaling",
        help="Benchmark large-tree validation, comparison, rendering, and reporting.",
    )
    benchmark_large_tree.add_argument("--replicates", type=int, default=1)
    benchmark_large_tree.add_argument(
        "--tip-count",
        action="append",
        dest="tip_counts",
        type=int,
        help="Add one governed tree size to benchmark. Repeat to benchmark multiple sizes.",
    )
    benchmark_large_tree.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree)

    benchmark_large_alignment = benchmark_subparsers.add_parser(
        "large-alignment-scaling",
        help="Benchmark large-alignment diagnostics, trimming, distance, and readiness.",
    )
    benchmark_large_alignment.add_argument("--replicates", type=int, default=1)
    benchmark_large_alignment.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_alignment.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one alignment length to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_alignment.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_alignment)

    benchmark_large_tree_set = benchmark_subparsers.add_parser(
        "large-tree-set-scaling",
        help="Benchmark large-tree-set consensus, RF diversity, clustering, and uncertainty summaries.",
    )
    benchmark_large_tree_set.add_argument("--replicates", type=int, default=1)
    benchmark_large_tree_set.add_argument(
        "--tree-count",
        action="append",
        dest="tree_counts",
        type=int,
        help="Add one posterior tree count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_tree_set.add_argument(
        "--tip-count",
        action="append",
        dest="tip_counts",
        type=int,
        help="Add one taxon count to benchmark. Repeat to benchmark multiple size classes.",
    )
    benchmark_large_tree_set.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree_set)

    benchmark_practical_limits = benchmark_subparsers.add_parser(
        "workflow-practical-limits",
        help="Report the largest governed workflow classes currently exercised in benchmark and stress lanes.",
    )
    benchmark_practical_limits.add_argument("--replicates", type=int, default=1)
    benchmark_practical_limits.add_argument(
        "--tree-tip-count",
        action="append",
        dest="tree_tip_counts",
        type=int,
        help="Add one large-tree taxon count. Repeat to override the governed tree-size classes.",
    )
    benchmark_practical_limits.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count for the large-alignment classes. Repeat alongside --alignment-length.",
    )
    benchmark_practical_limits.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one aligned-site count for the large-alignment classes. Repeat alongside --sequence-count.",
    )
    benchmark_practical_limits.add_argument(
        "--posterior-tree-count",
        action="append",
        dest="posterior_tree_counts",
        type=int,
        help="Add one posterior tree count for the tree-set classes. Repeat alongside --tree-set-tip-count.",
    )
    benchmark_practical_limits.add_argument(
        "--tree-set-tip-count",
        action="append",
        dest="tree_set_tip_counts",
        type=int,
        help="Add one taxon count for the tree-set classes. Repeat alongside --posterior-tree-count.",
    )
    benchmark_practical_limits.add_argument(
        "--stress-tier",
        action="append",
        dest="stress_tiers",
        choices=("small", "heavy"),
        help="Include one governed stress tier. Repeat to aggregate multiple tiers.",
    )
    benchmark_practical_limits.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_practical_limits)

    benchmark_alignment = benchmark_subparsers.add_parser(
        "alignment-diagnostics",
        help="Benchmark alignment diagnostics across increasing sequence counts.",
    )
    benchmark_alignment.add_argument("--replicates", type=int, default=3)
    benchmark_alignment.add_argument("--sequence-length", type=int, default=128)
    benchmark_alignment.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_alignment)

    benchmark_stress = benchmark_subparsers.add_parser(
        "stress-suite",
        help="Benchmark large-dataset stress workloads across governed tiers.",
    )
    benchmark_stress.add_argument(
        "--tier",
        choices=("small", "heavy"),
        default="small",
        help="Select the governed stress tier to execute.",
    )
    benchmark_stress.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_stress)

    benchmark_large_tree_model = benchmark_subparsers.add_parser(
        "large-tree-model-fitting",
        help="Benchmark 100+ taxon continuous and discrete model fitting with governed geiger comparison and heavy-tier review.",
    )
    benchmark_large_tree_model.add_argument(
        "--tier",
        choices=("small", "heavy"),
        default="small",
        help="Select the governed model-fitting tier to execute.",
    )
    benchmark_large_tree_model.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_large_tree_model)

    benchmark_real_dataset = benchmark_subparsers.add_parser(
        "real-dataset-macroevolution",
        help="Benchmark continuous and discrete macroevolution model fitting on the published Central European seashore flora dataset against stored local geiger references.",
    )
    benchmark_real_dataset.add_argument(
        "--json", action="store_true", help="Emit the benchmark report as JSON."
    )
    _add_manifest_argument(benchmark_real_dataset)


def run_benchmark_command(args: Any) -> int:
    command_line_api = _command_line_api()
    if args.benchmark_command == "tree-validation":
        report = command_line_api.benchmark_tree_validation(replicates=args.replicates)
    elif args.benchmark_command == "tree-comparison":
        report = command_line_api.benchmark_tree_comparison(replicates=args.replicates)
    elif args.benchmark_command == "large-tree-scaling":
        report = command_line_api.benchmark_large_tree_scaling(
            replicates=args.replicates,
            tip_counts=args.tip_counts,
        )
    elif args.benchmark_command == "large-alignment-scaling":
        classes = None
        if args.sequence_counts is not None or args.alignment_lengths is not None:
            sequence_counts = args.sequence_counts or []
            alignment_lengths = args.alignment_lengths or []
            if len(sequence_counts) != len(alignment_lengths):
                raise ValueError(
                    "large-alignment-scaling requires the same number of --sequence-count and --alignment-length values"
                )
            classes = [
                (
                    f"sequences-{sequence_count}-sites-{alignment_length}",
                    sequence_count,
                    alignment_length,
                )
                for sequence_count, alignment_length in zip(
                    sequence_counts,
                    alignment_lengths,
                    strict=True,
                )
            ]
        report = command_line_api.benchmark_large_alignment_scaling(
            replicates=args.replicates,
            size_classes=classes,
        )
    elif args.benchmark_command == "large-tree-set-scaling":
        classes = None
        if args.tree_counts is not None or args.tip_counts is not None:
            tree_counts = args.tree_counts or []
            tip_counts = args.tip_counts or []
            if len(tree_counts) != len(tip_counts):
                raise ValueError(
                    "large-tree-set-scaling requires the same number of --tree-count and --tip-count values"
                )
            classes = [
                (f"trees-{tree_count}-taxa-{tip_count}", tree_count, tip_count)
                for tree_count, tip_count in zip(
                    tree_counts,
                    tip_counts,
                    strict=True,
                )
            ]
        report = command_line_api.benchmark_large_tree_set_scaling(
            replicates=args.replicates,
            size_classes=classes,
        )
    elif args.benchmark_command == "workflow-practical-limits":
        alignment_classes = None
        if args.sequence_counts is not None or args.alignment_lengths is not None:
            sequence_counts = args.sequence_counts or []
            alignment_lengths = args.alignment_lengths or []
            if len(sequence_counts) != len(alignment_lengths):
                raise ValueError(
                    "workflow-practical-limits requires the same number of --sequence-count and --alignment-length values"
                )
            alignment_classes = [
                (
                    f"sequences-{sequence_count}-sites-{alignment_length}",
                    sequence_count,
                    alignment_length,
                )
                for sequence_count, alignment_length in zip(
                    sequence_counts,
                    alignment_lengths,
                    strict=True,
                )
            ]
        tree_set_classes = None
        if (
            args.posterior_tree_counts is not None
            or args.tree_set_tip_counts is not None
        ):
            posterior_tree_counts = args.posterior_tree_counts or []
            tree_set_tip_counts = args.tree_set_tip_counts or []
            if len(posterior_tree_counts) != len(tree_set_tip_counts):
                raise ValueError(
                    "workflow-practical-limits requires the same number of --posterior-tree-count and --tree-set-tip-count values"
                )
            tree_set_classes = [
                (
                    f"trees-{tree_count}-taxa-{tip_count}",
                    tree_count,
                    tip_count,
                )
                for tree_count, tip_count in zip(
                    posterior_tree_counts,
                    tree_set_tip_counts,
                    strict=True,
                )
            ]
        report = command_line_api.benchmark_workflow_practical_limits(
            replicates=args.replicates,
            tree_tip_counts=args.tree_tip_counts,
            alignment_size_classes=alignment_classes,
            tree_set_size_classes=tree_set_classes,
            stress_tiers=args.stress_tiers,
        )
    elif args.benchmark_command == "stress-suite":
        report = command_line_api.benchmark_large_dataset_stress_suite(tier=args.tier)
    elif args.benchmark_command == "large-tree-model-fitting":
        report = command_line_api.benchmark_large_tree_model_fitting(tier=args.tier)
    elif args.benchmark_command == "real-dataset-macroevolution":
        report = command_line_api.benchmark_real_dataset_macroevolution()
    else:
        report = command_line_api.benchmark_alignment_diagnostics(
            replicates=args.replicates,
            sequence_length=args.sequence_length,
        )

    outputs = _finalize_outputs(args, command="benchmark", inputs=[])
    if hasattr(report, "entries"):
        metrics = {
            "entry_count": len(report.entries),
        }
    elif hasattr(report, "summary_rows"):
        metrics = {
            "summary_row_count": len(report.summary_rows),
            "model_row_count": len(report.model_rows),
            "alignment_review_row_count": len(report.alignment_review_rows),
            "parity_row_count": len(report.parity_rows),
        }
    else:
        metrics = {
            "observation_count": (
                len(report.observations)
                if hasattr(report, "observations")
                else sum(len(row.observations) for row in report.workflows)
            ),
        }
    if hasattr(report, "replicates"):
        metrics["replicates"] = report.replicates
    if hasattr(report, "tier"):
        metrics["tier"] = report.tier
    if hasattr(report, "case_count"):
        metrics["case_count"] = report.case_count
    if hasattr(report, "geiger_match_case_count"):
        metrics["geiger_match_case_count"] = report.geiger_match_case_count
    if hasattr(report, "threshold_pass_case_count"):
        metrics["threshold_pass_case_count"] = report.threshold_pass_case_count
    if hasattr(report, "too_slow_case_count"):
        metrics["too_slow_case_count"] = report.too_slow_case_count
    if hasattr(report, "unstable_case_count"):
        metrics["unstable_case_count"] = report.unstable_case_count
    if hasattr(report, "stress_tiers"):
        metrics["stress_tier_count"] = len(report.stress_tiers)
    if hasattr(report, "workflows"):
        metrics["workflow_count"] = len(report.workflows)
        if hasattr(report, "tip_counts"):
            metrics["max_tip_count"] = max(report.tip_counts)
        if hasattr(report, "tree_counts"):
            metrics["max_tree_count"] = max(report.tree_counts)
        if hasattr(report, "alignment_lengths"):
            metrics["max_alignment_length"] = max(report.alignment_lengths)
            metrics["max_sequence_count"] = max(report.sequence_counts)
    if hasattr(report, "entries"):
        metrics["workflow_count"] = len(report.entries)
        taxon_limits = [
            row.tested_taxon_limit
            for row in report.entries
            if row.tested_taxon_limit is not None
        ]
        site_limits = [
            row.tested_site_limit
            for row in report.entries
            if row.tested_site_limit is not None
        ]
        tree_limits = [
            row.tested_tree_limit
            for row in report.entries
            if row.tested_tree_limit is not None
        ]
        posterior_limits = [
            row.tested_posterior_size
            for row in report.entries
            if row.tested_posterior_size is not None
        ]
        if taxon_limits:
            metrics["max_taxon_limit"] = max(taxon_limits)
        if site_limits:
            metrics["max_site_limit"] = max(site_limits)
        if tree_limits:
            metrics["max_tree_limit"] = max(tree_limits)
        if posterior_limits:
            metrics["max_posterior_size"] = max(posterior_limits)
    if hasattr(report, "observations") and report.observations:
        taxon_counts = [
            row.taxon_count
            for row in report.observations
            if hasattr(row, "taxon_count") and row.taxon_count is not None
        ]
        if taxon_counts:
            metrics["max_taxon_count"] = max(taxon_counts)
    _print_result(
        build_command_result(
            command="benchmark",
            inputs=[],
            outputs=outputs,
            metrics=metrics,
            data=report,
        ),
        json_output=args.json,
    )
    return 0
