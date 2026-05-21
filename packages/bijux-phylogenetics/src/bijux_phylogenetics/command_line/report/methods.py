from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
    write_ancestral_methods_summary_text,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import (
    write_alignment_filtering_methods_summary_text,
    write_tree_inference_methods_summary_text,
    write_tree_validation_methods_summary_text,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_methods_report_commands(report_subparsers: Any) -> None:
    report_tree_validation_methods_summary = report_subparsers.add_parser(
        "tree-validation-methods-summary",
        help="Write reviewer-facing methods-summary text for tree validation.",
    )
    report_tree_validation_methods_summary.add_argument("tree", type=Path)
    report_tree_validation_methods_summary.add_argument("--source-format")
    report_tree_validation_methods_summary.add_argument(
        "--out", required=True, type=Path
    )
    report_tree_validation_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(report_tree_validation_methods_summary)

    report_tree_inference_methods_summary = report_subparsers.add_parser(
        "tree-inference-methods-summary",
        help="Write reviewer-facing methods-summary text for a fasta-to-tree workflow manifest.",
    )
    report_tree_inference_methods_summary.add_argument("workflow_manifest", type=Path)
    report_tree_inference_methods_summary.add_argument(
        "--out", required=True, type=Path
    )
    report_tree_inference_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(report_tree_inference_methods_summary)

    report_ancestral_methods_summary = report_subparsers.add_parser(
        "ancestral-methods-summary",
        help="Write reviewer-facing methods-summary text for an ancestral reconstruction.",
    )
    report_ancestral_methods_summary.add_argument("tree", type=Path)
    report_ancestral_methods_summary.add_argument("traits", type=Path)
    report_ancestral_methods_summary.add_argument(
        "--trait",
        required=True,
        help="Trait column used for ancestral reconstruction.",
    )
    report_ancestral_methods_summary.add_argument(
        "--kind",
        choices=("continuous", "discrete"),
        required=True,
    )
    report_ancestral_methods_summary.add_argument("--taxon-column")
    report_ancestral_methods_summary.add_argument(
        "--model",
        help="Continuous model such as brownian or ou, or discrete model such as equal-rates, symmetric, all-rates-different, or fitch.",
    )
    report_ancestral_methods_summary.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Positive alpha value for continuous OU-style reconstruction.",
    )
    report_ancestral_methods_summary.add_argument(
        "--state-ordering",
        choices=("unordered", "ordered"),
        default="unordered",
    )
    report_ancestral_methods_summary.add_argument("--ordered-states")
    report_ancestral_methods_summary.add_argument(
        "--root-prior-mode",
        choices=("equal", "empirical", "fixed"),
        default="equal",
    )
    report_ancestral_methods_summary.add_argument("--fixed-root-state")
    report_ancestral_methods_summary.add_argument("--out", required=True, type=Path)
    report_ancestral_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(report_ancestral_methods_summary)

    report_alignment_filtering_methods_summary = report_subparsers.add_parser(
        "alignment-filtering-methods-summary",
        help="Write reviewer-facing methods-summary text for profile-driven alignment filtering.",
    )
    report_alignment_filtering_methods_summary.add_argument("alignment", type=Path)
    report_alignment_filtering_methods_summary.add_argument("--profile", required=True)
    report_alignment_filtering_methods_summary.add_argument("--group-table", type=Path)
    report_alignment_filtering_methods_summary.add_argument(
        "--group-column",
        dest="group_columns",
        action="append",
        default=None,
    )
    report_alignment_filtering_methods_summary.add_argument(
        "--out", required=True, type=Path
    )
    report_alignment_filtering_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(report_alignment_filtering_methods_summary)


def run_methods_report_command(args: Any) -> int | None:
    if args.report_command == "tree-validation-methods-summary":
        result = write_tree_validation_methods_summary_text(
            args.out,
            tree_path=args.tree,
            source_format=args.source_format,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree],
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=result.forensic.warnings,
                    metrics={
                        "warning_count": result.warning_count,
                        "blocked_context_count": result.blocked_context_count,
                        "repair_item_count": result.repair_item_count,
                        "validity_decision": result.validation.validity_decision,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "tree-inference-methods-summary":
        result = write_tree_inference_methods_summary_text(
            args.out,
            workflow_manifest_path=args.workflow_manifest,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.workflow_manifest],
            outputs=[result.output_path],
        )
        workflow_warnings = [
            str(item)
            for item in result.workflow_manifest.get("warnings", [])
            if str(item).strip()
        ]
        support_summary = result.workflow_manifest.get("support_summary")
        if isinstance(support_summary, dict):
            workflow_warnings.extend(
                str(item)
                for item in support_summary.get("warnings", [])
                if str(item).strip()
            )
        deduplicated_warnings: list[str] = []
        for warning in workflow_warnings:
            if warning not in deduplicated_warnings:
                deduplicated_warnings.append(warning)
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.workflow_manifest],
                    outputs=outputs,
                    warnings=deduplicated_warnings,
                    metrics={
                        "warning_count": result.warning_count,
                        "selected_model": result.selected_model,
                        "bootstrap_replicates": result.bootstrap_replicates,
                        "trimmed_alignment_length": result.trimmed_alignment_length,
                        "supported_node_count": result.supported_node_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "ancestral-methods-summary":
        if args.kind == "continuous":
            reconstruction = reconstruct_continuous_ancestral_states(
                args.tree,
                args.traits,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model or "brownian",
                alpha=args.alpha,
            )
        else:
            reconstruction = reconstruct_discrete_ancestral_states(
                args.tree,
                args.traits,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model or "equal-rates",
                state_ordering=args.state_ordering,
                ordered_states=_split_csv_values(args.ordered_states) or None,
                root_prior_mode=args.root_prior_mode,
                fixed_root_state=args.fixed_root_state,
            )
        result = write_ancestral_methods_summary_text(
            args.out,
            reconstruction_kind=args.kind,
            reconstruction=reconstruction,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree, args.traits],
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree, args.traits],
                    outputs=outputs,
                    warnings=reconstruction.warnings,
                    metrics={
                        "reconstruction_kind": result.reconstruction_kind,
                        "model": result.model,
                        "analyzed_taxon_count": result.analyzed_taxon_count,
                        "unstable_node_count": result.unstable_node_count,
                        "warning_count": result.warning_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "alignment-filtering-methods-summary":
        result = write_alignment_filtering_methods_summary_text(
            args.out,
            alignment_path=args.alignment,
            profile_name=args.profile,
            group_table_path=args.group_table,
            group_columns=args.group_columns,
        )
        inputs = (
            [args.alignment, args.group_table]
            if args.group_table is not None
            else [args.alignment]
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        warnings = [
            *(warning.message for warning in result.cleaning.signal_warnings),
            *result.cleaning.warnings,
            *result.cleaning.comparison.warnings,
        ]
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=warnings,
                    metrics={
                        "warning_count": result.warning_count,
                        "removed_site_count": result.removed_site_count,
                        "removed_sequence_count": result.removed_sequence_count,
                        "retained_sequence_count": result.retained_sequence_count,
                        "retained_alignment_length": result.retained_alignment_length,
                        "profile_name": result.cleaning.profile.name,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    return None
