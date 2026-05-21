from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports.service import (
    render_alignment_report,
    render_dataset_report,
    render_phylo_inputs_report,
    render_taxon_report,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_input_report_commands(report_subparsers: Any) -> None:
    report_alignment = report_subparsers.add_parser(
        "alignment", help="Render an alignment-only HTML diagnostic report."
    )
    report_alignment.add_argument("--alignment", required=True, type=Path)
    report_alignment.add_argument("--out", required=True, type=Path)
    report_alignment.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_alignment)

    report_dataset = report_subparsers.add_parser(
        "dataset", help="Render a tree plus table dataset HTML report."
    )
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--alignment", type=Path)
    report_dataset.add_argument("--tip-dates", type=Path)
    report_dataset.add_argument("--calibrations", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_dataset)

    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_phylo_inputs)

    report_taxonomy = report_subparsers.add_parser(
        "taxonomy", help="Render a reviewer-facing taxon audit HTML report."
    )
    report_taxonomy.add_argument("--tree", required=True, type=Path)
    report_taxonomy.add_argument("--synonym-table", type=Path)
    report_taxonomy.add_argument("--metadata", type=Path)
    report_taxonomy.add_argument("--traits", type=Path)
    report_taxonomy.add_argument("--alignment", type=Path)
    report_taxonomy.add_argument("--filtered-alignment", type=Path)
    report_taxonomy.add_argument("--inference-tree", type=Path)
    report_taxonomy.add_argument("--reported-taxa", type=Path)
    report_taxonomy.add_argument("--out", required=True, type=Path)
    report_taxonomy.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_taxonomy)


def run_input_report_command(args: Any) -> int | None:
    if args.report_command == "alignment":
        result = render_alignment_report(
            alignment_path=args.alignment, out_path=args.out
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.alignment],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.alignment],
                    outputs=outputs,
                    warnings=result.alignment_forensic.warnings,
                    metrics={
                        "sequence_count": result.alignment.sequence_count,
                        "alignment_length": result.alignment.alignment_length,
                        "quality_score": result.alignment_quality.quality_score,
                        "warning_count": len(result.alignment_forensic.warnings),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "dataset":
        result = render_dataset_report(
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
            out_path=args.out,
        )
        inputs = [args.tree, args.metadata]
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "linked_taxa": result.metadata_linkage.linked_taxa,
                        "readiness_decision": None
                        if result.dataset_audit is None
                        else result.dataset_audit.readiness_decision,
                        "excluded_taxa": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.exclusion_table.rows),
                        "blocked_analysis_count": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.blocked_analyses),
                        "risky_analysis_count": 0
                        if result.dataset_audit is None
                        else sum(
                            1
                            for row in result.dataset_audit.analysis_decisions
                            if row.decision == "risky"
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "phylo-inputs":
        result = render_phylo_inputs_report(
            tree_path=args.tree,
            alignment_path=args.alignment,
            out_path=args.out,
        )
        inputs = [args.tree, args.alignment]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "alignment_length": result.alignment.alignment_length,
                        "linked_taxa": result.alignment_linkage.linked_taxa,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "taxonomy":
        result = render_taxon_report(
            tree_path=args.tree,
            synonym_table_path=args.synonym_table,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            out_path=args.out,
        )
        inputs = [args.tree]
        if args.synonym_table is not None:
            inputs.append(args.synonym_table)
        if args.metadata is not None:
            inputs.append(args.metadata)
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.taxon_audit.warnings,
                    metrics={
                        "tree_tip_count": result.taxon_audit.tree_tip_count,
                        "status": result.taxon_audit.status,
                        "conflict_count": len(
                            result.taxon_audit.mapping_conflicts.rows
                        ),
                        "crosswalk_rows": 0
                        if result.taxon_crosswalk is None
                        else len(result.taxon_crosswalk.rows),
                        "excluded_taxa": 0
                        if result.taxon_exclusions is None
                        else len(result.taxon_exclusions.rows),
                        "loss_stage_count": 0
                        if result.taxon_workflow_loss is None
                        else len(result.taxon_workflow_loss.loss_stage_counts),
                        "unstable_taxa": 0
                        if result.taxon_stability is None
                        else len(result.taxon_stability.unstable_taxa),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    return None
