from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_tree_inspection_commands(subparsers: Any) -> None:
    validate = subparsers.add_parser(
        get_command_spec("validate").name,
        help=get_command_spec("validate").summary,
    )
    validate.add_argument("tree", type=Path)
    validate.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    validate.add_argument("--allow-duplicates", action="store_true")
    validate.add_argument("--allow-negative-branches", action="store_true")
    validate.add_argument("--require-rooted", action="store_true")
    validate.add_argument("--require-ultrametric", action="store_true")
    validate.add_argument("--strict", action="store_true")
    validate.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(validate)

    inspect = subparsers.add_parser(
        get_command_spec("inspect").name,
        help=get_command_spec("inspect").summary,
    )
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(inspect)


def run_validate_command(args: Any) -> int:
    report = validate_tree_path(
        args.tree,
        source_format=args.format,
        allow_duplicates=args.allow_duplicates,
        strict=args.strict,
        allow_negative_branch_lengths=args.allow_negative_branches,
        require_rooted=args.require_rooted,
        require_ultrametric=args.require_ultrametric,
    )
    outputs = _finalize_outputs(args, command="validate", inputs=[args.tree])
    _print_result(
        build_command_result(
            command="validate",
            inputs=[args.tree],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "tip_count": report.tip_count,
                "internal_node_count": report.internal_node_count,
                "validity_decision": report.validity_decision,
                "syntax_valid": report.syntax_valid,
                "biologically_safe": report.biologically_safe,
                "polytomy_count": report.polytomy_count,
                "missing_internal_branch_count": len(
                    report.missing_internal_branch_nodes
                ),
                "missing_terminal_branch_count": len(
                    report.missing_terminal_branch_taxa
                ),
                "singleton_internal_node_count": len(report.singleton_internal_nodes),
                "integrity_issue_count": len(report.integrity_issues),
                "unsafe_external_label_count": len(report.unsafe_external_labels),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def run_inspect_command(args: Any) -> int:
    report = inspect_tree_path(args.tree, source_format=args.format)
    outputs = _finalize_outputs(args, command="inspect", inputs=[args.tree])
    _print_result(
        build_command_result(
            command="inspect",
            inputs=[args.tree],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "tip_count": report.tip_count,
                "node_count": report.node_count,
                "internal_node_count": report.internal_node_count,
                "edge_count": report.edge_count,
                "clade_count": report.clade_count,
                "is_binary": report.is_binary,
                "polytomy_count": report.polytomy_count,
                "branch_length_status": report.branch_length_status,
                "is_ultrametric": report.is_ultrametric,
                "tree_diameter": report.tree_diameter,
                "colless_imbalance_index": report.colless_imbalance_index,
                "sackin_imbalance_index": report.sackin_imbalance_index,
                "tree_quality_score": report.tree_quality_score,
                "zero_length_branch_count": report.zero_length_branch_count,
                "cherry_count": report.cherry_count,
                "missing_internal_branch_count": len(
                    report.missing_internal_branch_nodes
                ),
                "missing_terminal_branch_count": len(
                    report.missing_terminal_branch_taxa
                ),
                "singleton_internal_node_count": len(report.singleton_internal_nodes),
                "long_branch_outlier_count": len(report.long_branch_outliers),
                "short_branch_outlier_count": len(report.short_branch_outliers),
                "likely_support_label_count": len(report.likely_support_labels),
                "likely_named_internal_label_count": len(
                    report.likely_named_internal_labels
                ),
                "suspicious_support_range_count": len(
                    report.suspicious_support_value_ranges
                ),
                "root_classification": report.root_state_confidence.classification,
                "internal_label_conflict_count": len(report.internal_label_conflicts),
                "unsafe_external_label_count": len(report.unsafe_external_labels),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
