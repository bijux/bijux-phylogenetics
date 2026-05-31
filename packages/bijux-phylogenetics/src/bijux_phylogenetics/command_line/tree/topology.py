from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.reference import (
    validate_support_reference_examples,
    validate_tree_distance_reference_examples,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.topology import (
    enumerate_rooted_nni_neighbors,
    enumerate_rooted_spr_neighbors,
    enumerate_rooted_tbr_neighbors,
    reroot_tree_by_midpoint,
    root_tree_on_outgroup,
    RootedSprEnumerationBudget,
    summarize_rooted_nni_move_application,
    summarize_rooted_spr_move_application,
    summarize_rooted_tbr_move_application,
    unroot_tree,
    write_rooted_nni_artifacts,
    write_rooted_nni_move_artifacts,
    write_rooted_spr_artifacts,
    write_rooted_spr_move_artifacts,
    write_rooted_tbr_artifacts,
    write_rooted_tbr_move_artifacts,
    write_tree_rooting_report,
)
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    extract_tree_clades,
    summarize_tree_shape,
    write_branch_length_table,
    write_clade_table,
    write_tree_shape_table,
)


def add_topology_commands(subparsers: Any) -> None:
    topology = subparsers.add_parser(
        get_command_spec("topology").name,
        help=get_command_spec("topology").summary,
    )
    topology_subparsers = topology.add_subparsers(
        dest="topology_command", required=True
    )

    topology_clades = topology_subparsers.add_parser(
        "clades",
        help="Extract one row per clade from a rooted tree.",
    )
    topology_clades.add_argument("tree", type=Path)
    topology_clades.add_argument("--metadata", type=Path)
    topology_clades.add_argument("--taxon-column", type=str)
    topology_clades.add_argument(
        "--metadata-column",
        dest="metadata_columns",
        action="append",
        default=[],
        help="Summarize this taxon-keyed metadata column for each clade.",
    )
    topology_clades.add_argument("--out", type=Path)
    topology_clades.add_argument(
        "--json", action="store_true", help="Emit the clade table report as JSON."
    )
    _add_manifest_argument(topology_clades)

    topology_shape = topology_subparsers.add_parser(
        "shape",
        help="Summarize rooted tree shape, balance, and ladderization.",
    )
    topology_shape.add_argument("tree", type=Path)
    topology_shape.add_argument(
        "--out", type=Path, help="Write one row summarizing the tree shape as TSV."
    )
    topology_shape.add_argument(
        "--json", action="store_true", help="Emit the tree-shape report as JSON."
    )
    _add_manifest_argument(topology_shape)

    topology_branch_lengths = topology_subparsers.add_parser(
        "branch-lengths",
        help="Summarize branch-length distribution and detect unusual edges.",
    )
    topology_branch_lengths.add_argument("tree", type=Path)
    topology_branch_lengths.add_argument(
        "--out", type=Path, help="Write one row per branch as TSV."
    )
    topology_branch_lengths.add_argument(
        "--json",
        action="store_true",
        help="Emit the branch-length distribution report as JSON.",
    )
    _add_manifest_argument(topology_branch_lengths)

    topology_distance_reference = topology_subparsers.add_parser(
        "distance-reference",
        help="Validate governed topology-distance reference examples.",
    )
    topology_distance_reference.add_argument(
        "--json", action="store_true", help="Emit the reference report as JSON."
    )
    _add_manifest_argument(topology_distance_reference)

    topology_support_reference = topology_subparsers.add_parser(
        "support-reference",
        help="Validate governed clade-support reference examples.",
    )
    topology_support_reference.add_argument(
        "--json", action="store_true", help="Emit the reference report as JSON."
    )
    _add_manifest_argument(topology_support_reference)

    topology_outgroup = topology_subparsers.add_parser(
        "root-outgroup", help="Root a tree on a named outgroup taxon set."
    )
    topology_outgroup.add_argument("tree", type=Path)
    topology_outgroup.add_argument("--taxa", nargs="+", required=True)
    topology_outgroup.add_argument("--out", required=True, type=Path)
    topology_outgroup.add_argument(
        "--report-out",
        type=Path,
        help="Write a rooting audit report as JSON.",
    )
    topology_outgroup.add_argument(
        "--json", action="store_true", help="Emit the rooting report as JSON."
    )
    _add_manifest_argument(topology_outgroup)

    topology_midpoint = topology_subparsers.add_parser(
        "reroot-midpoint", help="Reroot a tree at its midpoint."
    )
    topology_midpoint.add_argument("tree", type=Path)
    topology_midpoint.add_argument("--out", required=True, type=Path)
    topology_midpoint.add_argument(
        "--report-out",
        type=Path,
        help="Write a rooting audit report as JSON.",
    )
    topology_midpoint.add_argument(
        "--json", action="store_true", help="Emit the rerooting report as JSON."
    )
    _add_manifest_argument(topology_midpoint)

    topology_unroot = topology_subparsers.add_parser(
        "unroot", help="Remove rooting from a tree and write the unrooted topology."
    )
    topology_unroot.add_argument("tree", type=Path)
    topology_unroot.add_argument("--out", required=True, type=Path)
    topology_unroot.add_argument(
        "--json", action="store_true", help="Emit the unrooting report as JSON."
    )
    _add_manifest_argument(topology_unroot)

    topology_rooted_nni_neighbors = topology_subparsers.add_parser(
        "rooted-nni-neighbors",
        help="Enumerate every legal rooted NNI neighbor for one binary-root tree representation.",
    )
    topology_rooted_nni_neighbors.add_argument("tree", type=Path)
    topology_rooted_nni_neighbors.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_nni_neighbors.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted NNI neighborhood report as JSON.",
    )
    _add_manifest_argument(topology_rooted_nni_neighbors)

    topology_rooted_nni_apply = topology_subparsers.add_parser(
        "rooted-nni-apply",
        help="Apply one indexed rooted NNI move, derive its reverse, and write governed artifacts.",
    )
    topology_rooted_nni_apply.add_argument("tree", type=Path)
    topology_rooted_nni_apply.add_argument("--move-index", required=True, type=int)
    topology_rooted_nni_apply.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_nni_apply.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted NNI move-application report as JSON.",
    )
    _add_manifest_argument(topology_rooted_nni_apply)

    topology_rooted_spr_neighbors = topology_subparsers.add_parser(
        "rooted-spr-neighbors",
        help="Enumerate rooted SPR neighbors for one binary-root tree representation.",
    )
    topology_rooted_spr_neighbors.add_argument("tree", type=Path)
    topology_rooted_spr_neighbors.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_spr_neighbors.add_argument(
        "--max-pruned-clades",
        type=int,
        help="Limit rooted SPR enumeration to this many deterministically ordered pruned clades.",
    )
    topology_rooted_spr_neighbors.add_argument(
        "--max-regraft-targets-per-pruned-clade",
        type=int,
        help="Limit each retained pruned clade to this many deterministically ordered regraft targets.",
    )
    topology_rooted_spr_neighbors.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted SPR neighborhood report as JSON.",
    )
    _add_manifest_argument(topology_rooted_spr_neighbors)

    topology_rooted_spr_apply = topology_subparsers.add_parser(
        "rooted-spr-apply",
        help="Apply one indexed rooted SPR move and write governed artifacts.",
    )
    topology_rooted_spr_apply.add_argument("tree", type=Path)
    topology_rooted_spr_apply.add_argument("--move-index", required=True, type=int)
    topology_rooted_spr_apply.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_spr_apply.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted SPR move-application report as JSON.",
    )
    _add_manifest_argument(topology_rooted_spr_apply)

    topology_rooted_tbr_neighbors = topology_subparsers.add_parser(
        "rooted-tbr-neighbors",
        help="Enumerate rooted TBR neighbors for one binary-root tree representation.",
    )
    topology_rooted_tbr_neighbors.add_argument("tree", type=Path)
    topology_rooted_tbr_neighbors.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_tbr_neighbors.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted TBR neighborhood report as JSON.",
    )
    _add_manifest_argument(topology_rooted_tbr_neighbors)

    topology_rooted_tbr_apply = topology_subparsers.add_parser(
        "rooted-tbr-apply",
        help="Apply one indexed rooted TBR move and write governed artifacts.",
    )
    topology_rooted_tbr_apply.add_argument("tree", type=Path)
    topology_rooted_tbr_apply.add_argument("--move-index", required=True, type=int)
    topology_rooted_tbr_apply.add_argument("--out-dir", required=True, type=Path)
    topology_rooted_tbr_apply.add_argument(
        "--json",
        action="store_true",
        help="Emit the rooted TBR move-application report as JSON.",
    )
    _add_manifest_argument(topology_rooted_tbr_apply)


def run_topology_command(args: Any) -> int:
    if args.topology_command == "clades":
        report = extract_tree_clades(
            args.tree,
            metadata_path=args.metadata,
            taxon_column=args.taxon_column,
            metadata_columns=args.metadata_columns or None,
        )
        outputs = []
        if args.out is not None:
            outputs.append(write_clade_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "clade_count": len(report.rows),
                    "metadata_column_count": len(report.metadata_columns),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "shape":
        report = summarize_tree_shape(args.tree)
        outputs = []
        if args.out is not None:
            outputs.append(write_tree_shape_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=outputs,
        )
        row = report.rows[0]
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "cherry_count": row.cherry_count,
                    "sackin_imbalance_index": row.sackin_imbalance_index,
                    "colless_imbalance_index": row.colless_imbalance_index,
                    "tree_height_edges": row.tree_height_edges,
                    "imbalance_summary": row.imbalance_summary,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "branch-lengths":
        report = analyze_branch_length_distribution(args.tree)
        outputs = []
        if args.out is not None:
            outputs.append(write_branch_length_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "branch_count": report.aggregate.branch_count,
                    "zero_length_branch_count": report.aggregate.zero_length_branch_count,
                    "negative_branch_count": report.aggregate.negative_branch_count,
                    "long_outlier_count": report.aggregate.long_outlier_count,
                    "median_branch_length": report.aggregate.median_branch_length,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "distance-reference":
        report = validate_tree_distance_reference_examples()
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[],
            outputs=[],
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[],
                outputs=outputs,
                metrics={
                    "case_count": report.case_count,
                    "external_case_count": report.external_case_count,
                    "policy_case_count": report.policy_case_count,
                    "all_passed": report.all_passed,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "support-reference":
        report = validate_support_reference_examples()
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[],
            outputs=[],
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[],
                outputs=outputs,
                metrics={
                    "case_count": report.case_count,
                    "reference_case_count": report.reference_case_count,
                    "policy_case_count": report.policy_case_count,
                    "all_passed": report.all_passed,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-nni-neighbors":
        report = enumerate_rooted_nni_neighbors(args.tree)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_nni_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "tip_count": report.tip_count,
                    "internal_node_count": report.internal_node_count,
                    "expected_neighbor_count": report.expected_neighbor_count,
                    "generated_neighbor_count": report.generated_neighbor_count,
                    "unique_neighbor_topology_count": (
                        report.unique_neighbor_topology_count
                    ),
                    "duplicate_neighbor_topology_count": len(
                        report.duplicate_neighbor_topologies
                    ),
                    "missing_tip_taxa": len(report.missing_tip_taxa),
                    "unexpected_tip_taxa": len(report.unexpected_tip_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-nni-apply":
        report = summarize_rooted_nni_move_application(args.tree, args.move_index)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_nni_move_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "selected_move_index": report.selected_move_index,
                    "available_move_count": report.available_move_count,
                    "moved_topology_changed": report.moved_topology_changed,
                    "reverse_restores_original_topology": (
                        report.reverse_restores_original_topology
                    ),
                    "missing_tip_taxa": len(report.missing_tip_taxa),
                    "unexpected_tip_taxa": len(report.unexpected_tip_taxa),
                    "node_names_preserved": report.node_names_preserved,
                    "node_metadata_preserved": report.node_metadata_preserved,
                    "edge_metadata_preserved": report.edge_metadata_preserved,
                    "branch_lengths_preserved": report.branch_lengths_preserved,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-spr-neighbors":
        budget = RootedSprEnumerationBudget(
            max_pruned_clade_count=args.max_pruned_clades,
            max_regraft_target_count_per_pruned_clade=(
                args.max_regraft_targets_per_pruned_clade
            ),
        )
        report = enumerate_rooted_spr_neighbors(args.tree, budget=budget)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_spr_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "tip_count": report.tip_count,
                    "internal_node_count": report.internal_node_count,
                    "max_pruned_clade_count": report.max_pruned_clade_count,
                    "max_regraft_target_count_per_pruned_clade": (
                        report.max_regraft_target_count_per_pruned_clade
                    ),
                    "skipped_pruned_clade_count": report.skipped_pruned_clade_count,
                    "skipped_regraft_target_count": report.skipped_regraft_target_count,
                    "generated_move_candidate_count": (
                        report.generated_move_candidate_count
                    ),
                    "identity_move_candidate_count": (
                        report.identity_move_candidate_count
                    ),
                    "self_regraft_candidate_count": (
                        report.self_regraft_candidate_count
                    ),
                    "generated_neighbor_count": report.generated_neighbor_count,
                    "unique_neighbor_topology_count": (
                        report.unique_neighbor_topology_count
                    ),
                    "duplicate_move_neighbor_topology_count": len(
                        report.duplicate_move_neighbor_topologies
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-spr-apply":
        report = summarize_rooted_spr_move_application(args.tree, args.move_index)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_spr_move_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "selected_move_index": report.selected_move_index,
                    "available_move_count": report.available_move_count,
                    "moved_topology_changed": report.moved_topology_changed,
                    "affected_clade_count": len(report.affected_clade_ids),
                    "missing_tip_taxa": len(report.missing_tip_taxa),
                    "unexpected_tip_taxa": len(report.unexpected_tip_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-tbr-neighbors":
        report = enumerate_rooted_tbr_neighbors(args.tree)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_tbr_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "tip_count": report.tip_count,
                    "internal_node_count": report.internal_node_count,
                    "generated_cut_edge_count": report.generated_cut_edge_count,
                    "generated_reconnection_count": report.generated_reconnection_count,
                    "identity_reconnection_count": report.identity_reconnection_count,
                    "generated_neighbor_count": report.generated_neighbor_count,
                    "unique_neighbor_topology_count": (
                        report.unique_neighbor_topology_count
                    ),
                    "duplicate_reconnection_neighbor_topology_count": len(
                        report.duplicate_reconnection_neighbor_topologies
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "rooted-tbr-apply":
        report = summarize_rooted_tbr_move_application(args.tree, args.move_index)
        outputs = _finalize_outputs(
            args,
            command="topology",
            inputs=[args.tree],
            outputs=list(write_rooted_tbr_move_artifacts(args.out_dir, report).values()),
        )
        _print_result(
            build_command_result(
                command="topology",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "selected_move_index": report.selected_move_index,
                    "available_move_count": report.available_move_count,
                    "left_component_tip_count": report.left_component_tip_count,
                    "right_component_tip_count": report.right_component_tip_count,
                    "moved_topology_changed": report.moved_topology_changed,
                    "reverse_move_available": report.reverse_move_available,
                    "reverse_available_move_count": report.reverse_available_move_count,
                    "missing_tip_taxa": len(report.missing_tip_taxa),
                    "unexpected_tip_taxa": len(report.unexpected_tip_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.topology_command == "root-outgroup":
        tree, report = root_tree_on_outgroup(args.tree, outgroup_taxa=list(args.taxa))
    elif args.topology_command == "reroot-midpoint":
        tree, report = reroot_tree_by_midpoint(args.tree)
    else:
        tree, report = unroot_tree(args.tree)

    output_path = write_newick(args.out, tree)
    output_paths: list[Path | str] = [output_path]
    if (
        args.topology_command in {"root-outgroup", "reroot-midpoint"}
        and getattr(args, "report_out", None) is not None
    ):
        output_paths.append(write_tree_rooting_report(args.report_out, report))
    outputs = _finalize_outputs(
        args,
        command="topology",
        inputs=[args.tree],
        outputs=output_paths,
    )
    _print_result(
        build_command_result(
            command="topology",
            inputs=[args.tree],
            outputs=outputs,
            metrics={
                "tip_count": tree.tip_count,
                "matched_taxa": len(report.matched_taxa),
                "absent_taxa": len(report.absent_taxa),
                "ingroup_taxa": len(report.ingroup_taxa),
                "outgroup_monophyletic": report.outgroup_monophyletic,
                "outgroup_mrca_extra_taxa": len(report.outgroup_mrca_extra_taxa),
                "rooted_outgroup_taxa": len(report.rooted_outgroup_taxa),
                "rooted_ingroup_taxa": len(report.rooted_ingroup_taxa),
                "midpoint_anchor_taxa": len(report.midpoint_anchor_taxa),
                "midpoint_path_length": report.midpoint_path_length,
                "midpoint_position_kind": report.midpoint_position_kind,
                "midpoint_anchor_side_taxa": len(report.midpoint_anchor_side_taxa),
                "midpoint_opposite_side_taxa": len(report.midpoint_opposite_side_taxa),
                "midpoint_suitable": report.midpoint_suitable,
                "warning_count": len(report.warnings),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
