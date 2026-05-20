from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_distance_tree_method_argument,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_distance_method_maturity,
    bootstrap_distance_trees,
    build_distance_method_report,
    build_distance_tree,
    compare_distance_gap_policies,
    compare_distance_models,
    compare_distance_tree_to_reference_tree,
    compare_distance_tree_topologies,
    inspect_distance_matrix_quality,
    summarize_distance_bootstrap_support,
    write_distance_bootstrap_support,
    write_distance_reproducibility_bundle,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import write_tree_set

from .matrix import add_distance_matrix_command, run_distance_matrix_command
from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_alignment_distance_commands(alignment_subparsers: Any) -> None:
    add_distance_matrix_command(alignment_subparsers)

    alignment_distance_quality = alignment_subparsers.add_parser(
        "distance-quality",
        help="Inspect saturation, divergence, and low-information risks in a computed distance matrix.",
    )
    alignment_distance_quality.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_quality)
    add_gap_handling_option(alignment_distance_quality)
    add_ambiguity_policy_option(alignment_distance_quality)
    alignment_distance_quality.add_argument(
        "--json", action="store_true", help="Emit the diagnostics as JSON."
    )
    _add_manifest_argument(alignment_distance_quality)

    alignment_distance_assumptions = alignment_subparsers.add_parser(
        "distance-assumptions",
        help="Audit NJ and UPGMA assumptions, including UPGMA ultrametric compatibility.",
    )
    alignment_distance_assumptions.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_assumptions)
    add_gap_handling_option(alignment_distance_assumptions)
    add_ambiguity_policy_option(alignment_distance_assumptions)
    alignment_distance_assumptions.add_argument(
        "--json", action="store_true", help="Emit the assumption audit as JSON."
    )
    _add_manifest_argument(alignment_distance_assumptions)

    alignment_build_tree = alignment_subparsers.add_parser(
        "build-tree",
        help="Build a neighbor-joining or UPGMA tree from a DNA distance matrix.",
    )
    alignment_build_tree.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_build_tree)
    add_distance_model_option(alignment_build_tree)
    add_gap_handling_option(alignment_build_tree)
    add_ambiguity_policy_option(alignment_build_tree)
    alignment_build_tree.add_argument("--out", required=True, type=Path)
    alignment_build_tree.add_argument(
        "--json", action="store_true", help="Emit the build report as JSON."
    )
    _add_manifest_argument(alignment_build_tree)

    alignment_compare_distance_trees = alignment_subparsers.add_parser(
        "compare-distance-trees",
        help="Compare NJ and UPGMA topologies built from the same DNA alignment.",
    )
    alignment_compare_distance_trees.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_compare_distance_trees)
    add_gap_handling_option(alignment_compare_distance_trees)
    add_ambiguity_policy_option(alignment_compare_distance_trees)
    alignment_compare_distance_trees.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(alignment_compare_distance_trees)

    alignment_compare_distance_reference = alignment_subparsers.add_parser(
        "compare-distance-to-tree",
        help="Compare one built distance tree against an external inferred or reviewer reference tree.",
    )
    alignment_compare_distance_reference.add_argument("alignment", type=Path)
    alignment_compare_distance_reference.add_argument("reference_tree", type=Path)
    _add_distance_tree_method_argument(alignment_compare_distance_reference)
    add_distance_model_option(alignment_compare_distance_reference)
    add_gap_handling_option(alignment_compare_distance_reference)
    add_ambiguity_policy_option(alignment_compare_distance_reference)
    alignment_compare_distance_reference.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(alignment_compare_distance_reference)

    alignment_bootstrap_tree = alignment_subparsers.add_parser(
        "bootstrap-tree",
        help="Bootstrap a distance tree by resampling alignment sites with replacement.",
    )
    alignment_bootstrap_tree.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_bootstrap_tree)
    add_distance_model_option(alignment_bootstrap_tree)
    add_gap_handling_option(alignment_bootstrap_tree)
    add_ambiguity_policy_option(alignment_bootstrap_tree)
    alignment_bootstrap_tree.add_argument("--replicates", type=int, default=100)
    alignment_bootstrap_tree.add_argument("--seed", type=int, default=1)
    alignment_bootstrap_tree.add_argument(
        "--support-out", type=Path, help="Write bootstrap clade support as TSV."
    )
    alignment_bootstrap_tree.add_argument(
        "--tree-set-out", type=Path, help="Write bootstrap replicate trees as Newick."
    )
    alignment_bootstrap_tree.add_argument(
        "--json", action="store_true", help="Emit the bootstrap report as JSON."
    )
    _add_manifest_argument(alignment_bootstrap_tree)

    alignment_bootstrap_summary = alignment_subparsers.add_parser(
        "distance-support-summary",
        help="Summarize consensus clade support across distance-bootstrap replicates.",
    )
    alignment_bootstrap_summary.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_bootstrap_summary)
    add_distance_model_option(alignment_bootstrap_summary)
    add_gap_handling_option(alignment_bootstrap_summary)
    add_ambiguity_policy_option(alignment_bootstrap_summary)
    alignment_bootstrap_summary.add_argument("--replicates", type=int, default=25)
    alignment_bootstrap_summary.add_argument("--seed", type=int, default=1)
    alignment_bootstrap_summary.add_argument(
        "--json", action="store_true", help="Emit the support summary as JSON."
    )
    _add_manifest_argument(alignment_bootstrap_summary)

    alignment_distance_models = alignment_subparsers.add_parser(
        "distance-models",
        help="Compare all supported distance models on the same alignment.",
    )
    alignment_distance_models.add_argument("alignment", type=Path)
    add_gap_handling_option(alignment_distance_models)
    add_ambiguity_policy_option(alignment_distance_models)
    alignment_distance_models.add_argument(
        "--json", action="store_true", help="Emit the model comparison as JSON."
    )
    _add_manifest_argument(alignment_distance_models)

    alignment_distance_gap = alignment_subparsers.add_parser(
        "distance-gap-sensitivity",
        help="Compare pairwise versus complete deletion for the same distance workflow.",
    )
    alignment_distance_gap.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_gap)
    add_ambiguity_policy_option(alignment_distance_gap)
    alignment_distance_gap.add_argument(
        "--json", action="store_true", help="Emit the gap-policy sensitivity as JSON."
    )
    _add_manifest_argument(alignment_distance_gap)

    alignment_distance_suitability = alignment_subparsers.add_parser(
        "distance-suitability",
        help="Emit the explicit suitability decision for distance-method use on one alignment.",
    )
    alignment_distance_suitability.add_argument("alignment", type=Path)
    add_distance_model_option(alignment_distance_suitability)
    add_gap_handling_option(alignment_distance_suitability)
    add_ambiguity_policy_option(alignment_distance_suitability)
    alignment_distance_suitability.add_argument(
        "--json", action="store_true", help="Emit the suitability decision as JSON."
    )
    _add_manifest_argument(alignment_distance_suitability)

    alignment_distance_method_report = alignment_subparsers.add_parser(
        "distance-method-report",
        help="Build a structured distance-method report with support, model, and gap-sensitivity sections.",
    )
    alignment_distance_method_report.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_method_report)
    add_distance_model_option(alignment_distance_method_report)
    add_gap_handling_option(alignment_distance_method_report)
    add_ambiguity_policy_option(alignment_distance_method_report)
    alignment_distance_method_report.add_argument("--replicates", type=int, default=25)
    alignment_distance_method_report.add_argument("--seed", type=int, default=1)
    alignment_distance_method_report.add_argument(
        "--json", action="store_true", help="Emit the structured report as JSON."
    )
    _add_manifest_argument(alignment_distance_method_report)

    alignment_distance_maturity = alignment_subparsers.add_parser(
        "distance-maturity",
        help="Run the distance-method maturity gate for one alignment.",
    )
    alignment_distance_maturity.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_maturity)
    add_distance_model_option(alignment_distance_maturity)
    add_gap_handling_option(alignment_distance_maturity)
    add_ambiguity_policy_option(alignment_distance_maturity)
    alignment_distance_maturity.add_argument("--replicates", type=int, default=25)
    alignment_distance_maturity.add_argument("--seed", type=int, default=1)
    alignment_distance_maturity.add_argument(
        "--json", action="store_true", help="Emit the maturity gate as JSON."
    )
    _add_manifest_argument(alignment_distance_maturity)

    alignment_distance_bundle = alignment_subparsers.add_parser(
        "distance-bundle",
        help="Write a reproducibility bundle for one distance-analysis workflow.",
    )
    alignment_distance_bundle.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_bundle)
    add_distance_model_option(alignment_distance_bundle)
    add_gap_handling_option(alignment_distance_bundle)
    add_ambiguity_policy_option(alignment_distance_bundle)
    alignment_distance_bundle.add_argument("--replicates", type=int, default=100)
    alignment_distance_bundle.add_argument("--seed", type=int, default=1)
    alignment_distance_bundle.add_argument("--out-dir", required=True, type=Path)
    alignment_distance_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(alignment_distance_bundle)


def run_alignment_distance_command(args: Any) -> int | None:
    matrix_result = run_distance_matrix_command(args)
    if matrix_result is not None:
        return matrix_result

    if args.alignment_command == "distance-quality":
        report = inspect_distance_matrix_quality(
            args.alignment,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "taxon_count": report.taxon_count,
                    "saturated_pair_count": len(report.saturated_pairs),
                    "low_information_pair_count": len(report.low_information_pairs),
                    "decision": report.method_assessment.decision,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-suitability":
        report = inspect_distance_matrix_quality(
            args.alignment,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.method_assessment.reasons,
                metrics={
                    "decision": report.method_assessment.decision,
                    "reason_count": len(report.method_assessment.reasons),
                },
                data=report.method_assessment,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-assumptions":
        report = assess_distance_method_assumptions(
            args.alignment,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "taxon_count": report.taxon_count,
                    "ultrametric_compatible": report.ultrametric_compatible,
                    "upgma_violation_count": len(report.upgma_ultrametric_violations),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "build-tree":
        tree, report = build_distance_tree(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        output_path = write_newick(args.out, tree)
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "pair_count": report.pair_count,
                    "method": report.method,
                    "ambiguity_policy": report.ambiguity_policy,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "compare-distance-trees":
        report = compare_distance_tree_topologies(
            args.alignment,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "shared_taxa": len(report.shared_taxa),
                    "robinson_foulds_distance": report.robinson_foulds_distance,
                    "same_unrooted_topology": report.same_unrooted_topology,
                    "ambiguity_policy": report.ambiguity_policy,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "compare-distance-to-tree":
        report = compare_distance_tree_to_reference_tree(
            args.alignment,
            args.reference_tree,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment, args.reference_tree],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment, args.reference_tree],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "topology_equal": report.topology.topology_equal,
                    "same_unrooted_topology": report.topology.same_unrooted_topology,
                    "shared_taxa": len(report.topology.shared_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "bootstrap-tree":
        trees, report = bootstrap_distance_trees(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs: list[Path | str] = []
        if args.support_out is not None:
            outputs.append(write_distance_bootstrap_support(args.support_out, report))
        if args.tree_set_out is not None:
            outputs.append(write_tree_set(args.tree_set_out, trees))
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "replicate_count": report.tree_count,
                    "support_row_count": len(report.support),
                    "method": report.method,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-support-summary":
        report = summarize_distance_bootstrap_support(
            bootstrap_distance_trees(
                args.alignment,
                method=args.method,
                model=args.model,
                gap_handling=args.gap_handling,
                ambiguity_policy=args.ambiguity_policy,
                replicates=args.replicates,
                seed=args.seed,
            )[1]
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "clade_count": report.clade_count,
                    "weak_clade_count": report.weak_clade_count,
                    "replicates": report.replicates,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-models":
        report = compare_distance_models(
            args.alignment,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model_count": len(report.rows),
                    "alphabet": report.inferred_alphabet,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-gap-sensitivity":
        report = compare_distance_gap_policies(
            args.alignment,
            model=args.model,
            ambiguity_policy=args.ambiguity_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "changed_pair_count": report.changed_pair_count,
                    "pair_count": report.pair_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-method-report":
        report = build_distance_method_report(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            bootstrap_replicates=args.replicates,
            bootstrap_seed=args.seed,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.maturity_gate.warnings,
                metrics={
                    "method": report.method,
                    "decision": report.maturity_gate.decision,
                    "bootstrap_clade_count": report.bootstrap_summary.clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-maturity":
        report = assess_distance_method_maturity(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            bootstrap_replicates=args.replicates,
            bootstrap_seed=args.seed,
            validate_bundle=True,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "decision": report.decision,
                    "check_count": len(report.checks),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "distance-bundle":
        report = write_distance_reproducibility_bundle(
            args.out_dir,
            alignment_path=args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=list(report.files),
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "file_count": len(report.files),
                    "replicates": report.replicates,
                    "method": report.method,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
