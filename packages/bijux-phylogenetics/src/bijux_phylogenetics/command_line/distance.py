from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_distance_tree_method_argument,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    assess_imported_distance_method_assumptions,
    build_tree_from_imported_distance_matrix,
    fit_fitch_margoliash_tree_from_imported_distance_matrix,
    fit_minimum_evolution_tree_from_imported_distance_matrix,
    inspect_imported_distance_matrix_quality,
    search_balanced_minimum_evolution_nni_from_imported_distance_matrix,
    validate_distance_reference_examples,
    validate_imported_distance_matrix,
    write_balanced_minimum_evolution_nni_artifacts,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.reports.service import (
    distance_method_limitations,
    render_distance_report,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_distance_commands(subparsers: Any) -> None:
    distance = subparsers.add_parser(
        get_command_spec("distance").name, help=get_command_spec("distance").summary
    )
    distance_subparsers = distance.add_subparsers(
        dest="distance_command", required=True
    )

    distance_validate = distance_subparsers.add_parser(
        "validate", help="Validate an imported long-form distance matrix."
    )
    distance_validate.add_argument("matrix", type=Path)
    distance_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(distance_validate)

    distance_quality = distance_subparsers.add_parser(
        "quality",
        help="Audit structural, saturation, and low-information risks for an imported distance matrix.",
    )
    distance_quality.add_argument("matrix", type=Path)
    distance_quality.add_argument(
        "--json", action="store_true", help="Emit the quality report as JSON."
    )
    _add_manifest_argument(distance_quality)

    distance_assumptions = distance_subparsers.add_parser(
        "assumptions",
        help="Audit NJ and UPGMA assumptions for an imported distance matrix.",
    )
    distance_assumptions.add_argument("matrix", type=Path)
    distance_assumptions.add_argument(
        "--json", action="store_true", help="Emit the assumption audit as JSON."
    )
    _add_manifest_argument(distance_assumptions)

    distance_build_tree = distance_subparsers.add_parser(
        "build-tree",
        help="Build one owned distance tree from an imported distance matrix.",
    )
    distance_build_tree.add_argument("matrix", type=Path)
    _add_distance_tree_method_argument(distance_build_tree)
    distance_build_tree.add_argument("--out", required=True, type=Path)
    distance_build_tree.add_argument(
        "--json", action="store_true", help="Emit the build report as JSON."
    )
    _add_manifest_argument(distance_build_tree)

    distance_minimum_evolution = distance_subparsers.add_parser(
        "minimum-evolution",
        help="Fit one fixed tree topology to an imported distance matrix and score it by total fitted branch length.",
    )
    distance_minimum_evolution.add_argument("matrix", type=Path)
    distance_minimum_evolution.add_argument("tree", type=Path)
    distance_minimum_evolution.add_argument("--out", required=True, type=Path)
    distance_minimum_evolution.add_argument(
        "--json", action="store_true", help="Emit the minimum-evolution score as JSON."
    )
    _add_manifest_argument(distance_minimum_evolution)

    distance_fitch_margoliash = distance_subparsers.add_parser(
        "fitch-margoliash",
        help="Fit one fixed tree topology to an imported distance matrix by classical Fitch-Margoliash weighted least squares.",
    )
    distance_fitch_margoliash.add_argument("matrix", type=Path)
    distance_fitch_margoliash.add_argument("tree", type=Path)
    distance_fitch_margoliash.add_argument("--out", required=True, type=Path)
    distance_fitch_margoliash.add_argument(
        "--weighting-power",
        type=float,
        default=2.0,
        help="Power used in the classical Fitch-Margoliash pair weight d^(-p).",
    )
    distance_fitch_margoliash.add_argument(
        "--json", action="store_true", help="Emit the Fitch-Margoliash fit as JSON."
    )
    _add_manifest_argument(distance_fitch_margoliash)

    distance_bme_nni = distance_subparsers.add_parser(
        "bme-nni-search",
        help="Start from NJ or BIONJ and hill-climb one imported distance matrix by rooted NNI under the balanced minimum-evolution objective.",
    )
    distance_bme_nni.add_argument("matrix", type=Path)
    distance_bme_nni.add_argument(
        "--start-method",
        required=True,
        choices=("neighbor-joining", "bionj"),
        help="Distance-tree initializer used before the rooted NNI hill-climb.",
    )
    distance_bme_nni.add_argument("--out-dir", required=True, type=Path)
    distance_bme_nni.add_argument(
        "--json", action="store_true", help="Emit the BME NNI search report as JSON."
    )
    _add_manifest_argument(distance_bme_nni)

    distance_report = distance_subparsers.add_parser(
        "report",
        help="Render an HTML report for an imported distance matrix.",
    )
    distance_report.add_argument("matrix", type=Path)
    distance_report.add_argument("--out", required=True, type=Path)
    distance_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(distance_report)

    distance_explain = distance_subparsers.add_parser(
        "explain",
        help="Explain why distance-based tree building is approximate.",
    )
    distance_explain.add_argument("matrix", type=Path)
    distance_explain.add_argument(
        "--json", action="store_true", help="Emit the explanation as JSON."
    )
    _add_manifest_argument(distance_explain)

    distance_reference = distance_subparsers.add_parser(
        "reference",
        help="Validate built-in reference examples for core distance calculations.",
    )
    distance_reference.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(distance_reference)


def run_distance_command(args: Any) -> int:
    if args.distance_command == "reference":
        report = validate_distance_reference_examples()
        outputs = _finalize_outputs(args, command="distance", inputs=[])
        _print_result(
            build_command_result(
                command="distance",
                inputs=[],
                outputs=outputs,
                metrics={
                    "case_count": len(report.observations),
                    "all_passed": report.all_passed,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "assumptions":
        report = assess_imported_distance_method_assumptions(args.matrix)
        outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
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
    if args.distance_command == "validate":
        report = validate_imported_distance_matrix(args.matrix)
        outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "taxon_count": len(report.identifiers),
                    "pair_count": report.pair_count,
                    "complete": report.complete,
                    "symmetric": report.symmetric,
                    "nonmetric_observation_count": len(report.nonmetric_observations),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "quality":
        report = inspect_imported_distance_matrix_quality(args.matrix)
        outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "taxon_count": len(report.validation.identifiers),
                    "missing_pair_count": len(report.validation.missing_pairs),
                    "saturated_pair_count": len(report.saturated_pairs),
                    "low_information_pair_count": len(report.low_information_pairs),
                    "saturation_audit_scale": report.saturation_audit_scale,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "build-tree":
        tree, report = build_tree_from_imported_distance_matrix(
            args.matrix,
            method=args.method,
        )
        output_path = write_newick(args.out, tree)
        outputs = _finalize_outputs(
            args,
            command="distance",
            inputs=[args.matrix],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
                outputs=outputs,
                metrics={
                    "method": report.method,
                    "taxon_count": report.taxon_count,
                    "pair_count": report.pair_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "minimum-evolution":
        tree, report = fit_minimum_evolution_tree_from_imported_distance_matrix(
            args.matrix,
            args.tree,
        )
        output_path = write_newick(args.out, tree)
        outputs = _finalize_outputs(
            args,
            command="distance",
            inputs=[args.matrix, args.tree],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix, args.tree],
                outputs=outputs,
                metrics={
                    "criterion": "minimum-evolution",
                    "taxon_count": len(report.taxa),
                    "pair_count": report.pair_count,
                    "branch_count": report.branch_count,
                    "minimum_evolution_score": report.minimum_evolution_score,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "fitch-margoliash":
        tree, report = fit_fitch_margoliash_tree_from_imported_distance_matrix(
            args.matrix,
            args.tree,
            weighting_power=args.weighting_power,
        )
        output_path = write_newick(args.out, tree)
        outputs = _finalize_outputs(
            args,
            command="distance",
            inputs=[args.matrix, args.tree],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix, args.tree],
                outputs=outputs,
                metrics={
                    "criterion": "fitch-margoliash",
                    "taxon_count": len(report.taxa),
                    "pair_count": report.pair_count,
                    "branch_count": report.branch_count,
                    "weighting_power": report.weighting_power,
                    "residual_sum_squares": report.residual_sum_squares,
                    "weighted_residual_sum_squares": report.weighted_residual_sum_squares,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "bme-nni-search":
        report = search_balanced_minimum_evolution_nni_from_imported_distance_matrix(
            args.matrix,
            start_method=args.start_method,
        )
        artifact_paths = write_balanced_minimum_evolution_nni_artifacts(
            args.out_dir,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="distance",
            inputs=[args.matrix],
            outputs=list(artifact_paths.values()),
        )
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
                outputs=outputs,
                metrics={
                    "algorithm": report.algorithm,
                    "start_method": report.start_method,
                    "taxon_count": report.taxon_count,
                    "pair_count": report.pair_count,
                    "start_score": report.start_score,
                    "final_score": report.final_score,
                    "accepted_move_count": report.accepted_move_count,
                    "evaluated_neighbor_count": report.evaluated_neighbor_count,
                    "stopping_reason": report.stopping_reason,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.distance_command == "report":
        report = render_distance_report(out_path=args.out, matrix_path=args.matrix)
        outputs = _finalize_outputs(
            args,
            command="distance",
            inputs=[args.matrix],
            outputs=[args.out],
        )
        _print_result(
            build_command_result(
                command="distance",
                inputs=[args.matrix],
                outputs=outputs,
                metrics={"section_count": len(report.machine_manifest["sections"])},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    limitations = distance_method_limitations()
    outputs = _finalize_outputs(args, command="distance", inputs=[args.matrix])
    _print_result(
        build_command_result(
            command="distance",
            inputs=[args.matrix],
            outputs=outputs,
            metrics={"limitation_count": len(limitations)},
            data={"limitations": limitations},
        ),
        json_output=args.json,
    )
    return 0
