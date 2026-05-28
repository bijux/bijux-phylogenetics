from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.diagnostics.assumptions import assess_tree_assumptions
from bijux_phylogenetics.diagnostics.root_to_tip import (
    compute_root_to_tip_distances,
    diagnose_root_to_tip_regression,
    diagnose_tip_date_randomization,
    diagnose_ultrametricity,
    write_root_to_tip_regression_artifacts,
    write_root_to_tip_tsv,
    write_tip_date_randomization_artifacts,
)
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path
from bijux_phylogenetics.runtime.results import build_command_result


def add_diagnose_command(subparsers: Any) -> None:
    diagnose = subparsers.add_parser(
        get_command_spec("diagnose").name, help=get_command_spec("diagnose").summary
    )
    diagnose.add_argument("target")
    diagnose.add_argument("tree", nargs="?", type=Path)
    diagnose.add_argument("--metadata", type=Path)
    diagnose.add_argument("--taxon-column")
    diagnose.add_argument("--date-column", default="date")
    diagnose.add_argument("--out", type=Path)
    diagnose.add_argument("--out-dir", type=Path)
    diagnose.add_argument("--outlier-threshold", type=float, default=2.0)
    diagnose.add_argument("--permutations", type=int, default=99)
    diagnose.add_argument("--seed", type=int, default=17)
    diagnose.add_argument("--tolerance", type=float, default=1e-6)
    diagnose.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(diagnose)


def run_diagnose_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    if args.target == "distances":
        if args.tree is None:
            parser.exit(status=2, message="diagnose distances requires a tree path\n")
        report = compute_root_to_tip_distances(args.tree)
        outputs: list[Path | str] = []
        if args.out is not None:
            output_path = write_root_to_tip_tsv(args.out, report)
            outputs.append(output_path)
        outputs = _finalize_outputs(
            args, command="diagnose", inputs=[args.tree], outputs=outputs
        )
        _print_result(
            build_command_result(
                command="diagnose",
                inputs=[args.tree],
                outputs=outputs,
                metrics={"tip_count": len(report.distances)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.target == "ultrametric":
        if args.tree is None:
            parser.exit(status=2, message="diagnose ultrametric requires a tree path\n")
        report = diagnose_ultrametricity(args.tree, tolerance=args.tolerance)
        outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree])
        _print_result(
            build_command_result(
                command="diagnose",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "tolerance": report.tolerance,
                    "max_deviation": report.max_deviation,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.target == "root-to-tip-regression":
        if args.tree is None:
            parser.exit(
                status=2,
                message="diagnose root-to-tip-regression requires a tree path\n",
            )
        if args.metadata is None:
            parser.exit(
                status=2,
                message="diagnose root-to-tip-regression requires --metadata\n",
            )
        report = diagnose_root_to_tip_regression(
            args.tree,
            args.metadata,
            taxon_column=args.taxon_column,
            date_column=args.date_column,
            outlier_threshold=args.outlier_threshold,
        )
        diagnose_inputs: list[Path | str] = [args.tree, args.metadata]
        outputs: list[Path | str] = []
        if args.out_dir is not None:
            outputs.extend(
                write_root_to_tip_regression_artifacts(args.out_dir, report).values()
            )
        outputs = _finalize_outputs(
            args,
            command="diagnose",
            inputs=diagnose_inputs,
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="diagnose",
                inputs=diagnose_inputs,
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "slope": report.slope,
                    "intercept": report.intercept,
                    "r_squared": report.r_squared,
                    "outlier_count": len(report.outliers),
                    "outlier_threshold": report.outlier_threshold,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.target == "tip-date-randomization":
        if args.tree is None:
            parser.exit(
                status=2,
                message="diagnose tip-date-randomization requires a tree path\n",
            )
        if args.metadata is None:
            parser.exit(
                status=2,
                message="diagnose tip-date-randomization requires --metadata\n",
            )
        report = diagnose_tip_date_randomization(
            args.tree,
            args.metadata,
            taxon_column=args.taxon_column,
            date_column=args.date_column,
            outlier_threshold=args.outlier_threshold,
            permutations=args.permutations,
            seed=args.seed,
        )
        diagnose_inputs: list[Path | str] = [args.tree, args.metadata]
        outputs: list[Path | str] = []
        if args.out_dir is not None:
            outputs.extend(
                write_tip_date_randomization_artifacts(args.out_dir, report).values()
            )
        outputs = _finalize_outputs(
            args,
            command="diagnose",
            inputs=diagnose_inputs,
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="diagnose",
                inputs=diagnose_inputs,
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "observed_slope": report.observed_regression.slope,
                    "observed_r_squared": report.observed_regression.r_squared,
                    "permutations": report.permutations,
                    "seed": report.seed,
                    "p_value": report.p_value,
                    "null_distribution_mean": report.null_distribution_mean,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.target == "assumptions":
        if args.tree is None:
            parser.exit(status=2, message="diagnose assumptions requires a tree path\n")
        report = assess_tree_assumptions(
            args.tree,
            metadata_path=args.metadata,
            taxon_column=args.taxon_column,
        )
        diagnose_inputs: list[Path | str] = [args.tree]
        if args.metadata is not None:
            diagnose_inputs.append(args.metadata)
        outputs = _finalize_outputs(args, command="diagnose", inputs=diagnose_inputs)
        _print_result(
            build_command_result(
                command="diagnose",
                inputs=diagnose_inputs,
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "standardized_support_count": len(
                        report.standardized_support_labels
                    ),
                    "time_tree_compatible": report.time_tree_compatible,
                    "substitution_tree_compatible": report.substitution_tree_compatible,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    tree_path = Path(args.target) if args.tree is None else args.tree
    report = diagnose_tree_path(tree_path)
    outputs = _finalize_outputs(args, command="diagnose", inputs=[tree_path])
    _print_result(
        build_command_result(
            command="diagnose",
            inputs=[tree_path],
            outputs=outputs,
            warnings=report.forensic.warnings,
            metrics={
                "tip_count": report.inspection.tip_count,
                "validity_decision": report.validation.validity_decision,
                "polytomy_count": report.validation.polytomy_count,
                "cherry_count": report.inspection.cherry_count,
                "tree_diameter": report.inspection.tree_diameter,
                "tree_quality_score": report.inspection.tree_quality_score,
                "safe_for_topology_comparison": report.forensic.safe_for_topology_comparison,
                "safe_for_time_tree_analysis": report.forensic.safe_for_time_tree_analysis,
                "safe_for_comparative_methods": report.forensic.safe_for_comparative_methods,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
