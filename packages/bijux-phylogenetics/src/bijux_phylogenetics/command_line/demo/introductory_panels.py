from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_entries,
    emit_demo_result,
)
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.datasets import (
    run_avian_reproductive_trait_demo,
    run_central_european_seashore_flora_demo,
    run_primate_comparative_demo,
)


def add_introductory_demo_commands(demo_subparsers: Any) -> None:
    demo_run = demo_subparsers.add_parser(
        "run", help="Run the repository capability demo workflow."
    )
    demo_run.add_argument("--out", required=True, type=Path)
    demo_run.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_run)

    demo_primate = demo_subparsers.add_parser(
        "primate-comparative",
        help="Materialize the packaged primate dataset and comparative workflow outputs.",
    )
    demo_primate.add_argument("--out", required=True, type=Path)
    demo_primate.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_primate)

    demo_birds = demo_subparsers.add_parser(
        "avian-reproductive-traits",
        help="Materialize the packaged avian reproductive dataset and workflow outputs.",
    )
    demo_birds.add_argument("--out", required=True, type=Path)
    demo_birds.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_birds)

    demo_plants = demo_subparsers.add_parser(
        "central-european-seashore-flora",
        help="Materialize the packaged Central European plant dataset and workflow outputs.",
    )
    demo_plants.add_argument("--out", required=True, type=Path)
    demo_plants.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_plants)


def run_introductory_demo_command(args: Any) -> int | None:
    if args.demo_command == "run":
        result = run_capability_demo(args.out)
        return emit_demo_result(
            args,
            outputs=[
                result.tree_report,
                result.dataset_report,
                result.phylo_inputs_report,
                result.comparison_report,
                result.capability_summary,
            ],
            metrics={"artifact_count": 5},
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "primate-comparative":
        result = run_primate_comparative_demo(args.out)
        outputs = [
            result.dataset_export.readme_path,
            result.dataset_export.tree_path,
            result.dataset_export.traits_path,
            result.workflow_bundle.summary_path,
            result.workflow_bundle.pgls_lambda_profile_path,
            result.workflow_bundle.brownian_summary_path,
            result.workflow_bundle.ou_summary_path,
            result.workflow_bundle.signal_summary_path,
            result.workflow_bundle.signal_permutations_path,
            result.workflow_bundle.continuous_ancestral_summary_path,
            result.workflow_bundle.continuous_ancestral_uncertainty_path,
            result.workflow_bundle.discrete_ancestral_summary_path,
            result.workflow_bundle.discrete_ancestral_probability_path,
            result.overview_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "dataset_taxon_count": result.dataset.taxon_count,
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command == "avian-reproductive-traits":
        result = run_avian_reproductive_trait_demo(args.out)
        outputs = [
            result.dataset_export.readme_path,
            result.dataset_export.tree_path,
            result.dataset_export.traits_path,
            result.workflow_bundle.summary_path,
            result.workflow_bundle.pgls_lambda_profile_path,
            result.workflow_bundle.brownian_summary_path,
            result.workflow_bundle.ou_summary_path,
            result.workflow_bundle.signal_summary_path,
            result.workflow_bundle.signal_permutations_path,
            result.workflow_bundle.continuous_ancestral_summary_path,
            result.workflow_bundle.continuous_ancestral_uncertainty_path,
            result.workflow_bundle.discrete_ancestral_summary_path,
            result.workflow_bundle.discrete_ancestral_probability_path,
            result.workflow_bundle.clade_summary_path,
            result.workflow_bundle.clade_rows_path,
            result.overview_path,
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "dataset_taxon_count": result.dataset.taxon_count,
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command != "central-european-seashore-flora":
        return None

    result = run_central_european_seashore_flora_demo(args.out)
    outputs = [
        result.dataset_export.readme_path,
        result.dataset_export.tree_path,
        result.dataset_export.traits_path,
        result.workflow_bundle.summary_path,
        result.workflow_bundle.pgls_lambda_profile_path,
        result.workflow_bundle.brownian_summary_path,
        result.workflow_bundle.ou_summary_path,
        result.workflow_bundle.signal_summary_path,
        result.workflow_bundle.signal_permutations_path,
        result.workflow_bundle.continuous_ancestral_summary_path,
        result.workflow_bundle.continuous_ancestral_uncertainty_path,
        result.workflow_bundle.discrete_ancestral_summary_path,
        result.workflow_bundle.discrete_ancestral_probability_path,
        result.workflow_bundle.clade_summary_path,
        result.workflow_bundle.clade_rows_path,
        result.overview_path,
    ]
    return emit_demo_result(
        args,
        outputs=outputs,
        metrics={
            "artifact_count": len(outputs),
            "dataset_taxon_count": result.dataset.taxon_count,
            "reference_output_count": count_expected_output_entries(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
