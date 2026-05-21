from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_entries,
    emit_demo_result,
)
from bijux_phylogenetics.datasets.data_quality_stress import (
    run_catarrhine_data_quality_stress_panel_demo,
)


def add_quality_demo_commands(demo_subparsers: Any) -> None:
    demo_catarrhine_stress = demo_subparsers.add_parser(
        "catarrhine-data-quality-stress-panel",
        help="Materialize the packaged catarrhine dirty-data stress dataset and rerun the governed audit and cleanup outputs.",
    )
    demo_catarrhine_stress.add_argument("--out", required=True, type=Path)
    demo_catarrhine_stress.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_stress)


def run_quality_demo_command(args: Any) -> int | None:
    if args.demo_command != "catarrhine-data-quality-stress-panel":
        return None

    result = run_catarrhine_data_quality_stress_panel_demo(args.out)
    outputs = [
        result.dataset_export.readme_path,
        result.dataset_export.raw_alignment_path,
        result.dataset_export.raw_sequence_input_path,
        result.dataset_export.raw_coding_sequences_path,
        result.dataset_export.raw_tree_path,
        result.dataset_export.raw_traits_path,
        result.dataset_export.raw_trait_mismatch_path,
        result.workflow_bundle.workflow_summary_path,
        result.workflow_bundle.raw_sequence_findings_path,
        result.workflow_bundle.raw_sequence_repair_path,
        result.workflow_bundle.repaired_sequence_input_path,
        result.workflow_bundle.repaired_sequence_validation_path,
        result.workflow_bundle.coding_sequence_exclusions_path,
        result.workflow_bundle.prepared_coding_sequences_path,
        result.workflow_bundle.raw_trait_linkage_path,
        result.workflow_bundle.trait_duplicates_path,
        result.workflow_bundle.trait_missing_values_path,
        result.workflow_bundle.sequence_outliers_path,
        result.workflow_bundle.tree_issues_path,
        result.workflow_bundle.repair_actions_path,
        result.workflow_bundle.cleaned_traits_path,
        result.workflow_bundle.cleaned_alignment_path,
        result.workflow_bundle.cleaned_tree_path,
        result.workflow_bundle.cleaned_linkage_path,
        result.workflow_bundle.cleaned_validation_path,
        result.overview_path,
    ]
    return emit_demo_result(
        args,
        outputs=outputs,
        metrics={
            "artifact_count": len(outputs),
            "raw_taxon_count": result.workflow_bundle.raw_taxon_count,
            "cleaned_taxon_count": result.workflow_bundle.cleaned_taxon_count,
            "duplicate_sequence_identifier_count": (
                result.workflow_bundle.duplicate_sequence_identifier_count
            ),
            "illegal_character_count": (result.workflow_bundle.illegal_character_count),
            "empty_sequence_count": result.workflow_bundle.empty_sequence_count,
            "raw_sequence_length_outlier_count": (
                result.workflow_bundle.raw_sequence_length_outlier_count
            ),
            "duplicate_trait_taxon_count": (
                result.workflow_bundle.duplicate_trait_taxon_count
            ),
            "missing_trait_value_count": (
                result.workflow_bundle.missing_trait_value_count
            ),
            "sequence_outlier_count": result.workflow_bundle.sequence_outlier_count,
            "tree_zero_length_branch_count": (
                result.workflow_bundle.tree_zero_length_branch_count
            ),
            "tree_negative_branch_count": (
                result.workflow_bundle.tree_negative_branch_count
            ),
            "tree_long_branch_outlier_count": (
                result.workflow_bundle.tree_long_branch_outlier_count
            ),
            "coding_frame_error_count": (
                result.workflow_bundle.coding_frame_error_count
            ),
            "coding_internal_stop_count": (
                result.workflow_bundle.coding_internal_stop_count
            ),
            "raw_trait_missing_from_traits_count": (
                result.workflow_bundle.raw_trait_missing_from_traits_count
            ),
            "raw_trait_extra_taxon_count": (
                result.workflow_bundle.raw_trait_extra_taxon_count
            ),
            "dropped_taxon_count": result.workflow_bundle.dropped_taxon_count,
            "repaired_branch_count": result.workflow_bundle.repaired_branch_count,
            "reference_output_count": count_expected_output_entries(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
