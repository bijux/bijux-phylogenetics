from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.data_quality_stress as data_quality_stress_api
from bijux_phylogenetics.datasets.data_quality_stress import (
    export_catarrhine_data_quality_stress_panel_dataset,
    load_catarrhine_data_quality_stress_panel_dataset,
    run_catarrhine_data_quality_stress_panel_demo,
    run_catarrhine_data_quality_stress_panel_workflow,
    write_catarrhine_data_quality_stress_panel_workflow_bundle,
)

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_catarrhine_data_quality_stress_panel_dataset_exposes_packaged_surface() -> (
    None
):
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    assert dataset.dataset_id == "catarrhine_data_quality_stress_panel"
    assert dataset.label == "Catarrhine data quality stress panel"
    assert dataset.taxon_count == 6
    assert dataset.raw_trait_row_count == 7
    assert dataset.required_traits == ("body_mass_g", "gestation_days")
    assert dataset.raw_alignment_path.is_file()
    assert dataset.raw_sequence_input_path.is_file()
    assert dataset.raw_coding_sequences_path.is_file()
    assert dataset.raw_tree_path.is_file()
    assert dataset.raw_traits_path.is_file()
    assert dataset.raw_trait_mismatch_path.is_file()
    assert dataset.reference_output_root.is_dir()


def test_public_runtime_exports_include_catarrhine_data_quality_stress_panel_surface() -> (
    None
):
    assert (
        data_quality_stress_api.load_catarrhine_data_quality_stress_panel_dataset
        is load_catarrhine_data_quality_stress_panel_dataset
    )
    assert (
        data_quality_stress_api.export_catarrhine_data_quality_stress_panel_dataset
        is export_catarrhine_data_quality_stress_panel_dataset
    )
    assert (
        data_quality_stress_api.run_catarrhine_data_quality_stress_panel_workflow
        is run_catarrhine_data_quality_stress_panel_workflow
    )
    assert (
        data_quality_stress_api.write_catarrhine_data_quality_stress_panel_workflow_bundle
        is write_catarrhine_data_quality_stress_panel_workflow_bundle
    )
    assert (
        data_quality_stress_api.run_catarrhine_data_quality_stress_panel_demo
        is run_catarrhine_data_quality_stress_panel_demo
    )


def test_workflow_identifies_intended_stress_conditions_and_cleans_subset(
    tmp_path: Path,
) -> None:
    report = run_catarrhine_data_quality_stress_panel_workflow(tmp_path / "run")
    assert len(report.raw_sequence_input_validation.duplicate_identifiers) == 1
    assert len(report.raw_sequence_input_validation.illegal_characters) == 1
    assert len(report.raw_sequence_input_validation.empty_sequences) == 1
    assert [row.identifier for row in report.raw_sequence_length_outliers] == [
        "hylobates_lar",
        "macaca_mulatta",
    ]
    assert [
        (row.identifier, row.reason)
        for row in report.coding_sequence_preparation.excluded_sequences
    ] == [
        ("gorilla_gorilla", "internal-stop-codon"),
        ("pan_troglodytes", "frame-error"),
    ]
    assert report.raw_trait_mismatch_linkage.missing_from_traits == ["pan_troglodytes"]
    assert report.raw_trait_mismatch_linkage.extra_trait_taxa == ["nomascus_leucogenys"]
    assert report.raw_trait_mismatch_error == (
        "trait linkage mismatch: 1 tree taxa missing from traits and 1 trait taxa absent from tree"
    )
    assert [row.taxon for row in report.trait_duplicates] == ["hylobates_lar"]
    assert [(row.taxon, row.trait) for row in report.missing_traits] == [
        ("gorilla_gorilla", "habitat_note"),
        ("pongo_pygmaeus", "gestation_days"),
        ("hylobates_lar", "gestation_days"),
    ]
    assert [row.identifier for row in report.sequence_outliers] == ["macaca_mulatta"]
    assert report.raw_tree_validation.zero_length_branches == 1
    assert report.raw_tree_validation.negative_branch_lengths == 1
    assert [row.node for row in report.raw_tree_inspection.long_branch_outliers] == [
        "macaca_mulatta"
    ]
    assert report.dropped_taxa == ["macaca_mulatta", "pongo_pygmaeus"]
    assert report.cleaned_taxa == [
        "gorilla_gorilla",
        "homo_sapiens",
        "hylobates_lar",
        "pan_troglodytes",
    ]
    assert report.repaired_branch_nodes == ["gorilla_gorilla"]


def test_write_catarrhine_data_quality_stress_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_catarrhine_data_quality_stress_panel_workflow(tmp_path / "run")
    bundle = write_catarrhine_data_quality_stress_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.workflow_summary_path.name: bundle.workflow_summary_path,
        bundle.raw_sequence_findings_path.name: bundle.raw_sequence_findings_path,
        bundle.raw_sequence_repair_path.name: bundle.raw_sequence_repair_path,
        bundle.repaired_sequence_input_path.name: bundle.repaired_sequence_input_path,
        bundle.repaired_sequence_validation_path.name: (
            bundle.repaired_sequence_validation_path
        ),
        bundle.coding_sequence_exclusions_path.name: (
            bundle.coding_sequence_exclusions_path
        ),
        bundle.prepared_coding_sequences_path.name: (
            bundle.prepared_coding_sequences_path
        ),
        bundle.raw_trait_linkage_path.name: bundle.raw_trait_linkage_path,
        bundle.trait_duplicates_path.name: bundle.trait_duplicates_path,
        bundle.trait_missing_values_path.name: bundle.trait_missing_values_path,
        bundle.sequence_outliers_path.name: bundle.sequence_outliers_path,
        bundle.tree_issues_path.name: bundle.tree_issues_path,
        bundle.repair_actions_path.name: bundle.repair_actions_path,
        bundle.cleaned_traits_path.name: bundle.cleaned_traits_path,
        bundle.cleaned_alignment_path.name: bundle.cleaned_alignment_path,
        bundle.cleaned_tree_path.name: bundle.cleaned_tree_path,
        bundle.cleaned_linkage_path.name: bundle.cleaned_linkage_path,
        bundle.cleaned_validation_path.name: bundle.cleaned_validation_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


def test_demo_and_export_materialize_packaged_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    export_result = export_catarrhine_data_quality_stress_panel_dataset(
        tmp_path / "dataset"
    )
    assert export_result.readme_path.is_file()
    assert export_result.raw_alignment_path.is_file()
    assert export_result.raw_sequence_input_path.is_file()
    assert export_result.raw_coding_sequences_path.is_file()
    assert export_result.raw_tree_path.is_file()
    assert export_result.raw_traits_path.is_file()
    assert export_result.raw_trait_mismatch_path.is_file()
    assert len(list(export_result.expected_output_root.glob("*"))) == 18

    demo_result = run_catarrhine_data_quality_stress_panel_demo(tmp_path / "demo")
    assert demo_result.workflow_bundle.repaired_sequence_input_path.is_file()
    assert demo_result.workflow_bundle.prepared_coding_sequences_path.is_file()
    assert demo_result.workflow_bundle.cleaned_tree_path.is_file()
    assert demo_result.workflow_bundle.cleaned_alignment_path.is_file()
    assert demo_result.workflow_bundle.cleaned_traits_path.is_file()
    assert demo_result.overview_path.is_file()
    overview = demo_result.overview_path.read_text(encoding="utf-8")
    assert "repaired branch count" in overview
    assert "duplicate sequence identifiers" in overview
    assert "coding frame errors" in overview


def test_cli_demo_catarrhine_data_quality_stress_panel_json_output_reports_cleanup_review(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "stress-demo"
    exit_code = main(
        [
            "demo",
            "catarrhine-data-quality-stress-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 26
    assert payload["metrics"]["raw_taxon_count"] == 6
    assert payload["metrics"]["cleaned_taxon_count"] == 4
    assert payload["metrics"]["duplicate_sequence_identifier_count"] == 1
    assert payload["metrics"]["illegal_character_count"] == 1
    assert payload["metrics"]["empty_sequence_count"] == 1
    assert payload["metrics"]["raw_sequence_length_outlier_count"] == 2
    assert payload["metrics"]["duplicate_trait_taxon_count"] == 1
    assert payload["metrics"]["missing_trait_value_count"] == 3
    assert payload["metrics"]["sequence_outlier_count"] == 1
    assert payload["metrics"]["tree_zero_length_branch_count"] == 1
    assert payload["metrics"]["tree_negative_branch_count"] == 1
    assert payload["metrics"]["tree_long_branch_outlier_count"] == 1
    assert payload["metrics"]["coding_frame_error_count"] == 1
    assert payload["metrics"]["coding_internal_stop_count"] == 1
    assert payload["metrics"]["raw_trait_missing_from_traits_count"] == 1
    assert payload["metrics"]["raw_trait_extra_taxon_count"] == 1
    assert payload["metrics"]["dropped_taxon_count"] == 2
    assert payload["metrics"]["repaired_branch_count"] == 1
    assert payload["metrics"]["reference_output_count"] == 18
    assert payload["data"]["dataset"]["dataset_id"] == (
        "catarrhine_data_quality_stress_panel"
    )
    assert payload["data"]["dataset_export"]["raw_sequence_input_path"] == str(
        output / "dataset" / "raw" / "sequence-input.fasta"
    )
    assert payload["data"]["workflow_bundle"]["cleaned_tree_path"] == str(
        output / "workflow" / "cleaned-tree.nwk"
    )
    assert payload["data"]["workflow_bundle"]["raw_trait_linkage_path"] == str(
        output / "workflow" / "raw-trait-linkage.tsv"
    )
