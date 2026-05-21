from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.phylogenomics as phylogenomics_api
from bijux_phylogenetics.datasets.phylogenomics import (
    export_catarrhine_mitogenome_five_locus_panel_dataset,
    load_catarrhine_mitogenome_five_locus_panel_dataset,
    run_catarrhine_mitogenome_five_locus_panel_demo,
    run_catarrhine_mitogenome_five_locus_panel_workflow,
    write_catarrhine_mitogenome_five_locus_panel_workflow_bundle,
)

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_catarrhine_mitogenome_five_locus_panel_dataset_exposes_packaged_surface() -> (
    None
):
    dataset = load_catarrhine_mitogenome_five_locus_panel_dataset()
    assert dataset.dataset_id == "catarrhine_mitogenome_five_locus_panel"
    assert dataset.label == "Catarrhine mitogenome five-locus panel"
    assert dataset.taxon_count == 6
    assert dataset.locus_count == 5
    assert dataset.locus_names == (
        "mt-cox1",
        "mt-cox2",
        "mt-cox3",
        "mt-cytb",
        "mt-nd2",
    )
    assert dataset.sequence_type == "dna"
    assert dataset.taxa_path.is_file()
    assert dataset.locus_alignment_root.is_dir()
    assert dataset.reference_output_root.is_dir()
    assert "NC_012920.1" in dataset.source_accessions


@pytest.mark.slow
def test_write_catarrhine_mitogenome_five_locus_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_catarrhine_mitogenome_five_locus_panel_workflow(tmp_path / "run")
    bundle = write_catarrhine_mitogenome_five_locus_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.workflow_summary_path.name: bundle.workflow_summary_path,
        bundle.supermatrix_path.name: bundle.supermatrix_path,
        bundle.partitions_path.name: bundle.partitions_path,
        bundle.occupancy_taxa_path.name: bundle.occupancy_taxa_path,
        bundle.occupancy_loci_path.name: bundle.occupancy_loci_path,
        bundle.occupancy_matrix_path.name: bundle.occupancy_matrix_path,
        bundle.partition_summary_path.name: bundle.partition_summary_path,
        bundle.model_candidates_path.name: bundle.model_candidates_path,
        bundle.support_tree_path.name: bundle.support_tree_path,
        bundle.support_table_path.name: bundle.support_table_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_catarrhine_mitogenome_five_locus_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_catarrhine_mitogenome_five_locus_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 6
    assert result.dataset_export.taxa_path.is_file()
    assert result.dataset_export.locus_alignment_root.is_dir()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.supermatrix_path.is_file()
    assert result.workflow_bundle.partitions_path.is_file()
    assert result.workflow_bundle.support_tree_path.is_file()
    assert result.overview_path.is_file()
    assert "supported tree" in result.overview_path.read_text(encoding="utf-8")


def test_export_catarrhine_mitogenome_five_locus_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_catarrhine_mitogenome_five_locus_panel_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    locus_files = {path.name for path in result.locus_alignment_root.glob("*.fasta")}
    assert result.readme_path.is_file()
    assert result.taxa_path.is_file()
    assert result.locus_alignment_root.is_dir()
    assert len(locus_files) == 5
    assert len(expected_files) == 10
    assert "catarrhine-mitogenome-five-locus-panel.supported.tree" in expected_files


def test_public_runtime_exports_include_catarrhine_mitogenome_five_locus_panel_surface() -> (
    None
):
    assert (
        phylogenomics_api.load_catarrhine_mitogenome_five_locus_panel_dataset
        is load_catarrhine_mitogenome_five_locus_panel_dataset
    )
    assert (
        phylogenomics_api.export_catarrhine_mitogenome_five_locus_panel_dataset
        is export_catarrhine_mitogenome_five_locus_panel_dataset
    )
    assert (
        phylogenomics_api.run_catarrhine_mitogenome_five_locus_panel_workflow
        is run_catarrhine_mitogenome_five_locus_panel_workflow
    )
    assert (
        phylogenomics_api.write_catarrhine_mitogenome_five_locus_panel_workflow_bundle
        is write_catarrhine_mitogenome_five_locus_panel_workflow_bundle
    )
    assert (
        phylogenomics_api.run_catarrhine_mitogenome_five_locus_panel_demo
        is run_catarrhine_mitogenome_five_locus_panel_demo
    )


@pytest.mark.slow
def test_cli_demo_catarrhine_mitogenome_five_locus_panel_json_output_reports_multilocus_review(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "catarrhine-demo"
    exit_code = main(
        [
            "demo",
            "catarrhine-mitogenome-five-locus-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 18
    assert payload["metrics"]["taxon_count"] == 6
    assert payload["metrics"]["locus_count"] == 5
    assert payload["metrics"]["alignment_length"] == 5222
    assert payload["metrics"]["partition_count"] == 5
    assert payload["metrics"]["selected_model"] == "TIM2+F+G4"
    assert payload["metrics"]["minimum_support"] == 100.0
    assert payload["metrics"]["maximum_support"] == 100.0
    assert payload["metrics"]["weakly_supported_clade_count"] == 0
    assert payload["metrics"]["reference_output_count"] == 10
    assert payload["data"]["dataset"]["dataset_id"] == (
        "catarrhine_mitogenome_five_locus_panel"
    )
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
