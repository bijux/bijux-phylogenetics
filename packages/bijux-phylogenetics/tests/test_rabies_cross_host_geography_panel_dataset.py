from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.datasets import (
    export_rabies_cross_host_geography_panel_dataset,
    load_rabies_cross_host_geography_panel_dataset,
    run_rabies_cross_host_geography_panel_demo,
    run_rabies_cross_host_geography_panel_workflow,
    write_rabies_cross_host_geography_panel_workflow_bundle,
)
from tests.support.external_engines import require_alignment_engine_executables


def _stable_generated_outputs(bundle: object) -> dict[str, Path]:
    workflow_bundle = bundle
    return {
        workflow_bundle.workflow_summary_path.name: workflow_bundle.workflow_summary_path,
        workflow_bundle.alignment_path.name: workflow_bundle.alignment_path,
        workflow_bundle.trimmed_alignment_path.name: workflow_bundle.trimmed_alignment_path,
        workflow_bundle.tree_path.name: workflow_bundle.tree_path,
        workflow_bundle.rooting_report_path.name: workflow_bundle.rooting_report_path,
        workflow_bundle.model_table_path.name: workflow_bundle.model_table_path,
        workflow_bundle.support_table_path.name: workflow_bundle.support_table_path,
        workflow_bundle.host_switch_summary_path.name: workflow_bundle.host_switch_summary_path,
        workflow_bundle.host_state_nodes_path.name: workflow_bundle.host_state_nodes_path,
        workflow_bundle.host_switch_branches_path.name: workflow_bundle.host_switch_branches_path,
        workflow_bundle.host_switch_counts_path.name: workflow_bundle.host_switch_counts_path,
        workflow_bundle.host_switch_fits_path.name: workflow_bundle.host_switch_fits_path,
        workflow_bundle.host_switch_unsupported_path.name: workflow_bundle.host_switch_unsupported_path,
        workflow_bundle.host_switch_exclusions_path.name: workflow_bundle.host_switch_exclusions_path,
        workflow_bundle.final_report_path.name: workflow_bundle.final_report_path,
        workflow_bundle.final_manifest_path.name: workflow_bundle.final_manifest_path,
        "biogeography/biogeography-report.html": workflow_bundle.biogeography_report_path,
        "biogeography/ancestral-region-tree.svg": workflow_bundle.biogeography_tree_figure_path,
        "biogeography/geographic-region-map.html": workflow_bundle.biogeography_map_path,
        "biogeography/summary.tsv": workflow_bundle.biogeography_output_root / "summary.tsv",
        "biogeography/region-counts.tsv": workflow_bundle.biogeography_output_root / "region-counts.tsv",
        "biogeography/ancestral-regions.tsv": workflow_bundle.biogeography_output_root / "ancestral-regions.tsv",
        "biogeography/transition-matrix.tsv": workflow_bundle.biogeography_output_root / "transition-matrix.tsv",
        "biogeography/event-table.tsv": workflow_bundle.biogeography_output_root / "event-table.tsv",
        "biogeography/map-markers.tsv": workflow_bundle.biogeography_output_root / "map-markers.tsv",
        "biogeography/map-lines.tsv": workflow_bundle.biogeography_output_root / "map-lines.tsv",
        "biogeography/exclusions.tsv": workflow_bundle.biogeography_output_root / "exclusions.tsv",
    }


def test_load_rabies_cross_host_geography_panel_dataset_exposes_packaged_surface() -> None:
    dataset = load_rabies_cross_host_geography_panel_dataset()
    assert dataset.dataset_id == "rabies_cross_host_geography_panel"
    assert dataset.label == "Rabies cross-host geography panel"
    assert dataset.sequence_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_prefix == "rabies-cross-host-geography-panel"
    assert dataset.host_trait == "host_group"
    assert dataset.geography_trait == "region_group"
    assert dataset.host_model == "ard"
    assert dataset.geography_model == "ard"
    assert dataset.iqtree_seed == 1
    assert dataset.iqtree_threads == 1
    assert dataset.bootstrap_replicates == 1000
    assert dataset.outgroup_taxa == ("bat_chile_rv108",)
    assert dataset.observed_host_group_count == 3
    assert dataset.observed_region_group_count == 5
    assert dataset.sequences_path.is_file()
    assert dataset.metadata_path.is_file()
    assert dataset.centroids_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "MG458305" in dataset.source_accessions


def test_write_rabies_cross_host_geography_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_rabies_cross_host_geography_panel_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = _stable_generated_outputs(bundle)
    expected_files = {
        str(path.relative_to(expected_root))
        for path in expected_root.rglob("*")
        if path.is_file()
    }
    assert expected_files == set(generated)
    for relative_name, generated_path in generated.items():
        assert generated_path.read_text(encoding="utf-8") == (
            expected_root / relative_name
        ).read_text(encoding="utf-8")


def test_run_rabies_cross_host_geography_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    result = run_rabies_cross_host_geography_panel_demo(
        tmp_path / "demo",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    assert result.dataset.sequence_count == 9
    assert result.dataset_export.sequences_path.is_file()
    assert result.dataset_export.metadata_path.is_file()
    assert result.dataset_export.centroids_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.tree_path.is_file()
    assert result.workflow_bundle.host_switch_summary_path.is_file()
    assert result.workflow_bundle.biogeography_report_path.is_file()
    assert result.workflow_bundle.final_report_path.is_file()
    assert result.overview_path.is_file()
    assert "final report" in result.overview_path.read_text(encoding="utf-8")


def test_export_rabies_cross_host_geography_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_cross_host_geography_panel_dataset(tmp_path / "dataset")
    expected_files = {
        str(path.relative_to(result.expected_output_root))
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    }
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert result.metadata_path.is_file()
    assert result.centroids_path.is_file()
    assert "rabies-cross-host-geography-report.html" in expected_files
    assert "biogeography/biogeography-report.html" in expected_files


def test_public_runtime_exports_include_rabies_cross_host_geography_panel_surface() -> None:
    assert (
        bijux_phylogenetics.load_rabies_cross_host_geography_panel_dataset
        is load_rabies_cross_host_geography_panel_dataset
    )
    assert (
        bijux_phylogenetics.export_rabies_cross_host_geography_panel_dataset
        is export_rabies_cross_host_geography_panel_dataset
    )
    assert (
        bijux_phylogenetics.run_rabies_cross_host_geography_panel_workflow
        is run_rabies_cross_host_geography_panel_workflow
    )
    assert (
        bijux_phylogenetics.write_rabies_cross_host_geography_panel_workflow_bundle
        is write_rabies_cross_host_geography_panel_workflow_bundle
    )
    assert bijux_phylogenetics.run_rabies_cross_host_geography_panel_demo is (
        run_rabies_cross_host_geography_panel_demo
    )


def test_cli_demo_rabies_cross_host_geography_panel_json_output_reports_integrated_workflow(
    tmp_path: Path, capsys
) -> None:
    executables = require_alignment_engine_executables()
    output = tmp_path / "rabies-integrated-demo"
    exit_code = main(
        [
            "demo",
            "rabies-cross-host-geography-panel",
            "--out",
            str(output),
            "--mafft-executable",
            executables["mafft"],
            "--trimal-executable",
            executables["trimal"],
            "--iqtree-executable",
            executables["iqtree2"],
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 19
    assert payload["metrics"]["sequence_count"] == 9
    assert payload["metrics"]["host_trait"] == "host_group"
    assert payload["metrics"]["geography_trait"] == "region_group"
    assert payload["metrics"]["selected_model"] == "TPM2u+F+G4"
    assert payload["metrics"]["minimum_support"] == 84.0
    assert payload["metrics"]["maximum_support"] == 100.0
    assert payload["metrics"]["root_host"] == "bat"
    assert payload["metrics"]["root_region"] == "north_asia"
    assert payload["metrics"]["host_switch_count"] == 2
    assert payload["metrics"]["migration_event_count"] == 4
    assert payload["metrics"]["reference_output_count"] == 27
    assert payload["data"]["dataset"]["dataset_id"] == (
        "rabies_cross_host_geography_panel"
    )
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
