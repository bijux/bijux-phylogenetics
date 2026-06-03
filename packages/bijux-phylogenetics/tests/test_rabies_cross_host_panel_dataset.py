from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_rabies_cross_host_panel_dataset,
    load_rabies_cross_host_panel_dataset,
    run_rabies_cross_host_panel_demo,
    run_rabies_cross_host_panel_workflow,
    write_rabies_cross_host_panel_workflow_bundle,
)
import bijux_phylogenetics.datasets.pathogens as pathogens_api

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_rabies_cross_host_panel_dataset_exposes_packaged_surface() -> None:
    dataset = load_rabies_cross_host_panel_dataset()
    assert dataset.dataset_id == "rabies_cross_host_panel"
    assert dataset.label == "Rabies cross-host nucleoprotein panel"
    assert dataset.taxon_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_trait == "host_group"
    assert dataset.workflow_model == "ard"
    assert dataset.observed_host_group_count == 3
    assert dataset.sequences_path.is_file()
    assert dataset.tree_path.is_file()
    assert dataset.hosts_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "MG458305" in dataset.source_accessions


@pytest.mark.slow
def test_write_rabies_cross_host_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_rabies_cross_host_panel_workflow()
    bundle = write_rabies_cross_host_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.workflow_summary_path.name: bundle.workflow_summary_path,
        bundle.host_switch_summary_path.name: bundle.host_switch_summary_path,
        bundle.host_state_nodes_path.name: bundle.host_state_nodes_path,
        bundle.host_switch_branches_path.name: bundle.host_switch_branches_path,
        bundle.host_switch_counts_path.name: bundle.host_switch_counts_path,
        bundle.host_switch_fits_path.name: bundle.host_switch_fits_path,
        bundle.host_switch_unsupported_path.name: bundle.host_switch_unsupported_path,
        bundle.host_switch_exclusions_path.name: bundle.host_switch_exclusions_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_rabies_cross_host_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_rabies_cross_host_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 9
    assert result.dataset_export.sequences_path.is_file()
    assert result.dataset_export.tree_path.is_file()
    assert result.dataset_export.hosts_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.host_switch_summary_path.is_file()
    assert result.workflow_bundle.host_state_nodes_path.is_file()
    assert result.overview_path.is_file()
    assert "host-switch summary" in result.overview_path.read_text(encoding="utf-8")


def test_export_rabies_cross_host_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_cross_host_panel_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert result.tree_path.is_file()
    assert result.hosts_path.is_file()
    assert len(expected_files) == 8
    assert "host-switch-counts.tsv" in expected_files


def test_public_runtime_exports_include_rabies_cross_host_panel_surface() -> None:
    assert (
        pathogens_api.load_rabies_cross_host_panel_dataset
        is load_rabies_cross_host_panel_dataset
    )
    assert (
        pathogens_api.export_rabies_cross_host_panel_dataset
        is export_rabies_cross_host_panel_dataset
    )
    assert (
        pathogens_api.run_rabies_cross_host_panel_workflow
        is run_rabies_cross_host_panel_workflow
    )
    assert (
        pathogens_api.write_rabies_cross_host_panel_workflow_bundle
        is write_rabies_cross_host_panel_workflow_bundle
    )
    assert (
        pathogens_api.run_rabies_cross_host_panel_demo
        is run_rabies_cross_host_panel_demo
    )


@pytest.mark.slow
def test_cli_demo_rabies_cross_host_panel_json_output_reports_host_switch_review(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "rabies-demo"
    exit_code = main(
        [
            "demo",
            "rabies-cross-host-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 13
    assert payload["metrics"]["taxon_count"] == 9
    assert payload["metrics"]["workflow_trait"] == "host_group"
    assert payload["metrics"]["observed_host_group_count"] == 3
    assert payload["metrics"]["analysis_constraint_mode"] == "unconstrained"
    assert payload["metrics"]["root_host"] == "bat"
    assert payload["metrics"]["root_confidence"] == 1.0
    assert payload["metrics"]["host_switch_count"] == 2
    assert payload["metrics"]["certain_host_switch_count"] == 0
    assert payload["metrics"]["uncertain_host_switch_count"] == 2
    assert payload["metrics"]["reference_output_count"] == 8
    assert payload["data"]["dataset"]["dataset_id"] == "rabies_cross_host_panel"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
