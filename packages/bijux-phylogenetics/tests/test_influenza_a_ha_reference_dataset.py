from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_influenza_a_ha_reference_dataset,
    load_influenza_a_ha_reference_dataset,
    run_influenza_a_ha_reference_demo,
    run_influenza_a_ha_reference_workflow,
    write_influenza_a_ha_reference_workflow_bundle,
)
import bijux_phylogenetics.datasets.influenza_a_ha_reference as viruses_api

from .support.external_engines import require_alignment_engine_executables
from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)

pytestmark = [
    pytest.mark.real_local,
    pytest.mark.evaluation,
    pytest.mark.scientific_validation,
]


def test_load_influenza_a_ha_reference_dataset_exposes_packaged_viral_surface() -> None:
    dataset = load_influenza_a_ha_reference_dataset()
    assert dataset.dataset_id == "influenza_a_ha_reference_panel"
    assert dataset.label == "Influenza A hemagglutinin reference panel"
    assert dataset.sequence_count == 6
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_prefix == "influenza-a-ha-reference-panel"
    assert dataset.iqtree_seed == 1
    assert dataset.iqtree_threads == 1
    assert dataset.bootstrap_replicates == 1000
    assert dataset.sequences_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "NC_002017.1" in dataset.source_accessions


@pytest.mark.slow
def test_write_influenza_a_ha_reference_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_influenza_a_ha_reference_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_influenza_a_ha_reference_workflow_bundle(
        tmp_path / "workflow", report
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.summary_path.name: bundle.summary_path,
        bundle.alignment_path.name: bundle.alignment_path,
        bundle.trimmed_alignment_path.name: bundle.trimmed_alignment_path,
        bundle.tree_path.name: bundle.tree_path,
        bundle.model_table_path.name: bundle.model_table_path,
        bundle.support_table_path.name: bundle.support_table_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_influenza_a_ha_reference_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    result = run_influenza_a_ha_reference_demo(
        tmp_path / "demo",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    assert result.dataset.sequence_count == 6
    assert result.dataset_export.sequences_path.is_file()
    assert result.workflow_bundle.summary_path.is_file()
    assert result.workflow_bundle.tree_path.is_file()
    assert result.workflow_bundle.manifest_path.is_file()
    assert result.overview_path.is_file()
    assert "final supported tree" in result.overview_path.read_text(encoding="utf-8")


def test_export_influenza_a_ha_reference_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_influenza_a_ha_reference_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert len(expected_files) == 6
    assert "influenza-a-ha-reference-panel.tree" in expected_files


def test_public_runtime_exports_include_influenza_a_ha_reference_dataset_surface() -> (
    None
):
    assert (
        viruses_api.load_influenza_a_ha_reference_dataset
        is load_influenza_a_ha_reference_dataset
    )
    assert (
        viruses_api.export_influenza_a_ha_reference_dataset
        is export_influenza_a_ha_reference_dataset
    )
    assert (
        viruses_api.run_influenza_a_ha_reference_workflow
        is run_influenza_a_ha_reference_workflow
    )
    assert (
        viruses_api.write_influenza_a_ha_reference_workflow_bundle
        is write_influenza_a_ha_reference_workflow_bundle
    )
    assert (
        viruses_api.run_influenza_a_ha_reference_demo
        is run_influenza_a_ha_reference_demo
    )


@pytest.mark.slow
def test_cli_demo_influenza_a_ha_reference_panel_json_output_reports_dataset_and_workflow(
    tmp_path: Path, capsys
) -> None:
    executables = require_alignment_engine_executables()
    output = tmp_path / "viral-demo"
    exit_code = main(
        [
            "demo",
            "influenza-a-ha-reference-panel",
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
    assert payload["metrics"]["artifact_count"] == 11
    assert payload["metrics"]["sequence_count"] == 6
    assert payload["metrics"]["sequence_type"] == "dna"
    assert payload["metrics"]["selected_model"] == "TPM3u+F+R2"
    assert payload["metrics"]["minimum_support"] == 98.0
    assert payload["metrics"]["maximum_support"] == 99.0
    assert payload["metrics"]["weakly_supported_clade_count"] == 0
    assert payload["metrics"]["reference_output_count"] == 6
    assert payload["data"]["dataset"]["dataset_id"] == "influenza_a_ha_reference_panel"
    assert payload["data"]["workflow_bundle"]["summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
