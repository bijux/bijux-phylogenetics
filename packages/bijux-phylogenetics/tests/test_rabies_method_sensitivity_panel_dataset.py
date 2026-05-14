from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    export_rabies_method_sensitivity_panel_dataset,
    load_rabies_method_sensitivity_panel_dataset,
    run_rabies_method_sensitivity_panel_demo,
    run_rabies_method_sensitivity_panel_workflow,
    write_rabies_method_sensitivity_panel_workflow_bundle,
)


def test_load_rabies_method_sensitivity_panel_dataset_exposes_packaged_surface() -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    assert dataset.dataset_id == "rabies_method_sensitivity_panel"
    assert dataset.label == "Rabies method-sensitivity panel"
    assert dataset.taxon_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.outgroup_taxa == ("bat_chile_rv108",)
    assert dataset.parallel_workers == 2
    assert len(dataset.variants) == 4
    assert dataset.variants[0].variant_id == "auto-gap-threshold"
    assert dataset.sequences_path.is_file()
    assert dataset.metadata_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "MG458305" in dataset.source_accessions


def test_export_rabies_method_sensitivity_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_method_sensitivity_panel_dataset(tmp_path / "dataset")
    expected_files = {
        path.relative_to(result.expected_output_root)
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    }
    assert result.readme_path.is_file()
    assert result.config_path.is_file()
    assert result.sequences_path.is_file()
    assert result.metadata_path.is_file()
    assert len(expected_files) == 70
    assert Path("parallel-execution-summary.tsv") in expected_files
    assert Path("rabies-method-sensitivity.manifest.json") in expected_files
    assert Path("parallel-logs/auto-gap-threshold.log") in expected_files
    assert Path("workflow-summary.tsv") in expected_files
    assert Path("variants/auto-gap-threshold/unrooted-conclusions.tsv") in expected_files


@pytest.mark.slow
def test_write_rabies_method_sensitivity_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_rabies_method_sensitivity_panel_workflow(tmp_path / "run")
    bundle = write_rabies_method_sensitivity_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        path.relative_to(bundle.output_root): path
        for path in bundle.output_root.rglob("*")
        if path.is_file()
    }
    expected = {
        path.relative_to(expected_root): path
        for path in expected_root.rglob("*")
        if path.is_file()
    }
    assert set(generated) == set(expected)
    for relative_path, generated_path in generated.items():
        assert generated_path.read_text(encoding="utf-8") == expected[
            relative_path
        ].read_text(encoding="utf-8")


@pytest.mark.slow
def test_run_rabies_method_sensitivity_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_rabies_method_sensitivity_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 9
    assert result.dataset_export.config_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.variant_summary_path.is_file()
    assert result.workflow_bundle.parallel_summary_path.is_file()
    assert result.workflow_bundle.manifest_path.is_file()
    assert result.workflow_bundle.task_logs_root.is_dir()
    assert result.workflow_bundle.report_path.is_file()
    assert result.overview_path.is_file()
    assert "variants" in result.overview_path.read_text(encoding="utf-8")


def test_public_runtime_exports_include_rabies_method_sensitivity_surface() -> None:
    assert (
        bijux_phylogenetics.load_rabies_method_sensitivity_panel_dataset
        is load_rabies_method_sensitivity_panel_dataset
    )
    assert (
        bijux_phylogenetics.export_rabies_method_sensitivity_panel_dataset
        is export_rabies_method_sensitivity_panel_dataset
    )
    assert (
        bijux_phylogenetics.run_rabies_method_sensitivity_panel_workflow
        is run_rabies_method_sensitivity_panel_workflow
    )
    assert (
        bijux_phylogenetics.write_rabies_method_sensitivity_panel_workflow_bundle
        is write_rabies_method_sensitivity_panel_workflow_bundle
    )
    assert bijux_phylogenetics.run_rabies_method_sensitivity_panel_demo is (
        run_rabies_method_sensitivity_panel_demo
    )


@pytest.mark.slow
def test_cli_demo_rabies_method_sensitivity_panel_json_output_reports_method_review(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "rabies-demo"
    exit_code = main(
        [
            "demo",
            "rabies-method-sensitivity-panel",
            "--out",
            str(output),
            "--mafft-executable",
            "mafft",
            "--trimal-executable",
            "trimal",
            "--iqtree-executable",
            "iqtree2",
            "--fasttree-executable",
            "FastTree",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 15
    assert payload["metrics"]["taxon_count"] == 9
    assert payload["metrics"]["variant_count"] == 4
    assert payload["metrics"]["parallel_workers"] == 2
    assert payload["metrics"]["execution_mode"] == "parallel"
    assert payload["metrics"]["stable_clade_count"] == 2
    assert payload["metrics"]["changed_clade_count"] == 8
    assert payload["metrics"]["preprocessing_change_pair_count"] == 0
    assert payload["metrics"]["rooted_engine_change_variant_count"] == 0
    assert payload["metrics"]["serious_conflict_variant_count"] == 4
    assert payload["metrics"]["reference_output_count"] == 70
    assert payload["data"]["dataset"]["dataset_id"] == "rabies_method_sensitivity_panel"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
    assert payload["data"]["workflow_bundle"]["parallel_summary_path"] == str(
        output / "workflow" / "parallel-execution-summary.tsv"
    )
