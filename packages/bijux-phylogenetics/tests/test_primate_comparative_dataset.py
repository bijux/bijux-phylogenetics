from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_primate_comparative_dataset,
    load_primate_comparative_dataset,
    run_primate_comparative_demo,
    run_primate_comparative_workflow,
    write_primate_comparative_workflow_bundle,
)
import bijux_phylogenetics.datasets.primate_comparative as mammals_api

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_primate_comparative_dataset_exposes_packaged_mammal_surface() -> None:
    dataset = load_primate_comparative_dataset()
    assert dataset.dataset_id == "primate_comparative"
    assert dataset.label == "Primate comparative mammal dataset"
    assert dataset.taxon_column == "species"
    assert dataset.taxon_count == 75
    assert dataset.tree_path.is_file()
    assert dataset.traits_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert dataset.workflow_continuous_trait == "longevity"
    assert dataset.workflow_pgls_predictor == "social_group_size"
    assert dataset.workflow_discrete_trait == "mating_system"
    assert "body_mass" in dataset.continuous_traits
    assert "mating_system" in dataset.categorical_traits


@pytest.mark.slow
def test_write_primate_comparative_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_primate_comparative_workflow()
    bundle = write_primate_comparative_workflow_bundle(tmp_path / "workflow", report)
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.summary_path.name: bundle.summary_path,
        bundle.pgls_lambda_profile_path.name: bundle.pgls_lambda_profile_path,
        bundle.brownian_summary_path.name: bundle.brownian_summary_path,
        bundle.brownian_exclusion_path.name: bundle.brownian_exclusion_path,
        bundle.ou_summary_path.name: bundle.ou_summary_path,
        bundle.ou_exclusion_path.name: bundle.ou_exclusion_path,
        bundle.signal_summary_path.name: bundle.signal_summary_path,
        bundle.signal_permutations_path.name: bundle.signal_permutations_path,
        bundle.continuous_ancestral_summary_path.name: (
            bundle.continuous_ancestral_summary_path
        ),
        bundle.continuous_ancestral_uncertainty_path.name: (
            bundle.continuous_ancestral_uncertainty_path
        ),
        bundle.continuous_ancestral_exclusion_path.name: (
            bundle.continuous_ancestral_exclusion_path
        ),
        bundle.discrete_ancestral_summary_path.name: (
            bundle.discrete_ancestral_summary_path
        ),
        bundle.discrete_ancestral_probability_path.name: (
            bundle.discrete_ancestral_probability_path
        ),
        bundle.discrete_ancestral_exclusion_path.name: (
            bundle.discrete_ancestral_exclusion_path
        ),
    }
    assert {path.name for path in expected_root.glob("*.tsv")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_primate_comparative_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_primate_comparative_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 75
    assert result.dataset_export.tree_path.is_file()
    assert result.dataset_export.traits_path.is_file()
    assert result.workflow_bundle.summary_path.is_file()
    assert result.overview_path.is_file()
    assert "workflow summary" in result.overview_path.read_text(encoding="utf-8")


def test_export_primate_comparative_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_primate_comparative_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*.tsv")}
    assert result.readme_path.is_file()
    assert result.tree_path.is_file()
    assert result.traits_path.is_file()
    assert len(expected_files) == 14
    assert "workflow-summary.tsv" in expected_files


def test_public_runtime_exports_include_primate_comparative_dataset_surface() -> None:
    assert (
        mammals_api.load_primate_comparative_dataset is load_primate_comparative_dataset
    )
    assert (
        mammals_api.export_primate_comparative_dataset
        is export_primate_comparative_dataset
    )
    assert (
        mammals_api.run_primate_comparative_workflow is run_primate_comparative_workflow
    )
    assert (
        mammals_api.write_primate_comparative_workflow_bundle
        is write_primate_comparative_workflow_bundle
    )
    assert mammals_api.run_primate_comparative_demo is run_primate_comparative_demo


@pytest.mark.slow
def test_cli_demo_primate_comparative_json_output_reports_dataset_and_workflow(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "primate-demo"
    exit_code = main(["demo", "primate-comparative", "--out", str(output), "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 14
    assert payload["metrics"]["dataset_taxon_count"] == 75
    assert payload["metrics"]["reference_output_count"] == 14
    assert payload["data"]["dataset"]["dataset_id"] == "primate_comparative"
    assert payload["data"]["workflow_bundle"]["summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
