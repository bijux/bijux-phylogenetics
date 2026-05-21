from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_avian_reproductive_trait_dataset,
    load_avian_reproductive_trait_dataset,
    run_avian_reproductive_trait_demo,
    run_avian_reproductive_trait_workflow,
    write_avian_reproductive_trait_workflow_bundle,
)
import bijux_phylogenetics.datasets.avian_reproductive_traits as birds_api

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_avian_reproductive_trait_dataset_exposes_packaged_bird_surface() -> None:
    dataset = load_avian_reproductive_trait_dataset()
    assert dataset.dataset_id == "avian_reproductive_traits"
    assert dataset.label == "Avian reproductive trait dataset"
    assert dataset.taxon_column == "species"
    assert dataset.taxon_count == 94
    assert dataset.tree_path.is_file()
    assert dataset.traits_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert dataset.workflow_continuous_trait == "testes_mass"
    assert dataset.workflow_pgls_predictor == "body_mass"
    assert dataset.workflow_discrete_trait == "mating_system"
    assert dataset.workflow_clade_trait == "mating_system"
    assert "multiple_paternity_percentage" in dataset.continuous_traits
    assert "development_mode" in dataset.categorical_traits


@pytest.mark.slow
def test_write_avian_reproductive_trait_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_avian_reproductive_trait_workflow()
    bundle = write_avian_reproductive_trait_workflow_bundle(
        tmp_path / "workflow", report
    )
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
        bundle.clade_summary_path.name: bundle.clade_summary_path,
        bundle.clade_rows_path.name: bundle.clade_rows_path,
        bundle.clade_exclusion_path.name: bundle.clade_exclusion_path,
    }
    assert {path.name for path in expected_root.glob("*.tsv")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_avian_reproductive_trait_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_avian_reproductive_trait_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 94
    assert result.dataset_export.tree_path.is_file()
    assert result.dataset_export.traits_path.is_file()
    assert result.workflow_bundle.summary_path.is_file()
    assert result.workflow_bundle.clade_summary_path.is_file()
    assert result.overview_path.is_file()
    assert "workflow summary" in result.overview_path.read_text(encoding="utf-8")


@pytest.mark.slow
def test_export_avian_reproductive_trait_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_avian_reproductive_trait_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*.tsv")}
    assert result.readme_path.is_file()
    assert result.tree_path.is_file()
    assert result.traits_path.is_file()
    assert len(expected_files) == 17
    assert "clade-trait-summary.tsv" in expected_files


@pytest.mark.slow
def test_public_runtime_exports_include_avian_reproductive_trait_dataset_surface() -> (
    None
):
    assert (
        birds_api.load_avian_reproductive_trait_dataset
        is load_avian_reproductive_trait_dataset
    )
    assert (
        birds_api.export_avian_reproductive_trait_dataset
        is export_avian_reproductive_trait_dataset
    )
    assert (
        birds_api.run_avian_reproductive_trait_workflow
        is run_avian_reproductive_trait_workflow
    )
    assert (
        birds_api.write_avian_reproductive_trait_workflow_bundle
        is write_avian_reproductive_trait_workflow_bundle
    )
    assert (
        birds_api.run_avian_reproductive_trait_demo is run_avian_reproductive_trait_demo
    )


@pytest.mark.slow
def test_cli_demo_avian_reproductive_traits_json_output_reports_dataset_and_workflow(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "avian-demo"
    exit_code = main(
        ["demo", "avian-reproductive-traits", "--out", str(output), "--json"]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 16
    assert payload["metrics"]["dataset_taxon_count"] == 94
    assert payload["metrics"]["reference_output_count"] == 17
    assert payload["data"]["dataset"]["dataset_id"] == "avian_reproductive_traits"
    assert payload["data"]["workflow_bundle"]["summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
