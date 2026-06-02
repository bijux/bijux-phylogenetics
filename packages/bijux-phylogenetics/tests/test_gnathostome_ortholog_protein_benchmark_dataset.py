from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_gnathostome_ortholog_protein_benchmark_dataset,
    load_gnathostome_ortholog_protein_benchmark_dataset,
    run_gnathostome_ortholog_protein_benchmark_demo,
    run_gnathostome_ortholog_protein_benchmark_workflow,
    write_gnathostome_ortholog_protein_benchmark_workflow_bundle,
)
import bijux_phylogenetics.datasets.gnathostome_ortholog_protein_benchmark as vertebrates_api

from .support.external_engines import require_alignment_engine_executables
from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)

pytestmark = [
    pytest.mark.real_local,
    pytest.mark.evaluation,
    pytest.mark.scientific_validation,
]


def test_load_gnathostome_ortholog_protein_benchmark_dataset_exposes_packaged_surface() -> (
    None
):
    dataset = load_gnathostome_ortholog_protein_benchmark_dataset()
    assert dataset.dataset_id == "gnathostome_ortholog_protein_benchmark"
    assert dataset.label == "Gnathostome ortholog protein benchmark"
    assert dataset.sequence_count == 9
    assert dataset.sequence_type == "protein"
    assert dataset.workflow_prefix == "gnathostome-ortholog-protein-benchmark"
    assert dataset.iqtree_seed == 1
    assert dataset.iqtree_threads == 1
    assert dataset.bootstrap_replicates == 1000
    assert dataset.minimum_sequence_length == 142
    assert dataset.maximum_sequence_length == 183
    assert dataset.sequences_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "trimAl governed reference corpus" in dataset.source_reference


@pytest.mark.slow
def test_write_gnathostome_ortholog_protein_benchmark_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_gnathostome_ortholog_protein_benchmark_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_gnathostome_ortholog_protein_benchmark_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.summary_path.name: bundle.summary_path,
        bundle.assumptions_path.name: bundle.assumptions_path,
        bundle.alignment_path.name: bundle.alignment_path,
        bundle.trimmed_alignment_path.name: bundle.trimmed_alignment_path,
        bundle.tree_path.name: bundle.tree_path,
        bundle.model_table_path.name: bundle.model_table_path,
        bundle.support_table_path.name: bundle.support_table_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_gnathostome_ortholog_protein_benchmark_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    result = run_gnathostome_ortholog_protein_benchmark_demo(
        tmp_path / "demo",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    assert result.dataset.sequence_count == 9
    assert result.dataset_export.sequences_path.is_file()
    assert result.workflow_bundle.summary_path.is_file()
    assert result.workflow_bundle.assumptions_path.is_file()
    assert result.workflow_bundle.tree_path.is_file()
    assert result.workflow_bundle.manifest_path.is_file()
    assert result.overview_path.is_file()
    assert "Protein-specific assumptions" in result.overview_path.read_text(
        encoding="utf-8"
    )


def test_export_gnathostome_ortholog_protein_benchmark_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_gnathostome_ortholog_protein_benchmark_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert len(expected_files) == 7
    assert "molecular-assumptions.tsv" in expected_files
    assert "gnathostome-ortholog-protein-benchmark.tree" in expected_files


def test_public_runtime_exports_include_gnathostome_ortholog_protein_benchmark_surface() -> (
    None
):
    assert (
        vertebrates_api.load_gnathostome_ortholog_protein_benchmark_dataset
        is load_gnathostome_ortholog_protein_benchmark_dataset
    )
    assert (
        vertebrates_api.export_gnathostome_ortholog_protein_benchmark_dataset
        is export_gnathostome_ortholog_protein_benchmark_dataset
    )
    assert (
        vertebrates_api.run_gnathostome_ortholog_protein_benchmark_workflow
        is run_gnathostome_ortholog_protein_benchmark_workflow
    )
    assert (
        vertebrates_api.write_gnathostome_ortholog_protein_benchmark_workflow_bundle
        is write_gnathostome_ortholog_protein_benchmark_workflow_bundle
    )
    assert (
        vertebrates_api.run_gnathostome_ortholog_protein_benchmark_demo
        is run_gnathostome_ortholog_protein_benchmark_demo
    )


@pytest.mark.slow
def test_cli_demo_gnathostome_ortholog_protein_benchmark_json_output_reports_dataset_and_workflow(
    tmp_path: Path, capsys
) -> None:
    executables = require_alignment_engine_executables()
    output = tmp_path / "protein-benchmark"
    exit_code = main(
        [
            "demo",
            "gnathostome-ortholog-protein-benchmark",
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
    assert payload["metrics"]["artifact_count"] == 12
    assert payload["metrics"]["sequence_count"] == 9
    assert payload["metrics"]["sequence_type"] == "protein"
    assert payload["metrics"]["selected_model"] == "Q.insect+I"
    assert payload["metrics"]["alignment_length"] == 185
    assert payload["metrics"]["trimmed_alignment_length"] == 185
    assert payload["metrics"]["minimum_support"] == 45.0
    assert payload["metrics"]["maximum_support"] == 93.0
    assert payload["metrics"]["weakly_supported_clade_count"] == 2
    assert payload["metrics"]["state_space"] == "amino-acid"
    assert payload["metrics"]["model_selection_scope"] == "protein-models-only"
    assert payload["metrics"]["reference_output_count"] == 7
    assert (
        payload["data"]["dataset"]["dataset_id"]
        == "gnathostome_ortholog_protein_benchmark"
    )
    assert payload["data"]["workflow_bundle"]["assumptions_path"] == str(
        output / "workflow" / "molecular-assumptions.tsv"
    )
