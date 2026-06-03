from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    export_pleistocene_bear_cytb_fragment_dataset,
    load_pleistocene_bear_cytb_fragment_dataset,
    run_pleistocene_bear_cytb_fragment_demo,
    run_pleistocene_bear_cytb_fragment_workflow,
    write_pleistocene_bear_cytb_fragment_workflow_bundle,
)
import bijux_phylogenetics.datasets.pleistocene_bear_cytb_fragments as ancient_dna_api

from .support.external_engines import require_alignment_engine_executables
from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)

pytestmark = [
    pytest.mark.real_local,
    pytest.mark.evaluation,
    pytest.mark.scientific_validation,
]


def test_load_pleistocene_bear_cytb_fragment_dataset_exposes_packaged_surface() -> None:
    dataset = load_pleistocene_bear_cytb_fragment_dataset()
    assert dataset.dataset_id == "pleistocene_bear_cytb_fragments"
    assert dataset.label == "Pleistocene bear CYTB fragment panel"
    assert dataset.sequence_count == 5
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_prefix == "pleistocene-bear-cytb-fragments"
    assert dataset.iqtree_seed == 1
    assert dataset.iqtree_threads == 1
    assert dataset.bootstrap_replicates == 1000
    assert dataset.site_missingness_threshold == 0.15
    assert dataset.sequence_missingness_threshold == 0.15
    assert dataset.degraded_sequence_ids == (
        "cave_bear_ud1838_fragment",
        "cave_bear_wk01_fragment",
    )
    assert dataset.sequences_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "KX641337.1" in dataset.source_accessions


@pytest.mark.slow
def test_write_pleistocene_bear_cytb_fragment_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_pleistocene_bear_cytb_fragment_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_pleistocene_bear_cytb_fragment_workflow_bundle(
        tmp_path / "workflow", report
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.summary_path.name: bundle.summary_path,
        bundle.missingness_effects_path.name: bundle.missingness_effects_path,
        bundle.alignment_path.name: bundle.alignment_path,
        bundle.trimmed_alignment_path.name: bundle.trimmed_alignment_path,
        bundle.cleaned_alignment_path.name: bundle.cleaned_alignment_path,
        bundle.tree_path.name: bundle.tree_path,
        bundle.model_table_path.name: bundle.model_table_path,
        bundle.support_table_path.name: bundle.support_table_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_pleistocene_bear_cytb_fragment_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    result = run_pleistocene_bear_cytb_fragment_demo(
        tmp_path / "demo",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    assert result.dataset.sequence_count == 5
    assert result.dataset_export.sequences_path.is_file()
    assert result.workflow_bundle.summary_path.is_file()
    assert result.workflow_bundle.missingness_effects_path.is_file()
    assert result.workflow_bundle.cleaned_alignment_path.is_file()
    assert result.workflow_bundle.tree_path.is_file()
    assert result.workflow_bundle.manifest_path.is_file()
    assert result.overview_path.is_file()
    assert "cleaned alignment" in result.overview_path.read_text(encoding="utf-8")


def test_export_pleistocene_bear_cytb_fragment_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_pleistocene_bear_cytb_fragment_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert len(expected_files) == 8
    assert "missingness-effects.tsv" in expected_files
    assert "pleistocene-bear-cytb-fragments.cleaned.aln" in expected_files


def test_public_runtime_exports_include_pleistocene_bear_cytb_fragment_surface() -> (
    None
):
    assert (
        ancient_dna_api.load_pleistocene_bear_cytb_fragment_dataset
        is load_pleistocene_bear_cytb_fragment_dataset
    )
    assert (
        ancient_dna_api.export_pleistocene_bear_cytb_fragment_dataset
        is export_pleistocene_bear_cytb_fragment_dataset
    )
    assert (
        ancient_dna_api.run_pleistocene_bear_cytb_fragment_workflow
        is run_pleistocene_bear_cytb_fragment_workflow
    )
    assert (
        ancient_dna_api.write_pleistocene_bear_cytb_fragment_workflow_bundle
        is write_pleistocene_bear_cytb_fragment_workflow_bundle
    )
    assert (
        ancient_dna_api.run_pleistocene_bear_cytb_fragment_demo
        is run_pleistocene_bear_cytb_fragment_demo
    )


@pytest.mark.slow
def test_cli_demo_pleistocene_bear_cytb_fragments_json_output_reports_missingness_review(
    tmp_path: Path, capsys
) -> None:
    executables = require_alignment_engine_executables()
    output = tmp_path / "ancient-dna-demo"
    exit_code = main(
        [
            "demo",
            "pleistocene-bear-cytb-fragments",
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
    assert payload["metrics"]["sequence_count"] == 5
    assert payload["metrics"]["degraded_sequence_count"] == 2
    assert payload["metrics"]["selected_model"] == "HKY+F"
    assert payload["metrics"]["minimum_support"] == 80.0
    assert payload["metrics"]["maximum_support"] == 100.0
    assert payload["metrics"]["removed_column_count"] == 166
    assert payload["metrics"]["cleaned_missing_data_fraction"] == 0.0
    assert payload["metrics"]["reference_output_count"] == 8
    assert payload["data"]["dataset"]["dataset_id"] == "pleistocene_bear_cytb_fragments"
    assert payload["data"]["workflow_bundle"]["summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
