from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.api import (
    DiscreteAncestralReport,
    EngineWorkflowReport,
    FastaInputValidationReport,
    FastaToTreeWorkflowReport,
    PGLSResult,
    ReportBuildResult,
    TreeComparisonReport,
    WorkflowConfigRunReport,
    render_report_workflow,
    run_alignment_workflow,
    run_ancestral_reconstruction_workflow,
    run_comparative_model_workflow,
    run_configured_phylo_workflow,
    run_fasta_validation_workflow,
    run_sequence_to_tree_workflow,
    run_tree_comparison_workflow,
    run_tree_inference_workflow,
)

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_workflow_config(
    path: Path,
    *,
    input_path: Path,
    mafft_executable: Path,
    trimal_executable: Path,
    iqtree_executable: Path,
) -> Path:
    path.write_text(
        (
            "workflow: fasta-to-tree\n"
            "inputs:\n"
            f"  fasta: {input_path.name}\n"
            "engines:\n"
            f"  mafft_executable: {mafft_executable}\n"
            f"  trimal_executable: {trimal_executable}\n"
            f"  iqtree_executable: {iqtree_executable}\n"
            "alignment:\n"
            "  mode: auto\n"
            "  sequence_type: dna\n"
            "trimming:\n"
            "  mode: gap-threshold\n"
            "  gap_threshold: 0.2\n"
            "inference:\n"
            "  bootstrap_replicates: 1000\n"
            "  seed: 17\n"
            "  threads: 1\n"
            "output:\n"
            "  out_dir: artifacts/run\n"
            "  prefix: notebook-surface\n"
            "resources:\n"
            "  timeout_seconds: 45\n"
            "  resume: false\n"
            "  incomplete_run_policy: reject\n"
        ),
        encoding="utf-8",
    )
    return path


def test_python_workflow_api_runs_fasta_validation_from_python() -> None:
    report = run_fasta_validation_workflow(
        fixture("alignments/example_sequences_invalid_input.fasta")
    )

    assert isinstance(report, FastaInputValidationReport)
    assert len(report.empty_sequences) == 1
    assert len(report.duplicate_identifiers) == 1
    assert report.summary.sequence_count == 4


def test_python_workflow_api_runs_alignment_and_inference_workflows(
    tmp_path: Path,
) -> None:
    from tests.support.fake_external_engines import fake_iqtree, fake_mafft

    input_path = fixture("alignments/example_sequences_raw.fasta")
    alignment_path = tmp_path / "aligned.fasta"
    mafft_executable = fake_mafft(tmp_path / "mafft-fixture")
    iqtree_executable = fake_iqtree(tmp_path / "iqtree-fixture")

    alignment_report = run_alignment_workflow(
        input_path,
        alignment_path,
        executable=mafft_executable,
    )
    inference_report = run_tree_inference_workflow(
        alignment_path,
        out_dir=tmp_path / "inference",
        model="GTR+G",
        executable=iqtree_executable,
        sequence_type="dna",
    )

    assert isinstance(alignment_report, EngineWorkflowReport)
    assert alignment_report.workflow == "multiple-sequence-alignment"
    assert alignment_report.output_paths["alignment"].is_file()
    assert isinstance(inference_report, EngineWorkflowReport)
    assert inference_report.workflow == "maximum-likelihood-tree"
    assert inference_report.output_paths["tree"].is_file()
    assert inference_report.selected_model == "GTR+G"


def test_python_workflow_api_runs_sequence_to_tree_and_configured_workflows(
    tmp_path: Path,
) -> None:
    from tests.support.fake_external_engines import (
        fake_iqtree,
        fake_mafft,
        fake_trimal,
    )

    input_path = tmp_path / "input.fasta"
    input_path.write_text(
        ">A\nACTGACTG\n>B\nACTGACTA\n>C\nACTGACTC\n>D\nACTGACTT\n",
        encoding="utf-8",
    )
    mafft_executable = fake_mafft(tmp_path / "mafft-fixture")
    trimal_executable = fake_trimal(tmp_path / "trimal-fixture")
    iqtree_executable = fake_iqtree(tmp_path / "iqtree-fixture")
    config_path = _write_workflow_config(
        tmp_path / "workflow-config.yaml",
        input_path=input_path,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
    )

    workflow_report = run_sequence_to_tree_workflow(
        input_path,
        out_dir=tmp_path / "workflow",
        prefix="python-surface",
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        sequence_type="dna",
        bootstrap_replicates=1000,
        seed=17,
        threads=1,
    )
    configured_report = run_configured_phylo_workflow(config_path)

    assert isinstance(workflow_report, FastaToTreeWorkflowReport)
    assert workflow_report.output_paths["tree"].is_file()
    assert workflow_report.model_rows[0].selected_model == "GTR+G"
    assert isinstance(configured_report, WorkflowConfigRunReport)
    assert configured_report.bundle_validation.valid is True
    bundle_manifest = json.loads(
        configured_report.bundle_report.bundle_manifest_path.read_text(
            encoding="utf-8"
        )
    )
    assert bundle_manifest["workflow"] == "fasta-to-tree"


def test_python_workflow_api_runs_tree_comparative_and_ancestral_workflows() -> None:
    comparison = run_tree_comparison_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("trees/example_tree_alt.nwk"),
    )
    comparative = run_comparative_model_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("metadata/example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    ancestral = run_ancestral_reconstruction_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("metadata/example_traits_comparative.tsv"),
        trait="habitat",
        model="equal-rates",
    )

    assert isinstance(comparison, TreeComparisonReport)
    assert comparison.robinson_foulds_distance >= 0
    assert isinstance(comparative, PGLSResult)
    assert comparative.response == "response"
    assert comparative.coefficients
    assert isinstance(ancestral, DiscreteAncestralReport)
    assert ancestral.trait == "habitat"
    assert ancestral.estimates


def test_python_workflow_api_renders_report_with_cli_grade_artifacts(
    tmp_path: Path,
) -> None:
    report = render_report_workflow(
        tree_path=fixture("trees/example_tree.nwk"),
        alignment_path=fixture("alignments/example_alignment.fasta"),
        metadata_path=fixture("metadata/example_metadata.tsv"),
        traits_path=fixture("metadata/example_traits_comparative.tsv"),
        out_path=tmp_path / "phylogenetics-report.html",
    )

    assert isinstance(report, ReportBuildResult)
    assert report.output_path.is_file()
    assert report.machine_manifest_path.is_file()
    assert report.title == "Bijux Phylogenetics Report"
