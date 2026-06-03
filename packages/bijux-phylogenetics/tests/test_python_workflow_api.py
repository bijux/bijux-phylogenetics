from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.api import (
    AlignmentWorkflowResult,
    AncestralReconstructionWorkflowResult,
    ComparativeModelWorkflowResult,
    ConfiguredPhyloWorkflowResult,
    FastaValidationResult,
    InferenceWorkflowResult,
    ReportWorkflowResult,
    SequenceToTreeWorkflowResult,
    SupportWorkflowResult,
    TreeComparisonWorkflowResult,
    TrimmingWorkflowResult,
    render_report_workflow,
    run_alignment_workflow,
    run_ancestral_reconstruction_workflow,
    run_comparative_model_workflow,
    run_configured_phylo_workflow,
    run_fasta_validation_workflow,
    run_sequence_to_tree_workflow,
    run_support_workflow,
    run_tree_comparison_workflow,
    run_tree_inference_workflow,
    run_trimming_workflow,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

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


def _write_unsupported_iqtree(path: Path) -> Path:
    path.write_text(
        (
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "if '--version' in sys.argv:\n"
            "    print('IQ-TREE multicore version 1.6.0')\n"
            "    raise SystemExit(0)\n"
            "raise SystemExit(1)\n"
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def test_python_workflow_api_runs_fasta_validation_from_python() -> None:
    result = run_fasta_validation_workflow(
        fixture("alignments/example_sequences_invalid_input.fasta")
    )

    assert isinstance(result, FastaValidationResult)
    assert len(result.empty_sequences) == 1
    assert len(result.duplicate_identifiers) == 1
    assert result.summary.sequence_count == 4


def test_python_workflow_api_runs_alignment_trimming_inference_and_support_workflows(
    tmp_path: Path,
) -> None:
    from tests.support.fake_external_engines import (
        fake_iqtree,
        fake_mafft,
        fake_trimal,
    )

    input_path = fixture("alignments/example_sequences_raw.fasta")
    alignment_path = tmp_path / "aligned.fasta"
    trimmed_path = tmp_path / "trimmed.fasta"
    mafft_executable = fake_mafft(tmp_path / "mafft-fixture")
    trimal_executable = fake_trimal(tmp_path / "trimal-fixture")
    iqtree_executable = fake_iqtree(tmp_path / "iqtree-fixture")

    alignment_result = run_alignment_workflow(
        input_path,
        alignment_path,
        executable=mafft_executable,
    )
    trimming_result = run_trimming_workflow(
        alignment_path,
        trimmed_path,
        executable=trimal_executable,
    )
    inference_result = run_tree_inference_workflow(
        trimmed_path,
        out_dir=tmp_path / "inference",
        model="GTR+G",
        executable=iqtree_executable,
        sequence_type="dna",
    )
    support_result = run_support_workflow(
        trimmed_path,
        out_dir=tmp_path / "support",
        model="GTR+G",
        executable=iqtree_executable,
        sequence_type="dna",
        replicates=1000,
        seed=17,
    )

    assert isinstance(alignment_result, AlignmentWorkflowResult)
    assert alignment_result.workflow == "multiple-sequence-alignment"
    assert alignment_result.output_paths["alignment"].is_file()
    assert isinstance(trimming_result, TrimmingWorkflowResult)
    assert trimming_result.workflow == "alignment-trimming"
    assert trimming_result.output_paths["trimmed_alignment"].is_file()
    assert isinstance(inference_result, InferenceWorkflowResult)
    assert inference_result.workflow == "maximum-likelihood-tree"
    assert inference_result.output_paths["tree"].is_file()
    assert inference_result.selected_model == "GTR+G"
    assert isinstance(support_result, SupportWorkflowResult)
    assert support_result.workflow == "bootstrap-support"
    assert support_result.output_paths["support_table"].is_file()
    assert support_result.bootstrap_support_summary is not None


@pytest.mark.slow
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

    workflow_result = run_sequence_to_tree_workflow(
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
    configured_result = run_configured_phylo_workflow(config_path)

    assert isinstance(workflow_result, SequenceToTreeWorkflowResult)
    assert workflow_result.output_paths["tree"].is_file()
    assert workflow_result.model_rows[0].selected_model == "GTR+G"
    assert isinstance(configured_result, ConfiguredPhyloWorkflowResult)
    assert configured_result.bundle_validation.valid is True
    bundle_manifest = json.loads(
        configured_result.bundle_report.bundle_manifest_path.read_text(encoding="utf-8")
    )
    assert bundle_manifest["workflow"] == "fasta-to-tree"


def test_python_workflow_api_blocks_unsupported_iqtree_before_tree_inference(
    tmp_path: Path,
) -> None:
    iqtree_executable = _write_unsupported_iqtree(tmp_path / "iqtree-unsupported")

    with pytest.raises(
        EngineWorkflowError,
        match="workflow 'maximum-likelihood-tree' is blocked",
    ) as error:
        run_tree_inference_workflow(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "inference",
            model="GTR+G",
            executable=iqtree_executable,
            sequence_type="dna",
        )

    assert error.value.code == "engine_preflight_workflow_blocked"
    assert error.value.details["blocking_engines"] == ["IQ-TREE"]
    assert (tmp_path / "inference").exists() is False


def test_python_workflow_api_blocks_missing_iqtree_before_sequence_to_tree_starts(
    tmp_path: Path,
) -> None:
    from tests.support.fake_external_engines import fake_mafft, fake_trimal

    input_path = fixture("alignments/example_sequences_raw.fasta")
    mafft_executable = fake_mafft(tmp_path / "mafft-fixture")
    trimal_executable = fake_trimal(tmp_path / "trimal-fixture")
    missing_iqtree = tmp_path / "missing-iqtree"

    with pytest.raises(
        EngineWorkflowError,
        match="workflow 'fasta-to-tree' is blocked",
    ) as error:
        run_sequence_to_tree_workflow(
            input_path,
            out_dir=tmp_path / "workflow",
            prefix="python-surface",
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=missing_iqtree,
            sequence_type="dna",
        )

    assert error.value.code == "engine_preflight_workflow_blocked"
    assert error.value.details["blocking_engines"] == ["IQ-TREE"]
    assert (tmp_path / "workflow").exists() is False


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

    assert isinstance(comparison, TreeComparisonWorkflowResult)
    assert comparison.robinson_foulds_distance >= 0
    assert isinstance(comparative, ComparativeModelWorkflowResult)
    assert comparative.response == "response"
    assert comparative.coefficients
    assert isinstance(ancestral, AncestralReconstructionWorkflowResult)
    assert ancestral.trait == "habitat"
    assert ancestral.estimates


def test_python_workflow_api_renders_report_with_cli_grade_artifacts(
    tmp_path: Path,
) -> None:
    result = render_report_workflow(
        tree_path=fixture("trees/example_tree.nwk"),
        alignment_path=fixture("alignments/example_alignment.fasta"),
        metadata_path=fixture("metadata/example_metadata.tsv"),
        traits_path=fixture("metadata/example_traits_comparative.tsv"),
        out_path=tmp_path / "phylogenetics-report.html",
    )

    assert isinstance(result, ReportWorkflowResult)
    assert result.output_path.is_file()
    assert result.machine_manifest_path.is_file()
    assert result.title == "Bijux Phylogenetics Report"
