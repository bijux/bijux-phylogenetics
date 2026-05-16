from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.api import (
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


@pytest.mark.slow
def test_python_workflow_results_write_json_and_tsv_for_major_workflows(
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

    alignment_path = tmp_path / "aligned.fasta"
    trimmed_path = tmp_path / "trimmed.fasta"

    results = {
        "validation": run_fasta_validation_workflow(input_path),
        "alignment": run_alignment_workflow(
            input_path,
            alignment_path,
            executable=mafft_executable,
        ),
        "trimming": run_trimming_workflow(
            alignment_path,
            trimmed_path,
            executable=trimal_executable,
        ),
        "inference": run_tree_inference_workflow(
            trimmed_path,
            out_dir=tmp_path / "inference",
            model="GTR+G",
            executable=iqtree_executable,
            sequence_type="dna",
        ),
        "support": run_support_workflow(
            trimmed_path,
            out_dir=tmp_path / "support",
            model="GTR+G",
            executable=iqtree_executable,
            sequence_type="dna",
            replicates=1000,
            seed=17,
        ),
        "sequence_to_tree": run_sequence_to_tree_workflow(
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
        ),
        "tree_comparison": run_tree_comparison_workflow(
            fixture("trees/example_tree.nwk"),
            fixture("trees/example_tree_alt.nwk"),
        ),
        "comparative": run_comparative_model_workflow(
            fixture("trees/example_tree.nwk"),
            fixture("metadata/example_traits_comparative.tsv"),
            response="response",
            predictors=["predictor_one"],
            lambda_value=1.0,
        ),
        "ancestral": run_ancestral_reconstruction_workflow(
            fixture("trees/example_tree.nwk"),
            fixture("metadata/example_traits_comparative.tsv"),
            trait="habitat",
            model="equal-rates",
        ),
        "report": render_report_workflow(
            tree_path=fixture("trees/example_tree.nwk"),
            alignment_path=fixture("alignments/example_alignment.fasta"),
            metadata_path=fixture("metadata/example_metadata.tsv"),
            traits_path=fixture("metadata/example_traits_comparative.tsv"),
            out_path=tmp_path / "phylogenetics-report.html",
        ),
        "configured": run_configured_phylo_workflow(config_path),
    }

    for label, result in results.items():
        json_path = tmp_path / f"{label}.json"
        payload_path = result.write_json(json_path)
        assert payload_path == json_path
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        assert payload["result_type"].endswith("Result")
        assert "report" in payload
        if hasattr(result, "write_tsv"):
            tsv_path = tmp_path / f"{label}.tsv"
            written_tsv = result.write_tsv(tsv_path)
            assert written_tsv == tsv_path
            lines = tsv_path.read_text(encoding="utf-8").splitlines()
            assert lines
            assert len(lines) >= 2


def test_python_workflow_results_emit_distinct_tsv_shapes_for_serialized_domains(
    tmp_path: Path,
) -> None:
    support_result = run_tree_comparison_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("trees/example_tree_alt.nwk"),
    )
    comparative_result = run_comparative_model_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("metadata/example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=1.0,
    )
    ancestral_result = run_ancestral_reconstruction_workflow(
        fixture("trees/example_tree.nwk"),
        fixture("metadata/example_traits_comparative.tsv"),
        trait="habitat",
        model="equal-rates",
    )

    comparison_tsv = support_result.write_tsv(tmp_path / "tree-comparison.tsv")
    comparative_tsv = comparative_result.write_tsv(tmp_path / "comparative.tsv")
    ancestral_tsv = ancestral_result.write_tsv(tmp_path / "ancestral.tsv")

    assert "robinson_foulds_distance" in comparison_tsv.read_text(encoding="utf-8")
    assert "coefficient" in comparative_tsv.read_text(encoding="utf-8")
    assert "node\t" in ancestral_tsv.read_text(encoding="utf-8")
