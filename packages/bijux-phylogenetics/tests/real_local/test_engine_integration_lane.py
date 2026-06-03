from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines import (
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
    run_tree_inference_comparison,
)
from bijux_phylogenetics.engines.inference import (
    run_inference_reproducibility_check,
    run_large_alignment_inference,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment

from ..support.external_engines import (
    real_fasttree_executable,
    real_iqtree_executable,
    real_mafft_executable,
    real_trimal_executable,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_run_multiple_sequence_alignment_with_real_mafft_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = real_mafft_executable()
    if executable is None:
        pytest.skip("real MAFFT executable is not available for integration coverage")

    output_path = tmp_path / "real-mafft-alignment.fasta"
    report = run_multiple_sequence_alignment(
        fixture("alignments/example_sequences_raw.fasta"),
        output_path,
        executable=executable,
        mode="linsi",
    )

    alignment = load_fasta_alignment(output_path)
    widths = {len(record.sequence) for record in alignment}

    assert report.run.exit_code == 0
    assert report.run.command[1:-1] == ["--localpair", "--maxiterate", "1000"]
    assert "v7." in report.run.version.text
    assert report.manifest_path.exists()
    assert len(alignment) == 4
    assert len(widths) == 1


def test_run_alignment_trimming_with_real_trimal_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = real_trimal_executable()
    if executable is None:
        pytest.skip("real trimAl executable is not available for integration coverage")

    input_path = fixture("alignments/example_alignment_trim.fasta")
    output_path = tmp_path / "real-trimal-trimmed.fasta"
    report = run_alignment_trimming(
        input_path,
        output_path,
        executable=executable,
        mode="gap-threshold",
        gap_threshold=0.8,
    )

    input_records = load_fasta_alignment(input_path)
    trimmed_records = load_fasta_alignment(output_path)

    assert report.run.exit_code == 0
    assert report.run.command[5:] == ["-gt", "0.800000"]
    assert "trimal v" in report.run.version.text.lower()
    assert report.manifest_path.exists()
    assert report.trimming_summary is not None
    assert report.trimming_summary.removed_site_count >= 1
    assert report.trimming_summary.trimmed_alignment_length == len(
        trimmed_records[0].sequence
    )
    assert len(trimmed_records[0].sequence) < len(input_records[0].sequence)


@pytest.mark.slow
def test_run_iqtree_backend_with_real_executable_on_small_dataset(
    tmp_path: Path,
) -> None:
    executable = real_iqtree_executable()
    if executable is None:
        pytest.skip("real IQ-TREE executable is not available for integration coverage")

    input_path = fixture("alignments/example_alignment.fasta")
    model_report = run_model_selection(
        input_path,
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="real",
        sequence_type="dna",
    )
    assert model_report.selected_model is not None
    ml_report = run_maximum_likelihood_tree_inference(
        input_path,
        out_dir=tmp_path / "ml",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        input_path,
        out_dir=tmp_path / "bootstrap",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
        replicates=1000,
    )
    sh_alrt_report = run_sh_alrt_support_estimation(
        input_path,
        out_dir=tmp_path / "sh-alrt",
        model=model_report.selected_model,
        executable=executable,
        prefix="real",
        sequence_type="dna",
        sh_alrt_replicates=1000,
        bootstrap_replicates=1000,
    )
    consensus_report = run_bootstrap_consensus_tree(
        bootstrap_report.output_paths["bootstrap_trees"],
        out_dir=tmp_path / "consensus",
        executable=executable,
        prefix="real",
    )

    assert model_report.output_paths["iqtree_report"].exists()
    assert model_report.output_paths["iqtree_log"].exists()
    assert model_report.output_paths["model_candidates"].exists()
    assert model_report.log_likelihood is not None
    assert model_report.model_selection_summary is not None
    assert model_report.model_selection_summary.candidate_count >= 1
    assert model_report.model_selection_summary.best_model_bic is not None
    assert ml_report.output_paths["tree"].exists()
    assert ml_report.output_paths["iqtree_log"].exists()
    assert ml_report.log_likelihood is not None
    assert bootstrap_report.output_paths["bootstrap_trees"].exists()
    assert bootstrap_report.output_paths["support_tree"].exists()
    assert bootstrap_report.output_paths["support_table"].exists()
    assert bootstrap_report.output_paths["low_support_branches"].exists()
    assert bootstrap_report.output_paths["support_histogram"].exists()
    assert bootstrap_report.output_paths["iqtree_log"].exists()
    assert bootstrap_report.log_likelihood is not None
    assert bootstrap_report.iqtree_summary is not None
    assert bootstrap_report.iqtree_summary.support_value_count >= 1
    assert bootstrap_report.bootstrap_support_summary is not None
    assert bootstrap_report.weak_backbone_report is not None
    assert sh_alrt_report.output_paths["support_tree"].exists()
    assert sh_alrt_report.output_paths["bootstrap_trees"].exists()
    assert sh_alrt_report.output_paths["support_table"].exists()
    assert sh_alrt_report.output_paths["conflicting_support_branches"].exists()
    assert sh_alrt_report.iqtree_summary is not None
    assert sh_alrt_report.iqtree_summary.support_value_count >= 1
    assert sh_alrt_report.sh_alrt_support_summary is not None
    assert sh_alrt_report.sh_alrt_support_summary.annotated_node_count >= 1
    assert consensus_report.output_paths["consensus_tree"].exists()
    assert consensus_report.output_paths["iqtree_log"].exists()
    assert consensus_report.iqtree_summary is not None
    assert consensus_report.iqtree_summary.support_value_count >= 1


@pytest.mark.parametrize(
    ("input_path", "sequence_type", "expected_command_prefix"),
    (
        (fixture("alignments/example_alignment.fasta"), "dna", ["-gtr", "-nt"]),
        (
            REPOSITORY_ROOT
            / "packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree/strnog-enog411bqtj-proteins/strnog-enog411bqtj-proteins.trimmed.aln",
            "protein",
            ["-lg"],
        ),
    ),
)
def test_run_fasttree_backend_with_real_executable_on_supported_alignments(
    tmp_path: Path,
    input_path: Path,
    sequence_type: str,
    expected_command_prefix: list[str],
) -> None:
    executable = real_fasttree_executable()
    if executable is None:
        pytest.skip(
            "real FastTree executable is not available for integration coverage"
        )

    output_path = tmp_path / f"real-fasttree-{sequence_type}.nwk"
    report = run_fast_tree_inference(
        input_path,
        output_path,
        executable=executable,
        sequence_type=sequence_type,
    )

    assert report.run.exit_code == 0
    assert report.run.command[1 : 1 + len(expected_command_prefix)] == (
        expected_command_prefix
    )
    assert "FastTree" in report.run.version.text
    assert report.fasttree_support_summary is not None
    assert report.fasttree_support_summary.annotated_node_count > 0
    assert report.fasttree_support_summary.approximate_method is True
    assert report.output_paths["support_table"].exists()
    assert report.output_paths["low_support_branches"].exists()
    assert report.output_paths["support_histogram"].exists()


def test_run_tree_inference_comparison_with_real_executables_on_small_alignment(
    tmp_path: Path,
) -> None:
    iqtree_executable = real_iqtree_executable()
    fasttree_executable = real_fasttree_executable()
    if iqtree_executable is None or fasttree_executable is None:
        pytest.skip(
            "real IQ-TREE and FastTree executables are required for comparison integration coverage"
        )

    report = run_tree_inference_comparison(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "real-engine-comparison",
        prefix="real",
        sequence_type="dna",
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=1,
        iqtree_threads=1,
        bootstrap_replicates=1000,
    )

    assert report.output_paths["comparison_report"].exists()
    assert report.output_paths["comparison_table"].exists()
    assert report.engine_comparison.topology.shared_taxa
    assert report.engine_comparison.support.shared_taxa


@pytest.mark.slow
def test_run_inference_reproducibility_check_with_real_iqtree_on_small_alignment(
    tmp_path: Path,
) -> None:
    executable = real_iqtree_executable()
    if executable is None:
        pytest.skip(
            "real IQ-TREE executable is required for reproducibility integration coverage"
        )

    report = run_inference_reproducibility_check(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "reproducibility",
        executable=executable,
        repeats=2,
        bootstrap_replicates=1000,
        seed=1,
        threads=1,
    )

    assert report.output_paths["comparison_table"].exists()
    assert report.output_paths["support_delta_table"].exists()
    assert report.overall_status != "unstable"
    assert all(row.classification != "unstable" for row in report.comparison_rows)


def test_run_large_alignment_inference_with_real_fasttree_on_stress_fixture(
    tmp_path: Path,
) -> None:
    executable = real_fasttree_executable()
    if executable is None:
        pytest.skip(
            "real FastTree executable is not available for stress-fixture coverage"
        )

    input_path = (
        REPOSITORY_ROOT
        / "packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree/strnog-enog411bqtj-proteins/strnog-enog411bqtj-proteins.trimmed.aln"
    )
    report = run_large_alignment_inference(
        input_path,
        out_dir=tmp_path / "large-inference",
        prefix="strnog-enog411bqtj-proteins",
        sequence_type="protein",
        executable=executable,
        timeout_seconds=120.0,
    )

    assert report.input_summary.sequence_count >= 30
    assert report.input_summary.alignment_length >= 100
    assert report.input_summary.total_site_cells >= 4000
    assert report.output_paths["tree"].exists()
    assert report.output_paths["resource_table"].exists()
    assert any(row.elapsed_seconds >= 0.0 for row in report.resource_rows)
