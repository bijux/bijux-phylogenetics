from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.validation import (
    AlignmentValidationMatrixInputs,
    BayesianValidationMatrixInputs,
    run_alignment_engine_validation_matrix,
    run_bayesian_engine_validation_matrix,
    run_external_engine_validation_matrix,
)

from .support.fake_bayesian_engines import (
    fake_beast,
    fake_mrbayes,
)
from .support.fake_external_engines import (
    fake_fasttree,
    fake_iqtree,
    fake_mafft,
    fake_trimal,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def _alignment_inputs() -> AlignmentValidationMatrixInputs:
    return AlignmentValidationMatrixInputs(
        raw_sequence_path=fixture("alignments/example_sequences_raw.fasta"),
        trimming_alignment_path=fixture("alignments/example_alignment_trim.fasta"),
        inference_alignment_path=fixture("alignments/example_alignment.fasta"),
    )


def _bayesian_inputs() -> BayesianValidationMatrixInputs:
    return BayesianValidationMatrixInputs(
        mrbayes_alignment_path=fixture("alignments/example_multilocus_alignment.fasta"),
        mrbayes_partition_path=fixture("alignments/example_multilocus_partitions.txt"),
        beast_alignment_path=fixture("example_alignment.fasta"),
    )


def test_run_alignment_engine_validation_matrix_collects_governed_cases(
    tmp_path: Path,
) -> None:
    matrix = run_alignment_engine_validation_matrix(
        inputs=_alignment_inputs(),
        out_dir=tmp_path / "alignment-matrix",
        mafft_executable=fake_mafft(tmp_path / "mafft"),
        trimal_executable=fake_trimal(tmp_path / "trimal"),
        iqtree_executable=fake_iqtree(tmp_path / "iqtree"),
        fasttree_executable=fake_fasttree(tmp_path / "fasttree"),
    )

    assert len(matrix.cases) == 5
    assert [case.engine_name for case in matrix.cases] == [
        "MAFFT",
        "trimAl",
        "IQ-TREE",
        "IQ-TREE",
        "FastTree",
    ]
    assert all(case.output_checksums for case in matrix.cases)


def test_run_bayesian_engine_validation_matrix_supports_live_and_governed_beast(
    tmp_path: Path,
) -> None:
    fixture_matrix = run_bayesian_engine_validation_matrix(
        inputs=_bayesian_inputs(),
        out_dir=tmp_path / "bayesian-fixture",
        mrbayes_executable=fake_mrbayes(tmp_path / "mrbayes-fixture"),
        beast_executable=None,
    )
    live_matrix = run_bayesian_engine_validation_matrix(
        inputs=_bayesian_inputs(),
        out_dir=tmp_path / "bayesian-live",
        mrbayes_executable=fake_mrbayes(tmp_path / "mrbayes-live"),
        beast_executable=fake_beast(tmp_path / "beast-live"),
    )

    assert [case.validation_mode for case in fixture_matrix.cases] == [
        "workflow-run",
        "fixture-parse",
    ]
    assert [case.validation_mode for case in live_matrix.cases] == [
        "workflow-run",
        "workflow-run",
    ]
    assert live_matrix.cases[-1].engine_name == "BEAST"


def test_run_external_engine_validation_matrix_merges_alignment_and_bayesian_cases(
    tmp_path: Path,
) -> None:
    matrix = run_external_engine_validation_matrix(
        alignment_inputs=_alignment_inputs(),
        bayesian_inputs=_bayesian_inputs(),
        out_dir=tmp_path / "external-matrix",
        mafft_executable=fake_mafft(tmp_path / "mafft"),
        trimal_executable=fake_trimal(tmp_path / "trimal"),
        iqtree_executable=fake_iqtree(tmp_path / "iqtree"),
        fasttree_executable=fake_fasttree(tmp_path / "fasttree"),
        mrbayes_executable=fake_mrbayes(tmp_path / "mrbayes"),
        beast_executable=fake_beast(tmp_path / "beast"),
    )

    assert len(matrix.cases) == 7
    assert matrix.cases[0].engine_name == "MAFFT"
    assert matrix.cases[-1].engine_name == "BEAST"
    assert all(
        case.command or case.validation_mode == "fixture-parse" for case in matrix.cases
    )
