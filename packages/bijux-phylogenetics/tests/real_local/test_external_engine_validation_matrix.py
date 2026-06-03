from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation import (
    AlignmentValidationMatrixInputs,
    BayesianValidationMatrixInputs,
    run_external_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import (
    real_beast_executable,
    require_alignment_validation_matrix_executables,
    require_bayesian_validation_matrix_executables,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real, pytest.mark.slow]

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in ("trees", "alignments", "metadata", "expected"):
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_external_engine_validation_matrix_collects_all_governed_engines(
    tmp_path: Path,
) -> None:
    alignment_executables = require_alignment_validation_matrix_executables()
    bayesian_executables = require_bayesian_validation_matrix_executables()
    matrix = run_external_engine_validation_matrix(
        alignment_inputs=AlignmentValidationMatrixInputs(
            raw_sequence_path=fixture("alignments/example_sequences_raw.fasta"),
            trimming_alignment_path=fixture("alignments/example_alignment_trim.fasta"),
            inference_alignment_path=fixture("alignments/example_alignment.fasta"),
        ),
        bayesian_inputs=BayesianValidationMatrixInputs(
            mrbayes_alignment_path=fixture(
                "alignments/example_multilocus_alignment.fasta"
            ),
            mrbayes_partition_path=fixture(
                "alignments/example_multilocus_partitions.txt"
            ),
            beast_alignment_path=fixture("example_alignment.fasta"),
        ),
        out_dir=tmp_path,
        mafft_executable=alignment_executables["mafft"],
        trimal_executable=alignment_executables["trimal"],
        iqtree_executable=alignment_executables["iqtree"],
        fasttree_executable=alignment_executables["fasttree"],
        mrbayes_executable=bayesian_executables["mrbayes"],
        beast_executable=real_beast_executable(),
    )

    output_path = write_external_engine_validation_matrix(
        tmp_path / "external-engine-validation-matrix.json",
        matrix,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert payload["case_count"] == 7
    assert payload["engine_names"] == [
        "BEAST",
        "FastTree",
        "IQ-TREE",
        "MAFFT",
        "MrBayes",
        "trimAl",
    ]
    assert payload["cases"][0]["engine_name"] == "MAFFT"
    assert payload["cases"][-1]["engine_name"] == "BEAST"

    for case in payload["cases"]:
        assert case["output_checksums"]
        assert case["output_paths"]
        if case["engine_name"] == "BEAST" and real_beast_executable() is None:
            assert case["validation_mode"] == "fixture-parse"
            assert case["command"] == []
            assert case["exit_code"] is None
            assert case["runtime_seconds"] is None
        else:
            assert case["validation_mode"] == "workflow-run"
            assert case["version_text"]
            assert case["command"]
            assert case["exit_code"] == 0
            assert case["runtime_seconds"] is not None
            assert case["runtime_seconds"] >= 0.0
