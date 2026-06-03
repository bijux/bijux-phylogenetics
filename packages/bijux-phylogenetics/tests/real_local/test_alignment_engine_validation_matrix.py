from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation import (
    AlignmentValidationMatrixInputs,
    run_alignment_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import (
    require_alignment_validation_matrix_executables,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]

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


@pytest.mark.slow
def test_alignment_engine_validation_matrix_collects_real_run_metadata(
    tmp_path: Path,
) -> None:
    executables = require_alignment_validation_matrix_executables()
    matrix = run_alignment_engine_validation_matrix(
        inputs=AlignmentValidationMatrixInputs(
            raw_sequence_path=fixture("alignments/example_sequences_raw.fasta"),
            trimming_alignment_path=fixture("alignments/example_alignment_trim.fasta"),
            inference_alignment_path=fixture("alignments/example_alignment.fasta"),
        ),
        out_dir=tmp_path,
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree"],
        fasttree_executable=executables["fasttree"],
    )
    output_path = write_external_engine_validation_matrix(
        tmp_path / "alignment-engine-validation-matrix.json",
        matrix,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert payload["case_count"] == 5
    assert payload["engine_names"] == ["FastTree", "IQ-TREE", "MAFFT", "trimAl"]
    for case in payload["cases"]:
        assert case["validation_mode"] == "workflow-run"
        assert case["version_text"]
        assert case["command"]
        assert case["exit_code"] == 0
        assert case["runtime_seconds"] is not None
        assert case["runtime_seconds"] >= 0.0
        assert case["output_checksums"]
    assert any(
        case["validation_name"] == "iqtree bootstrap support"
        for case in payload["cases"]
    )
