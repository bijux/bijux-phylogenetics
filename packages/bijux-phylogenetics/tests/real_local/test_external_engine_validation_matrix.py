from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation_matrix import (
    build_external_engine_validation_matrix,
    merge_external_engine_validation_matrices,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import real_beast_executable
from ..support.engine_validation_matrix_cases import (
    build_real_alignment_validation_cases,
    build_real_bayesian_validation_cases,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]


def test_external_engine_validation_matrix_collects_all_governed_engines(
    tmp_path: Path,
) -> None:
    alignment_matrix = build_external_engine_validation_matrix(
        build_real_alignment_validation_cases(tmp_path / "alignment")
    )
    bayesian_matrix = build_external_engine_validation_matrix(
        build_real_bayesian_validation_cases(tmp_path / "bayesian")
    )
    matrix = merge_external_engine_validation_matrices(
        [alignment_matrix, bayesian_matrix]
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
