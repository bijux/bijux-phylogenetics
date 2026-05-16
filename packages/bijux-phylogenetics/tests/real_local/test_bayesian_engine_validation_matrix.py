from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation_matrix import (
    build_external_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import (
    real_beast_executable,
)
from ..support.engine_validation_matrix_cases import (
    build_real_bayesian_validation_cases,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]


def test_bayesian_engine_validation_matrix_collects_real_and_governed_cases(
    tmp_path: Path,
) -> None:
    matrix = build_external_engine_validation_matrix(
        build_real_bayesian_validation_cases(tmp_path)
    )
    output_path = write_external_engine_validation_matrix(
        tmp_path / "bayesian-engine-validation-matrix.json",
        matrix,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert payload["case_count"] == 2
    assert payload["engine_names"] == ["BEAST", "MrBayes"]

    mrbayes_case = next(
        case for case in payload["cases"] if case["engine_name"] == "MrBayes"
    )
    assert mrbayes_case["validation_mode"] == "workflow-run"
    assert mrbayes_case["version_text"]
    assert mrbayes_case["command"]
    assert mrbayes_case["exit_code"] == 0
    assert mrbayes_case["runtime_seconds"] is not None
    assert mrbayes_case["runtime_seconds"] >= 0.0
    assert mrbayes_case["output_checksums"]

    beast_payload = next(
        case for case in payload["cases"] if case["engine_name"] == "BEAST"
    )
    if real_beast_executable() is None:
        assert beast_payload["validation_mode"] == "fixture-parse"
        assert beast_payload["version_text"] == "2.7"
        assert beast_payload["command"] == []
        assert beast_payload["exit_code"] is None
        assert beast_payload["output_checksums"]
    else:
        assert beast_payload["validation_mode"] == "workflow-run"
        assert beast_payload["version_text"]
        assert beast_payload["command"]
        assert beast_payload["exit_code"] == 0
        assert beast_payload["runtime_seconds"] is not None
