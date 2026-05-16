from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines.validation_matrix import (
    build_external_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.engine_validation_matrix_cases import (
    build_real_alignment_validation_cases,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]


def test_alignment_engine_validation_matrix_collects_real_run_metadata(
    tmp_path: Path,
) -> None:
    matrix = build_external_engine_validation_matrix(
        build_real_alignment_validation_cases(tmp_path)
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
