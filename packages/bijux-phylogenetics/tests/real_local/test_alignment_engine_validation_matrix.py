from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.engines import (
    run_alignment_trimming,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.engines.validation_matrix import (
    build_external_engine_validation_case,
    build_external_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import require_alignment_validation_matrix_executables

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_alignment_engine_validation_matrix_collects_real_run_metadata(
    tmp_path: Path,
) -> None:
    executables = require_alignment_validation_matrix_executables()

    mafft_report = run_multiple_sequence_alignment(
        fixture("alignments/example_sequences_raw.fasta"),
        tmp_path / "alignment" / "real-mafft-alignment.fasta",
        executable=executables["mafft"],
        mode="linsi",
    )
    trimal_report = run_alignment_trimming(
        fixture("alignments/example_alignment_trim.fasta"),
        tmp_path / "trim" / "real-trimal-trimmed.fasta",
        executable=executables["trimal"],
        mode="gap-threshold",
        gap_threshold=0.8,
    )
    model_report = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executables["iqtree"],
        prefix="real",
        sequence_type="dna",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model=model_report.selected_model or "GTR+G",
        executable=executables["iqtree"],
        prefix="real",
        sequence_type="dna",
        replicates=1000,
    )
    fasttree_report = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "fasttree" / "real-fasttree.nwk",
        executable=executables["fasttree"],
        sequence_type="dna",
    )

    matrix = build_external_engine_validation_matrix(
        [
            build_external_engine_validation_case("mafft alignment", mafft_report),
            build_external_engine_validation_case("trimal trimming", trimal_report),
            build_external_engine_validation_case(
                "iqtree model selection",
                model_report,
                notes=[
                    f"selected model: {model_report.selected_model}",
                ],
            ),
            build_external_engine_validation_case(
                "iqtree bootstrap support",
                bootstrap_report,
                notes=[
                    f"support value count: {bootstrap_report.iqtree_summary.support_value_count if bootstrap_report.iqtree_summary is not None else 0}",
                ],
            ),
            build_external_engine_validation_case("fasttree inference", fasttree_report),
        ]
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
    assert any(case["validation_name"] == "iqtree bootstrap support" for case in payload["cases"])
