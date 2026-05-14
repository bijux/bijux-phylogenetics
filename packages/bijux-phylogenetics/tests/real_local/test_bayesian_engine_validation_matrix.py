from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.beast import (
    prepare_beast_time_tree_analysis,
    run_beast_posterior_inference,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
)
from bijux_phylogenetics.engines.validation_matrix import (
    build_beast_artifact_validation_case,
    build_external_engine_validation_case,
    build_external_engine_validation_matrix,
    write_external_engine_validation_matrix,
)

from ..support.external_engines import (
    real_beast_executable,
    require_bayesian_validation_matrix_executables,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
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


def _build_beast_validation_case(tmp_path: Path):
    executable = real_beast_executable()
    if executable is None:
        return build_beast_artifact_validation_case(
            "beast fixture parser acceptance",
            xml_path=fixture("beast2_strict_yule_posterior.xml"),
            log_path=fixture("beast2_strict_yule_posterior.log"),
            tree_path=fixture("beast2_strict_yule_posterior.trees"),
            burnin_fraction=0.1,
        )
    xml_path = tmp_path / "live-strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    report = run_beast_posterior_inference(
        xml_path,
        executable=executable,
        threads=1,
        seed=1,
    )
    return build_external_engine_validation_case(
        "beast posterior inference",
        report,
    )


def test_bayesian_engine_validation_matrix_collects_real_and_governed_cases(
    tmp_path: Path,
) -> None:
    executables = require_bayesian_validation_matrix_executables()

    nexus_path = tmp_path / "partitioned-analysis.nex"
    prepare_mrbayes_analysis(
        fixture("alignments/example_multilocus_alignment.fasta"),
        nexus_path,
        partition_path=fixture("alignments/example_multilocus_partitions.txt"),
        model="gtr",
        rates="gamma",
        ngen=20,
        samplefreq=10,
        printfreq=10,
        burnin_fraction=0.25,
    )
    mrbayes_report = run_mrbayes_posterior_inference(
        nexus_path,
        executable=executables["mrbayes"],
        resume=False,
    )
    beast_case = _build_beast_validation_case(tmp_path)

    matrix = build_external_engine_validation_matrix(
        [
            build_external_engine_validation_case(
                "mrbayes posterior inference",
                mrbayes_report,
            ),
            beast_case,
        ]
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
