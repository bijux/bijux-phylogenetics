from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
)
from bijux_phylogenetics.engines import run_model_selection
from bijux_phylogenetics.engines.validation.matrix import (
    ExternalEngineValidationCase,
    ExternalEngineValidationMatrixReport,
    build_external_engine_validation_case,
    build_external_engine_validation_matrix,
    build_governed_beast_fixture_validation_case,
    merge_external_engine_validation_matrices,
    write_external_engine_validation_matrix,
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


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_iqtree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".iqtree").write_text(
    " No. Model         -LnL         df  AIC          AICc         BIC\\n"
    "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
    "Best-fit model according to BIC: GTR+G\\n"
    "Log-likelihood of the tree: -123.456\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture model-selection log\\nBEST SCORE FOUND : -123.456\\n",
    encoding="utf-8",
)
prefix.with_suffix(".model").write_text("Best-fit model: GTR+G\\n", encoding="utf-8")
raise SystemExit(0)
""",
    )


def test_build_external_engine_validation_case_from_workflow_records_run_metadata(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )

    case = build_external_engine_validation_case(
        "iqtree model selection",
        workflow,
    )
    matrix = build_external_engine_validation_matrix([case])
    output_path = write_external_engine_validation_matrix(
        tmp_path / "engine-validation-matrix.json",
        matrix,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert case.engine_name == "IQ-TREE"
    assert case.validation_mode == "workflow-run"
    assert case.manifest_path == workflow.manifest_path
    assert case.executable == str(executable)
    assert case.version_text == "IQ-TREE multicore version 2.9.9"
    assert case.command[0] == str(executable)
    assert case.exit_code == 0
    assert case.runtime_seconds is not None
    assert case.runtime_seconds >= 0.0
    assert case.output_checksums
    assert sorted(case.output_checksums) == sorted(case.output_paths)
    assert payload["case_count"] == 1
    assert payload["engine_names"] == ["IQ-TREE"]
    assert payload["cases"][0]["manifest_path"] == str(workflow.manifest_path)


def test_build_governed_beast_fixture_validation_case_parses_real_fixture_outputs(
    tmp_path: Path,
) -> None:
    fixture_entry = get_shared_beast_posterior_fixture("strict_yule_real_posterior")
    case = build_governed_beast_fixture_validation_case(
        "beast fixture parser acceptance",
        fixture_entry,
    )
    matrix = build_external_engine_validation_matrix([case])
    output_path = write_external_engine_validation_matrix(
        tmp_path / "beast-validation-matrix.json",
        matrix,
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))

    assert case.engine_name == "BEAST"
    assert case.validation_mode == "fixture-parse"
    assert case.manifest_path is None
    assert case.executable is None
    assert case.version_text == "2.7"
    assert case.command == []
    assert case.exit_code is None
    assert case.runtime_seconds is None
    assert sorted(case.output_paths) == [
        "analysis_xml",
        "posterior_log",
        "posterior_trees",
    ]
    assert len(case.output_checksums) == 3
    assert any("kept rows" in note for note in case.notes)
    assert any("kept after burn-in" in note for note in case.notes)
    assert any("consensus and maximum clade credibility" in note for note in case.notes)
    assert any("ESS and interval summaries" in note for note in case.notes)
    assert payload["engine_names"] == ["BEAST"]
    assert payload["cases"][0]["output_paths"]["analysis_xml"].endswith(
        "beast2_strict_yule_posterior.xml"
    )


def test_merge_external_engine_validation_matrices_preserves_case_order() -> None:
    first = ExternalEngineValidationMatrixReport(
        generated_at_utc="2026-05-16T00:00:00Z",
        cases=[
            ExternalEngineValidationCase(
                engine_name="MAFFT",
                validation_name="mafft alignment",
                validation_mode="workflow-run",
                manifest_path=None,
                executable="/tmp/mafft",
                version_text="MAFFT v7",
                command=["/tmp/mafft", "--auto", "input.fasta"],
                exit_code=0,
                runtime_seconds=1.0,
                output_paths={},
                output_checksums={},
            )
        ],
    )
    second = ExternalEngineValidationMatrixReport(
        generated_at_utc="2026-05-16T00:00:01Z",
        cases=[
            ExternalEngineValidationCase(
                engine_name="BEAST",
                validation_name="beast fixture parser acceptance",
                validation_mode="fixture-parse",
                manifest_path=None,
                executable=None,
                version_text="2.7",
                command=[],
                exit_code=None,
                runtime_seconds=None,
                output_paths={},
                output_checksums={},
            )
        ],
    )

    merged = merge_external_engine_validation_matrices([first, second])

    assert len(merged.cases) == 2
    assert [case.engine_name for case in merged.cases] == ["MAFFT", "BEAST"]
