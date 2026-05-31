from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest

from bijux_phylogenetics.bayesian.posterior_execution import (
    run_bayesian_posterior_execution,
)
from bijux_phylogenetics.engines.common import (
    EngineVersionInfo,
    read_engine_version,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

pytestmark = pytest.mark.engine_contract


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_posterior_engine(
    path: Path, *, version_text: str = "posterior-engine v1.0"
) -> Path:
    return _write_executable(
        path,
        f"""#!{sys.executable}
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print({version_text!r})
    raise SystemExit(0)

input_path = Path(args[-1])
output_dir = input_path.parent
(output_dir / "posterior.log").write_text("posterior\\n-100.0\\n", encoding="utf-8")
(output_dir / "posterior.trees").write_text("(A:1,B:1);\\n", encoding="utf-8")
print("warning: posterior fixture completed", file=sys.stderr)
""",
    )


def _build_execution(
    tmp_path: Path,
    *,
    executable: Path,
    validate_outputs,
    resume: bool = False,
    incomplete_run_policy: str = "reject",
) -> tuple[Path, Path, Path]:
    input_path = tmp_path / "analysis.input"
    input_path.write_text("alignment\n", encoding="utf-8")
    manifest_path = tmp_path / "analysis.manifest.json"
    version = read_engine_version(
        "PosteriorEngine",
        executable,
        version_args=("--version",),
    )
    run_bayesian_posterior_execution(
        engine_name="PosteriorEngine",
        executable=str(executable),
        version=version,
        command=[str(executable), input_path.name],
        input_paths=[input_path],
        output_paths={
            "posterior_log": tmp_path / "posterior.log",
            "posterior_trees": tmp_path / "posterior.trees",
        },
        manifest_path=manifest_path,
        stdout_path=tmp_path / "analysis.stdout.log",
        stderr_path=tmp_path / "analysis.stderr.log",
        work_dir=tmp_path,
        timeout_seconds=None,
        config={"timeout_seconds": None},
        notes=["posterior outputs validated"],
        resume=resume,
        incomplete_run_policy=incomplete_run_policy,
        validate_outputs=validate_outputs,
    )
    return (
        input_path,
        manifest_path,
        manifest_path.with_suffix(".incomplete.json"),
    )


def test_run_bayesian_posterior_execution_persists_complete_manifest(
    tmp_path: Path,
) -> None:
    executable = _fake_posterior_engine(tmp_path / "posterior-engine")
    _input_path, manifest_path, marker_path = _build_execution(
        tmp_path,
        executable=executable,
        validate_outputs=lambda: None,
    )

    assert manifest_path.exists()
    assert marker_path.exists() is False


def test_run_bayesian_posterior_execution_marks_validation_failures_incomplete(
    tmp_path: Path,
) -> None:
    executable = _fake_posterior_engine(tmp_path / "posterior-engine")

    with pytest.raises(EngineWorkflowError, match="posterior outputs are malformed"):
        _build_execution(
            tmp_path,
            executable=executable,
            validate_outputs=lambda: (_ for _ in ()).throw(
                EngineWorkflowError(
                    "posterior outputs are malformed",
                    code="posterior_outputs_malformed",
                )
            ),
        )

    marker_path = tmp_path / "analysis.manifest.incomplete.json"
    assert marker_path.exists()
    payload = json.loads(marker_path.read_text(encoding="utf-8"))
    assert payload["failure_reason"] == "posterior_outputs_malformed"
    assert payload["missing_output_names"] == []


def test_run_bayesian_posterior_execution_reuses_only_verified_completed_manifest(
    tmp_path: Path,
) -> None:
    executable = _fake_posterior_engine(tmp_path / "posterior-engine")
    _input_path, manifest_path, marker_path = _build_execution(
        tmp_path,
        executable=executable,
        validate_outputs=lambda: None,
    )

    report = run_bayesian_posterior_execution(
        engine_name="PosteriorEngine",
        executable=str(executable),
        version=EngineVersionInfo(
            engine_name="PosteriorEngine",
            executable=str(executable),
            command=[str(executable), "--version"],
            text="posterior-engine v1.0",
        ),
        command=[str(executable), "analysis.input"],
        input_paths=[tmp_path / "analysis.input"],
        output_paths={
            "posterior_log": tmp_path / "posterior.log",
            "posterior_trees": tmp_path / "posterior.trees",
        },
        manifest_path=manifest_path,
        stdout_path=tmp_path / "analysis.stdout.log",
        stderr_path=tmp_path / "analysis.stderr.log",
        work_dir=tmp_path,
        timeout_seconds=None,
        config={"timeout_seconds": None},
        notes=["posterior outputs validated"],
        resume=True,
        incomplete_run_policy="reject",
        validate_outputs=lambda: None,
    )

    assert report.resumed is True
    assert marker_path.exists() is False
