from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines.inference import (
    replay_workflow_manifest,
    run_large_alignment_inference,
)

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_fasttree(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

print("((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);")
print("warning: fasttree fixture approximate support only", file=sys.stderr)
""",
    )


def test_replay_workflow_manifest_reuses_recorded_command_executable(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree(tmp_path / "FastTree-fixture")
    workflow = run_large_alignment_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "large-alignment",
        prefix="example",
        sequence_type="dna",
        executable=executable,
    )

    report = replay_workflow_manifest(
        workflow.manifest_path,
        out_dir=tmp_path / "large-alignment-replay",
    )

    assert report.workflow == "large-alignment-inference"
    assert report.engine_version_drift_detected is False
    assert report.outputs_equivalent is True
    assert report.replay_manifest_path.exists()
