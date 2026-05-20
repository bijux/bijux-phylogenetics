from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines.inference import (
    run_large_alignment_inference,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


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


def _fake_fasttree_streaming(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

input_path = Path(args[-1])
identifiers = []
with input_path.open(encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if line.startswith(">"):
            identifiers.append(line[1:])

tips = [f"{identifier}:0.1" for identifier in identifiers]
while len(tips) > 1:
    left = tips.pop(0)
    right = tips.pop(0)
    tips.append(f"({left},{right})0.95:0.1")
if tips:
    print(tips[0] + ";")
else:
    raise SystemExit(2)
print("warning: fasttree fixture streamed identifiers only", file=sys.stderr)
""",
    )


def _fake_fasttree_slow(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys
import time

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

time.sleep(2.0)
print("((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);")
""",
    )


def _write_large_alignment(
    path: Path, *, sequence_count: int, sequence_length: int
) -> Path:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    with path.open("w", encoding="utf-8") as handle:
        for index in range(sequence_count):
            sequence = "".join(
                alphabet[(index + offset) % len(alphabet)]
                for offset in range(sequence_length)
            )
            handle.write(f">taxon_{index:04d}\n{sequence}\n")
    return path


@pytest.mark.slow
def test_run_large_alignment_inference_streams_many_sequences_and_reports_resources(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree_streaming(tmp_path / "FastTree-streaming-fixture")
    input_path = _write_large_alignment(
        tmp_path / "large-alignment.fasta",
        sequence_count=1200,
        sequence_length=24,
    )

    report = run_large_alignment_inference(
        input_path,
        out_dir=tmp_path / "large-inference",
        prefix="stress",
        sequence_type="protein",
        executable=executable,
    )

    assert report.input_summary.sequence_count == 1200
    assert report.input_summary.alignment_length == 24
    assert report.input_summary.total_site_cells == 28800
    assert [row.stage for row in report.resource_rows] == [
        "preflight-scan",
        "fasttree-inference",
    ]
    assert report.output_paths["resource_table"].exists()
    assert report.output_paths["tree"].exists()
    assert report.output_paths["log"].exists()
    assert report.workflow == "large-alignment-inference"
    assert report.runtime_seconds >= 0.0
    assert report.config["sequence_type"] == "protein"
    assert any("scanned linearly before inference" in note for note in report.notes)
    assert any(
        "does not create an intermediate alignment copy" in note
        for note in report.notes
    )


def test_run_large_alignment_inference_resume_reuses_verified_outputs(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree(tmp_path / "FastTree-fixture")
    input_path = fixture("alignments/example_alignment.fasta")

    first = run_large_alignment_inference(
        input_path,
        out_dir=tmp_path / "large-inference",
        prefix="example",
        sequence_type="dna",
        executable=executable,
        resume=False,
    )
    second = run_large_alignment_inference(
        input_path,
        out_dir=tmp_path / "large-inference",
        prefix="example",
        sequence_type="dna",
        executable=executable,
        resume=True,
    )

    assert first.resumed is False
    assert second.resumed is True
    assert second.output_checksums == first.output_checksums
    assert second.workflow == "large-alignment-inference"


def test_run_large_alignment_inference_honors_timeout_seconds(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree_slow(tmp_path / "FastTree-slow-fixture")

    with pytest.raises(
        EngineWorkflowError,
        match="timed out after 0.100 seconds",
    ):
        run_large_alignment_inference(
            fixture("alignments/example_alignment.fasta"),
            out_dir=tmp_path / "large-inference",
            prefix="timeout",
            sequence_type="dna",
            executable=executable,
            timeout_seconds=0.1,
        )


def test_run_large_alignment_inference_can_clean_incomplete_outputs(
    tmp_path: Path,
) -> None:
    slow_executable = _fake_fasttree_slow(tmp_path / "FastTree-slow-fixture")
    out_dir = tmp_path / "large-inference"
    prefix = "timeout"

    with pytest.raises(EngineWorkflowError, match="timed out after 0.100 seconds"):
        run_large_alignment_inference(
            fixture("alignments/example_alignment.fasta"),
            out_dir=out_dir,
            prefix=prefix,
            sequence_type="dna",
            executable=slow_executable,
            timeout_seconds=0.1,
        )

    marker_path = out_dir / f"{prefix}.manifest.incomplete.json"
    assert marker_path.exists()

    executable = _fake_fasttree(tmp_path / "FastTree-fixture")
    report = run_large_alignment_inference(
        fixture("alignments/example_alignment.fasta"),
        out_dir=out_dir,
        prefix=prefix,
        sequence_type="dna",
        executable=executable,
        timeout_seconds=0.1,
        incomplete_run_policy="clean",
    )

    assert report.output_paths["tree"].exists()
    assert marker_path.exists() is False
