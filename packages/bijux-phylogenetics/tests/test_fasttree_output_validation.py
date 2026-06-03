from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.engines import run_fast_tree_inference

pytestmark = pytest.mark.engine_contract

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fake_fasttree_partial_support_labels(path: Path) -> Path:
    return _write_executable(
        path,
        """#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

print("((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1):0.3);")
""",
    )


def test_run_fast_tree_inference_preserves_partial_support_annotations_in_report(
    tmp_path: Path,
) -> None:
    executable = _fake_fasttree_partial_support_labels(
        tmp_path / "fasttree-partial-support"
    )

    report = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=executable,
    )

    assert report.fasttree_support_summary is not None
    assert report.fasttree_support_summary.annotated_node_count == 1
    assert any(
        "did not expose parsable FastTree local support labels" in warning
        for warning in report.fasttree_support_summary.warnings
    )
