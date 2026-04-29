from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines import run_model_selection
from bijux_phylogenetics.engines.validation import (
    audit_alignment_inference_readiness,
    validate_model_selection_against_engine_outputs,
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
prefix.with_suffix(".iqtree").write_text("Best-fit model according to BIC: GTR+G\\n", encoding="utf-8")
prefix.with_suffix(".model").write_text("Best-fit model: GTR+G\\n", encoding="utf-8")
raise SystemExit(0)
""",
    )


def test_audit_alignment_inference_readiness_prefers_ml_for_aligned_variable_data() -> None:
    report = audit_alignment_inference_readiness(fixture("example_alignment.fasta"))
    assert report.overall_decision == "ready"
    assert report.recommended_workflow == "maximum_likelihood"
    assert any(decision.workflow == "bayesian" and decision.ready for decision in report.decisions)


def test_audit_alignment_inference_readiness_blocks_unaligned_raw_sequences() -> None:
    report = audit_alignment_inference_readiness(fixture("example_sequences_raw.fasta"))
    assert report.overall_decision == "blocked"
    assert report.recommended_workflow == "unsuitable"
    assert any("not yet aligned" in blocker for decision in report.decisions for blocker in decision.blockers)


def test_validate_model_selection_against_engine_outputs_requires_exact_match(tmp_path: Path) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )
    report = validate_model_selection_against_engine_outputs(workflow.manifest_path)
    assert report.valid is True
    assert report.manifest_selected_model == "GTR+G"
    assert report.report_selected_model == "GTR+G"
    assert report.artifact_selected_model == "GTR+G"
