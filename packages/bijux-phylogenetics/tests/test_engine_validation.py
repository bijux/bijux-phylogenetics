from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines import run_model_selection
from bijux_phylogenetics.engines.validation import (
    audit_alignment_inference_readiness,
    classify_inference_workflow_failure,
    compare_inferred_tree_to_taxon_metadata,
    validate_ml_tree_contains_expected_taxa,
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


def _fake_iqtree_tree(path: Path) -> Path:
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
prefix.with_suffix(".treefile").write_text("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n", encoding="utf-8")
prefix.with_suffix(".iqtree").write_text("Tree inference completed\\n", encoding="utf-8")
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


def test_validate_ml_tree_contains_expected_taxa_matches_alignment_ids(tmp_path: Path) -> None:
    executable = _fake_iqtree_tree(tmp_path / "iqtree-tree-fixture")
    from bijux_phylogenetics.engines import run_maximum_likelihood_tree_inference

    workflow = run_maximum_likelihood_tree_inference(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=executable,
        prefix="example",
    )
    report = validate_ml_tree_contains_expected_taxa(workflow.manifest_path)
    assert report.valid is True
    assert report.expected_taxa == ["A", "B", "C", "D"]
    assert report.observed_taxa == ["A", "B", "C", "D"]


def test_compare_inferred_tree_to_taxon_metadata_reports_monophyly_and_splits(tmp_path: Path) -> None:
    metadata_path = tmp_path / "metadata.tsv"
    metadata_path.write_text(
        "taxon\tgroup\nA\tleft\nB\tleft\nC\tright\nD\tright\n",
        encoding="utf-8",
    )
    report = compare_inferred_tree_to_taxon_metadata(
        fixture("example_tree.nwk"),
        metadata_path,
        group_column="group",
    )
    assert report.monophyletic_group_count == 2

    metadata_path.write_text(
        "taxon\tgroup\nA\tmixed\nC\tmixed\nB\tother\nD\tother\n",
        encoding="utf-8",
    )
    split_report = compare_inferred_tree_to_taxon_metadata(
        fixture("example_tree.nwk"),
        metadata_path,
        group_column="group",
    )
    assert split_report.split_group_count >= 1
    assert any(row.status == "split_unexpectedly" for row in split_report.observations)


def test_classify_inference_workflow_failure_distinguishes_engine_and_parse_failures(tmp_path: Path) -> None:
    engine_failure = classify_inference_workflow_failure(
        workflow="maximum-likelihood-tree",
        input_paths=[fixture("example_alignment.fasta")],
        output_paths={},
        run_exit_code=3,
    )
    assert engine_failure.failure_category == "engine_failure"

    invalid_tree = tmp_path / "broken.treefile"
    invalid_tree.write_text("(A,B\n", encoding="utf-8")
    parse_failure = classify_inference_workflow_failure(
        workflow="maximum-likelihood-tree",
        input_paths=[fixture("example_alignment.fasta")],
        output_paths={"tree": invalid_tree},
        run_exit_code=0,
    )
    assert parse_failure.failure_category == "parse_failure"


def test_classify_inference_workflow_failure_detects_input_failures(tmp_path: Path) -> None:
    missing_input = tmp_path / "missing.fasta"
    report = classify_inference_workflow_failure(
        workflow="maximum-likelihood-tree",
        input_paths=[missing_input],
        output_paths={},
    )
    assert report.failure_category == "input_failure"
