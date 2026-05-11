from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines import run_model_selection
from bijux_phylogenetics.engines.validation import (
    audit_alignment_inference_readiness,
    classify_inference_workflow_failure,
    compare_inferred_tree_to_taxon_metadata,
    compare_inferred_trees_across_engines,
    compare_ml_trees_across_models,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    validate_bootstrap_tree_set,
    validate_inference_engine_outputs,
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
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model according to BIC: GTR+G\\nLog-likelihood of the tree: -123.456\\n",
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
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: GTR+G\\nLog-likelihood of the tree: -345.678\\nTree inference completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -345.678\\n",
    encoding="utf-8",
)
raise SystemExit(0)
""",
    )


def _fake_iqtree_tree_alt(path: Path, *, model: str, tree_newick: str) -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("IQ-TREE multicore version 2.9.9")
    raise SystemExit(0)

prefix = Path(args[args.index("-pre") + 1])
prefix.parent.mkdir(parents=True, exist_ok=True)
prefix.with_suffix(".treefile").write_text({tree_newick!r} + "\\n", encoding="utf-8")
prefix.with_suffix(".iqtree").write_text(
    "Best-fit model: {model}\\nLog-likelihood of the tree: -456.789\\nTree inference completed\\n",
    encoding="utf-8",
)
prefix.with_suffix(".log").write_text(
    "IQ-TREE fixture inference log\\nBEST SCORE FOUND : -456.789\\n",
    encoding="utf-8",
)
prefix.with_suffix(".model").write_text("Best-fit model: {model}\\n", encoding="utf-8")
raise SystemExit(0)
""",
    )


def _fake_fasttree_tree(path: Path, *, tree_newick: str) -> Path:
    return _write_executable(
        path,
        f"""#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 fixture")
    raise SystemExit(0)

print({tree_newick!r})
""",
    )


def test_audit_alignment_inference_readiness_prefers_ml_for_aligned_variable_data() -> (
    None
):
    report = audit_alignment_inference_readiness(fixture("example_alignment.fasta"))
    assert report.overall_decision == "ready"
    assert report.recommended_workflow == "maximum_likelihood"
    assert any(
        decision.workflow == "bayesian" and decision.ready
        for decision in report.decisions
    )


def test_audit_alignment_inference_readiness_blocks_unaligned_raw_sequences() -> None:
    report = audit_alignment_inference_readiness(fixture("example_sequences_raw.fasta"))
    assert report.overall_decision == "blocked"
    assert report.recommended_workflow == "unsuitable"
    assert any(
        "not yet aligned" in blocker
        for decision in report.decisions
        for blocker in decision.blockers
    )


def test_validate_model_selection_against_engine_outputs_requires_exact_match(
    tmp_path: Path,
) -> None:
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


def test_validate_ml_tree_contains_expected_taxa_matches_alignment_ids(
    tmp_path: Path,
) -> None:
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


def test_compare_inferred_tree_to_taxon_metadata_reports_monophyly_and_splits(
    tmp_path: Path,
) -> None:
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


def test_classify_inference_workflow_failure_distinguishes_engine_and_parse_failures(
    tmp_path: Path,
) -> None:
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


def test_classify_inference_workflow_failure_detects_input_failures(
    tmp_path: Path,
) -> None:
    missing_input = tmp_path / "missing.fasta"
    report = classify_inference_workflow_failure(
        workflow="maximum-likelihood-tree",
        input_paths=[missing_input],
        output_paths={},
    )
    assert report.failure_category == "input_failure"


def test_validate_bootstrap_tree_set_requires_consistent_taxa(tmp_path: Path) -> None:
    valid_path = tmp_path / "valid.ufboot"
    valid_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n",
        encoding="utf-8",
    )
    valid_report = validate_bootstrap_tree_set(valid_path)
    assert valid_report.valid is True
    assert valid_report.tree_count == 2

    invalid_path = tmp_path / "invalid.ufboot"
    invalid_path.write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\n",
        encoding="utf-8",
    )
    invalid_report = validate_bootstrap_tree_set(invalid_path)
    assert invalid_report.valid is False
    assert any("same taxon set" in issue for issue in invalid_report.issues)


def test_summarize_bootstrap_support_distribution_reports_range_median_and_histogram(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "supported.nwk"
    tree_path.write_text(
        "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)68:0.2)72:0.3;\n", encoding="utf-8"
    )

    report = summarize_bootstrap_support_distribution(tree_path)

    assert report.supported_node_count == 3
    assert report.minimum_support == 68.0
    assert report.maximum_support == 95.0
    assert report.median_support == 72.0
    assert report.weakly_supported_clade_count == 1
    assert report.support_histogram == {"lt50": 0, "50to69": 1, "70to89": 1, "ge90": 1}


def test_detect_weakly_supported_backbone_flags_major_internal_branches(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "weak-backbone.nwk"
    tree_path.write_text(
        "(((A:0.1,B:0.1)95:0.2,C:0.1)62:0.3,(D:0.1,E:0.1)91:0.2)58:0.4;\n",
        encoding="utf-8",
    )

    report = detect_weakly_supported_backbone(tree_path, threshold=70.0)

    assert report.evaluated_backbone_node_count >= 2
    assert report.weak_backbone_node_count == 2
    assert any(node.node == "A|B|C|D|E" for node in report.weak_nodes)
    assert any("backbone" in warning for warning in report.warnings)


def test_compare_ml_trees_across_models_reports_topology_and_branch_length_differences(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import run_maximum_likelihood_tree_inference

    left_executable = _fake_iqtree_tree_alt(
        tmp_path / "iqtree-left",
        model="GTR+G",
        tree_newick="((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);",
    )
    right_executable = _fake_iqtree_tree_alt(
        tmp_path / "iqtree-right",
        model="HKY+G",
        tree_newick="((A:0.1,C:0.1):0.25,(B:0.1,D:0.1):0.25);",
    )
    left = run_maximum_likelihood_tree_inference(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "left-ml",
        model="GTR+G",
        executable=left_executable,
        prefix="left",
    )
    right = run_maximum_likelihood_tree_inference(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "right-ml",
        model="HKY+G",
        executable=right_executable,
        prefix="right",
    )

    report = compare_ml_trees_across_models(left.manifest_path, right.manifest_path)

    assert report.comparison_kind == "model"
    assert report.left_selected_model == "GTR+G"
    assert report.right_selected_model == "HKY+G"
    assert report.topology.topology_equal is False
    assert any("topologies differ" in warning for warning in report.warnings)


def test_compare_inferred_trees_across_engines_reports_engine_labels(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import (
        run_fast_tree_inference,
        run_maximum_likelihood_tree_inference,
    )

    iqtree_executable = _fake_iqtree_tree_alt(
        tmp_path / "iqtree-engine",
        model="GTR+G",
        tree_newick="((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);",
    )
    fasttree_executable = _fake_fasttree_tree(
        tmp_path / "fasttree-engine",
        tree_newick="((A:0.1,B:0.1):0.3,(C:0.1,D:0.1):0.3);",
    )
    ml = run_maximum_likelihood_tree_inference(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=iqtree_executable,
        prefix="engine",
    )
    fast = run_fast_tree_inference(
        fixture("example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=fasttree_executable,
    )

    report = compare_inferred_trees_across_engines(fast.manifest_path, ml.manifest_path)

    assert report.comparison_kind == "engine"
    assert report.left_label == "FastTree"
    assert report.right_label == "IQ-TREE"
    assert report.topology.topology_equal is True
    assert report.branch_lengths.shared_splits
    assert any("branch-length" in warning for warning in report.warnings)


def test_validate_inference_engine_outputs_checks_manifest_consistency(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree_tree(tmp_path / "iqtree-tree-fixture")
    from bijux_phylogenetics.engines import run_maximum_likelihood_tree_inference

    workflow = run_maximum_likelihood_tree_inference(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "ml",
        model="GTR+G",
        executable=executable,
        prefix="example",
    )
    report = validate_inference_engine_outputs(workflow.manifest_path)
    assert report.valid is True
    assert report.current_output_checksum_match is True

    workflow.output_paths["tree"].write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,X:0.1):0.2);\n", encoding="utf-8"
    )
    drift_report = validate_inference_engine_outputs(workflow.manifest_path)
    assert drift_report.valid is False
    assert any(
        "checksums" in issue or "expected taxa" in issue
        for issue in drift_report.issues
    )
