from __future__ import annotations

import json
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
    summarize_sh_alrt_support_distribution,
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
    " No. Model         -LnL         df  AIC          AICc         BIC\\n"
    "  1  GTR+G         123.456      12  270.912      330.912      272.912\\n"
    "  2  HKY+G         124.000      10  268.000      320.000      269.000\\n"
    "  3  JC            130.500      5   271.000      300.000      271.500\\n"
    "Akaike Information Criterion:           HKY+G\\n"
    "Corrected Akaike Information Criterion: JC\\n"
    "Bayesian Information Criterion:         GTR+G\\n"
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
if "-alrt" in args:
    prefix.with_suffix(".treefile").write_text(
        "((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1)79/96:0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -222.222\\nSH-aLRT and ultrafast bootstrap analysis completed\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture sh-alrt log\\nBEST SCORE FOUND : -222.222\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)
if "-bb" in args:
    support_tree = "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\\n"
    prefix.with_suffix(".treefile").write_text(support_tree, encoding="utf-8")
    prefix.with_suffix(".contree").write_text(support_tree, encoding="utf-8")
    prefix.with_suffix(".ufboot").write_text(
        "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".iqtree").write_text(
        "Best-fit model: GTR+G\\nLog-likelihood of the tree: -234.567\\nBootstrap analysis completed\\n",
        encoding="utf-8",
    )
    prefix.with_suffix(".log").write_text(
        "IQ-TREE fixture bootstrap log\\nBEST SCORE FOUND : -234.567\\n",
        encoding="utf-8",
    )
    raise SystemExit(0)
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
    assert report.manifest_selected_criterion == "BIC"
    assert report.report_selected_model == "GTR+G"
    assert report.report_selected_criterion == "BIC"
    assert report.artifact_selected_model == "GTR+G"
    assert report.candidate_model_count == 3
    assert report.best_model_aic == "HKY+G"
    assert report.best_model_aicc == "JC"
    assert report.best_model_bic == "GTR+G"


def test_validate_model_selection_against_engine_outputs_requires_manifest_likelihood(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )
    payload = json.loads(workflow.manifest_path.read_text(encoding="utf-8"))
    payload["log_likelihood"] = None
    workflow.manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_model_selection_against_engine_outputs(workflow.manifest_path)

    assert report.valid is False
    assert "manifest log_likelihood field is missing" in report.issues


def test_validate_model_selection_against_engine_outputs_requires_candidate_summary(
    tmp_path: Path,
) -> None:
    executable = _fake_iqtree(tmp_path / "iqtree-fixture")
    workflow = run_model_selection(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executable,
        prefix="example",
    )
    payload = json.loads(workflow.manifest_path.read_text(encoding="utf-8"))
    payload["model_selection_summary"]["candidate_count"] = 0
    workflow.manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_model_selection_against_engine_outputs(workflow.manifest_path)

    assert report.valid is False
    assert "manifest does not record any candidate substitution models" in report.issues


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


def test_classify_inference_workflow_failure_explains_invalid_fasta_records() -> None:
    report = classify_inference_workflow_failure(
        workflow="multiple-sequence-alignment",
        input_paths=[fixture("example_sequences_invalid_input.fasta")],
        output_paths={},
    )

    assert report.failure_category == "input_failure"
    assert report.failure_reason == "invalid_fasta_input"
    assert "record-level biological data problems" in report.scientific_explanation
    assert report.evidence["duplicate_identifier_count"] == 1
    assert report.evidence["empty_sequence_count"] == 1
    assert report.evidence["illegal_character_count"] == 1


def test_classify_inference_workflow_failure_explains_missing_tree_outputs(
    tmp_path: Path,
) -> None:
    report = classify_inference_workflow_failure(
        workflow="maximum-likelihood-tree",
        input_paths=[fixture("example_alignment.fasta")],
        output_paths={"tree": tmp_path / "missing.treefile"},
        run_exit_code=0,
    )

    assert report.failure_category == "missing_output"
    assert report.failure_reason == "tree_output_missing"
    assert "without writing a tree artifact" in report.scientific_explanation
    assert report.evidence["missing_outputs"] == [
        {
            "output_name": "tree",
            "path": str(tmp_path / "missing.treefile"),
        }
    ]


def test_classify_inference_workflow_failure_explains_empty_trimmed_alignment(
    tmp_path: Path,
) -> None:
    trimmed_path = tmp_path / "trimmed.fasta"
    trimmed_path.write_text("", encoding="utf-8")

    report = classify_inference_workflow_failure(
        workflow="alignment-trimming",
        input_paths=[fixture("example_alignment_trim.fasta")],
        output_paths={"trimmed_alignment": trimmed_path},
        run_exit_code=0,
    )

    assert report.failure_category == "invalid_output"
    assert report.failure_reason == "trimmed_alignment_empty"
    assert "removed all usable alignment signal" in report.scientific_explanation
    assert report.evidence["invalid_outputs"] == [
        {
            "output_name": "trimmed_alignment",
            "path": str(trimmed_path),
        }
    ]


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

    assert report.internal_node_count == 2
    assert report.supported_node_count == 2
    assert report.minimum_support == 68.0
    assert report.maximum_support == 95.0
    assert report.median_support == 81.5
    assert report.weakly_supported_clade_count == 1
    assert report.support_histogram == {"lt50": 0, "50to69": 1, "70to89": 0, "ge90": 1}


def test_summarize_bootstrap_support_distribution_treats_collapsed_unlabeled_branches_as_zero_support(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "collapsed-unlabeled.nwk"
    tree_path.write_text(
        "(((A:0.0,B:0.0):0.000001,C:0.1)78:0.2,D:0.1)90:0.3;\n",
        encoding="utf-8",
    )

    report = summarize_bootstrap_support_distribution(tree_path)

    assert report.internal_node_count == 2
    assert report.supported_node_count == 2
    assert report.minimum_support == 0.0
    assert report.maximum_support == 78.0
    assert report.median_support == 39.0
    assert report.weakly_supported_clade_count == 1
    assert report.support_histogram == {"lt50": 1, "50to69": 0, "70to89": 1, "ge90": 0}
    assert any(
        node.descendant_taxa == ["A", "B"] and node.support == 0.0
        for node in report.nodes
    )
    assert (
        "one or more internal nodes did not expose numeric support labels"
        not in report.warnings
    )


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


def test_summarize_sh_alrt_support_distribution_tracks_conflicting_support_signals(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "sh-alrt.treefile"
    tree_path.write_text(
        "((A:0.1,B:0.1)82/97:0.2,(C:0.1,D:0.1)79/96:0.2);\n",
        encoding="utf-8",
    )

    report = summarize_sh_alrt_support_distribution(tree_path)

    assert report.annotated_node_count == 2
    assert report.fully_scored_node_count == 2
    assert report.minimum_sh_alrt_support == 79.0
    assert report.maximum_sh_alrt_support == 82.0
    assert report.minimum_ufboot_support == 96.0
    assert report.maximum_ufboot_support == 97.0
    assert report.weak_sh_alrt_clade_count == 1
    assert report.weak_ufboot_clade_count == 0
    assert report.conflicting_support_signal_count == 1
    assert [node.support_agreement for node in report.nodes] == [
        "both_strong",
        "ufboot_only",
    ]
    assert (
        "one or more internal clades show conflicting sh-alrt and ufboot support signals"
        in report.warnings
    )


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
        tree_newick="((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);",
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


def test_compare_inferred_trees_across_engines_accepts_support_tree_manifests(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import (
        run_bootstrap_support_estimation,
        run_fast_tree_inference,
    )

    iqtree_executable = _fake_iqtree_tree(tmp_path / "iqtree-engine")
    fasttree_executable = _fake_fasttree_tree(
        tmp_path / "fasttree-engine",
        tree_newick="((A:0.1,B:0.1)0.98:0.3,(C:0.1,D:0.1)0.62:0.3);",
    )
    iqtree = run_bootstrap_support_estimation(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=iqtree_executable,
        prefix="engine",
        replicates=1000,
    )
    fast = run_fast_tree_inference(
        fixture("example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=fasttree_executable,
    )

    report = compare_inferred_trees_across_engines(
        fast.manifest_path, iqtree.manifest_path
    )

    assert report.comparison_kind == "engine"
    assert report.left_label == "FastTree"
    assert report.right_label == "IQ-TREE"
    assert report.topology.topology_equal is True
    assert report.support.shared_clades


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


def test_validate_inference_engine_outputs_requires_fasttree_support_evidence(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import run_fast_tree_inference

    executable = _fake_fasttree_tree(
        tmp_path / "fasttree-engine",
        tree_newick="((A:0.1,B:0.1)0.96:0.3,(C:0.1,D:0.1)0.64:0.3);",
    )
    workflow = run_fast_tree_inference(
        fixture("example_alignment.fasta"),
        tmp_path / "fasttree.nwk",
        executable=executable,
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


def test_validate_inference_engine_outputs_requires_bootstrap_support_summary(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import run_bootstrap_support_estimation

    executable = _fake_iqtree_tree(tmp_path / "iqtree-tree-fixture")
    workflow = run_bootstrap_support_estimation(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=executable,
        prefix="example",
        replicates=1000,
    )
    payload = json.loads(workflow.manifest_path.read_text(encoding="utf-8"))
    payload["iqtree_summary"]["support_value_count"] = 0
    payload["iqtree_summary"]["support_values"] = []
    workflow.manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_inference_engine_outputs(workflow.manifest_path)

    assert report.valid is False
    assert (
        "bootstrap-support manifest does not record parsed support values"
        in report.issues
    )


def test_validate_inference_engine_outputs_requires_bootstrap_review_artifacts(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import run_bootstrap_support_estimation

    executable = _fake_iqtree_tree(tmp_path / "iqtree-tree-fixture")
    workflow = run_bootstrap_support_estimation(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model="GTR+G",
        executable=executable,
        prefix="example",
        replicates=1000,
    )
    payload = json.loads(workflow.manifest_path.read_text(encoding="utf-8"))
    del payload["output_paths"]["support_table"]
    del payload["output_paths"]["low_support_branches"]
    del payload["output_paths"]["support_histogram"]
    payload["bootstrap_support_summary"] = None
    payload["weak_backbone_report"] = None
    workflow.manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_inference_engine_outputs(workflow.manifest_path)

    assert report.valid is False
    assert (
        "bootstrap-support manifest is missing the support_table output"
        in report.issues
    )
    assert (
        "bootstrap-support manifest is missing the low_support_branches output"
        in report.issues
    )
    assert (
        "bootstrap-support manifest is missing the support_histogram output"
        in report.issues
    )
    assert (
        "bootstrap-support manifest is missing the bootstrap_support_summary"
        in report.issues
    )
    assert (
        "bootstrap-support manifest is missing the weak_backbone_report"
        in report.issues
    )


def test_validate_inference_engine_outputs_requires_sh_alrt_review_evidence(
    tmp_path: Path,
) -> None:
    from bijux_phylogenetics.engines import run_sh_alrt_support_estimation

    executable = _fake_iqtree_tree(tmp_path / "iqtree-tree-fixture")
    workflow = run_sh_alrt_support_estimation(
        fixture("example_alignment.fasta"),
        out_dir=tmp_path / "sh-alrt",
        model="GTR+G",
        executable=executable,
        prefix="example",
        sh_alrt_replicates=1000,
        bootstrap_replicates=1000,
    )
    payload = json.loads(workflow.manifest_path.read_text(encoding="utf-8"))
    del payload["output_paths"]["support_table"]
    del payload["output_paths"]["conflicting_support_branches"]
    payload["sh_alrt_support_summary"] = None
    workflow.manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_inference_engine_outputs(workflow.manifest_path)

    assert report.valid is False
    assert (
        "sh-alrt-support manifest is missing the support_table output" in report.issues
    )
    assert (
        "sh-alrt-support manifest is missing the conflicting_support_branches output"
        in report.issues
    )
    assert (
        "sh-alrt-support manifest is missing the sh_alrt_support_summary"
        in report.issues
    )
