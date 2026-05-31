from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_likelihood_multi_start_candidates,
    search_nucleotide_likelihood_multi_start_from_alignment,
    validate_likelihood_multi_start_evaluation_budget,
    write_nucleotide_likelihood_multi_start_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_multi_start_search_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_multi_start_from_alignment
        is search_nucleotide_likelihood_multi_start_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_multi_start_artifacts
        is write_nucleotide_likelihood_multi_start_artifacts
    )
    assert (
        likelihood_api.validate_likelihood_multi_start_evaluation_budget
        is validate_likelihood_multi_start_evaluation_budget
    )


def test_likelihood_multi_start_candidates_include_input_tree_and_distinct_random_starts() -> None:
    start_tree = load_tree(fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"))

    candidates = build_likelihood_multi_start_candidates(
        start_tree,
        start_tree_count=4,
        start_tree_source_policy="input-tree-plus-random-tree",
        start_tree_seed=17,
    )

    assert [candidate.source_kind for candidate in candidates] == [
        "input-tree",
        "random-tree",
        "random-tree",
        "random-tree",
    ]
    assert [candidate.source_label for candidate in candidates] == [
        "input-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
        "random-tree-seed-19",
    ]
    assert [candidate.generation_seed for candidate in candidates] == [None, 17, 18, 19]
    assert len({candidate.tree.to_newick() for candidate in candidates}) == 4


def test_likelihood_multi_start_search_reports_each_source_likelihood_hash_and_best_run() -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-multi-start-search"
    assert report.model_name == "JC69"
    assert report.local_search_method == "nni"
    assert report.input_tree_included is True
    assert report.generated_start_tree_count == 2
    assert report.start_tree_count == 3
    assert report.start_tree_seed == 17
    assert report.best_run_source_label == "input-tree"
    assert math.isclose(
        report.best_final_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert len(report.best_final_topology_fingerprint) == 64
    assert [row.start_tree_source_label for row in report.run_summaries] == [
        "input-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row.search_algorithm for row in report.run_summaries] == [
        "nucleotide-likelihood-nni-search",
        "nucleotide-likelihood-nni-search",
        "nucleotide-likelihood-nni-search",
    ]
    assert len({row.start_tree_newick for row in report.run_summaries}) == 3
    assert {
        row.final_topology_fingerprint for row in report.run_summaries
    } == {report.best_final_topology_fingerprint}
    assert [row.best_run for row in report.run_summaries] == [True, False, False]


def test_likelihood_multi_start_spr_forwards_budget_to_each_restart() -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        local_search_method="spr",
        start_tree_count=2,
        start_tree_seed=23,
        evaluation_budget=2,
        upper_branch_length_bound=1.0,
    )

    assert report.local_search_method == "spr"
    assert report.evaluation_budget == 2
    assert [row.search_algorithm for row in report.run_summaries] == [
        "nucleotide-likelihood-spr-search",
        "nucleotide-likelihood-spr-search",
    ]
    assert len({row.start_tree_source_label for row in report.run_summaries}) == 2
    assert sum(1 for row in report.run_summaries if row.best_run) == 1


def test_write_nucleotide_likelihood_multi_start_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_multi_start_artifacts(
        tmp_path / "likelihood-multi-start-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "start_tree_path",
        "best_tree_path",
        "summary_path",
        "run_json_path",
    }
    assert outputs["summary_path"].read_text(encoding="utf-8").startswith(
        "start_tree_source_kind\tstart_tree_source_label\tstart_tree_generation_seed\tsearch_algorithm\tstart_log_likelihood\tfinal_log_likelihood\tfinal_likelihood_rank\tfinal_topology_fingerprint\tsearch_iteration_count\taccepted_move_count\tevaluated_neighbor_count\tbranch_reoptimization_policy\tsubstitution_parameter_policy\tstopping_reason\tbest_run\tstart_tree_newick\tfinal_tree_newick\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-multi-start-search"
    assert payload["best_run_source_label"] == "input-tree"
    assert payload["generated_start_tree_count"] == 2
    assert len(payload["run_summaries"]) == 3
    assert [row["best_run"] for row in payload["run_summaries"]] == [True, False, False]


def test_likelihood_multi_start_rejects_nni_neighbor_budget() -> None:
    try:
        validate_likelihood_multi_start_evaluation_budget("nni", 1)
    except ValueError as error:
        assert str(error) == (
            "evaluation_budget is supported only when local_search_method is 'spr'"
        )
    else:
        raise AssertionError("NNI multi-start searches must reject SPR-only budgets")
