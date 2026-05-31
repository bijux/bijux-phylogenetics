from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
    search_nucleotide_likelihood_spr_from_alignment,
    validate_likelihood_spr_evaluation_budget,
    validate_nucleotide_likelihood_spr_search_budget,
    write_nucleotide_likelihood_spr_artifacts,
)
from bijux_phylogenetics.phylo.likelihood.validation import validate_explicit_branch_lengths
from bijux_phylogenetics.phylo.topology import (
    apply_rooted_spr_move,
    descendant_taxa,
    iter_rooted_spr_move_candidates,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_spr_search_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_spr_from_alignment
        is search_nucleotide_likelihood_spr_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_spr_artifacts
        is write_nucleotide_likelihood_spr_artifacts
    )
    assert (
        likelihood_api.validate_likelihood_spr_evaluation_budget
        is validate_likelihood_spr_evaluation_budget
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_spr_search_budget
        is validate_nucleotide_likelihood_spr_search_budget
    )


def test_rooted_spr_candidates_preserve_taxa_and_explicit_branch_lengths() -> None:
    start_tree = load_tree(fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"))
    candidates = list(iter_rooted_spr_move_candidates(start_tree))

    assert len(candidates) == 52
    assert len(
        {
            (candidate.pruned_clade_id, candidate.regraft_target_branch_id)
            for candidate in candidates
        }
    ) == len(candidates)

    start_taxa = descendant_taxa(start_tree.root)
    for candidate in candidates:
        moved_tree = apply_rooted_spr_move(start_tree, candidate)
        assert descendant_taxa(moved_tree.root) == start_taxa
        validate_explicit_branch_lengths(moved_tree, model_name="JC69")


def test_likelihood_spr_search_improves_when_likelihood_nni_stalls() -> None:
    tree_path = fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk")
    alignment_path = fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta")

    nni_report = search_nucleotide_likelihood_nni_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )
    spr_report = search_nucleotide_likelihood_spr_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
        upper_branch_length_bound=1.0,
        evaluation_budget=2,
    )

    assert nni_report.accepted_move_count == 0
    assert nni_report.evaluated_neighbor_count == 6
    assert nni_report.stopping_reason == "no-improving-neighbor"

    assert spr_report.algorithm == "nucleotide-likelihood-spr-search"
    assert spr_report.model_name == "JC69"
    assert spr_report.accepted_move_count == 1
    assert spr_report.iteration_count == 1
    assert spr_report.evaluated_neighbor_count == 2
    assert spr_report.evaluation_budget == 2
    assert spr_report.search_budget.max_candidate_count == 2
    assert spr_report.search_budget.max_iteration_count is None
    assert spr_report.search_budget.max_elapsed_seconds is None
    assert spr_report.search_budget.max_accepted_move_count is None
    assert spr_report.branch_reoptimization_policy == "coordinate-branch-lengths"
    assert spr_report.substitution_parameter_policy == "fixed-from-model"
    assert spr_report.substitution_parameter_values == {}
    assert spr_report.substitution_parameter_warnings == []
    assert spr_report.stopping_reason == "candidate-budget-exhausted"
    assert spr_report.unsearched_candidate_count == 43
    assert spr_report.final_tree_newick != nni_report.final_tree_newick
    assert math.isclose(
        spr_report.start_log_likelihood,
        -8.174114933271198,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        spr_report.final_log_likelihood,
        -8.17411493278412,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert spr_report.final_log_likelihood > nni_report.final_log_likelihood
    assert [row.event_kind for row in spr_report.trace_rows] == [
        "start",
        "accepted-move",
        "final",
    ]
    assert spr_report.trace_rows[0].optimized_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "E",
        "A|D",
        "A|B|D",
        "A|B|C|D",
    ]
    assert spr_report.trace_rows[1].pruned_clade_id == "A"
    assert spr_report.trace_rows[1].regraft_target_branch_id == "B"
    assert spr_report.trace_rows[1].affected_branch_clade_ids == ["A", "A|B", "A|D", "B"]
    assert spr_report.trace_rows[1].optimized_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "E",
        "A|B",
        "A|B|D",
        "A|B|C|D",
    ]
    assert spr_report.trace_rows[1].branch_reoptimization_converged is True
    assert spr_report.trace_rows[2].stopping_reason == "candidate-budget-exhausted"
    assert spr_report.trace_rows[2].unsearched_candidate_count == 43
    assert any(
        "search boundary" in warning
        for warning in spr_report.trace_rows[0].boundary_warning_messages
    )
    assert spr_report.trace_rows[-1].boundary_warning_messages == []


def test_likelihood_spr_search_local_branch_reoptimization_reports_affected_scope() -> None:
    tree_path = fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk")
    alignment_path = fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta")

    full_report = search_nucleotide_likelihood_spr_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
        upper_branch_length_bound=1.0,
        evaluation_budget=2,
    )
    local_report = search_nucleotide_likelihood_spr_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
        upper_branch_length_bound=1.0,
        evaluation_budget=2,
        branch_reoptimization_policy="spr-local-affected-branches",
    )

    assert local_report.branch_reoptimization_policy == "spr-local-affected-branches"
    assert local_report.trace_rows[0].branch_reoptimization_scope == "all-branches"
    assert local_report.trace_rows[0].branch_reoptimization_converged is True
    accepted_rows = [
        row for row in local_report.trace_rows if row.event_kind == "accepted-move"
    ]
    assert accepted_rows
    assert accepted_rows[0].branch_reoptimization_scope == "local-affected-branches"
    assert accepted_rows[0].affected_branch_clade_ids == ["A", "A|B", "A|D", "B"]
    assert accepted_rows[0].optimized_branch_clade_ids == ["A", "B", "A|B"]
    assert accepted_rows[0].optimized_branch_count == 3
    assert accepted_rows[0].branch_reoptimization_converged is True
    assert (
        accepted_rows[0].optimized_branch_count
        < full_report.trace_rows[1].optimized_branch_count
    )
    assert local_report.total_branch_function_evaluation_count < full_report.total_branch_function_evaluation_count
    assert math.isclose(
        local_report.final_log_likelihood,
        full_report.final_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_write_nucleotide_likelihood_spr_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_spr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
        evaluation_budget=2,
    )

    outputs = write_nucleotide_likelihood_spr_artifacts(
        tmp_path / "likelihood-spr-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "start_tree_path",
        "final_tree_path",
        "best_tree_path",
        "trace_path",
        "run_json_path",
    }
    assert outputs["best_tree_path"].name == "best_trees.nwk"
    assert len(load_newick_tree_set(outputs["best_tree_path"])) == 1
    assert outputs["trace_path"].read_text(encoding="utf-8").startswith(
        "event_index\tevent_kind\titeration\tmove_type\tcandidate_topology_fingerprint\tlog_likelihood_before\tlog_likelihood_after\tlog_likelihood_delta\taccepted_move\ttrace_reason\ttree_before_newick\ttree_after_newick\tpruned_clade_id\tregraft_target_branch_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\taffected_branch_clade_ids\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\tstopping_reason\tunsearched_candidate_count\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-spr-search"
    assert payload["evaluation_budget"] == 2
    assert payload["search_budget"] == {
        "max_accepted_move_count": None,
        "max_candidate_count": 2,
        "max_elapsed_seconds": None,
        "max_iteration_count": None,
    }
    assert payload["accepted_move_count"] == 1
    assert payload["iteration_count"] == 1
    assert payload["evaluated_neighbor_count"] == 2
    assert payload["stopping_reason"] == "candidate-budget-exhausted"
    assert payload["unsearched_candidate_count"] == 43
    assert payload["trace_rows"][0]["move_type"] == "spr"
    assert payload["trace_rows"][0]["accepted_move"] is False
    assert payload["trace_rows"][1]["accepted_move"] is True
    assert payload["trace_rows"][1]["trace_reason"] == "accepted-improving-move"
    assert payload["trace_rows"][-1]["trace_reason"] == "candidate-budget-exhausted"
    assert payload["trace_rows"][1]["branch_reoptimization_scope"] == "all-branches"
    assert payload["trace_rows"][1]["affected_branch_clade_ids"] == ["A", "A|B", "A|D", "B"]
    assert payload["trace_rows"][1]["optimized_branch_count"] == 8
    assert payload["trace_rows"][1]["branch_reoptimization_converged"] is True
    assert payload["trace_rows"][1]["boundary_warning_messages"]
    assert payload["trace_rows"][2]["unsearched_candidate_count"] == 43
