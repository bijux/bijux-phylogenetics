from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
    write_nucleotide_likelihood_nni_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_search_uses_likelihood_objective_and_reaches_local_optimum(
    monkeypatch,
) -> None:
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError(
            "parsimony NNI must not be reused by likelihood NNI search"
        )

    monkeypatch.setattr(parsimony_api, "search_parsimony_nni", fail_if_called)

    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-nni-search"
    assert report.model_name == "JC69"
    assert report.taxon_count == 4
    assert report.site_count == 12
    assert report.pattern_count == 2
    assert report.input_tree_newick == "(((A:0.1,C:0.1):0.1,B:0.1):0.1,D:0.1);"
    assert report.branch_reoptimization_policy == "coordinate-branch-lengths"
    assert report.substitution_parameter_policy == "fixed-from-model"
    assert report.substitution_parameter_values == {}
    assert report.substitution_parameter_warnings == []
    assert report.start_tree_newick == (
        "(((A:2.43539016857288e-10,C:0.999999999756461):2.43539016857288e-10,"
        "B:2.43539016857288e-10):0.999999999756461,D:0.999999999756461);"
    )
    assert math.isclose(
        report.start_log_likelihood,
        -54.442517349645854,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.final_tree_newick == (
        "((A:2.43539016857288e-10,B:2.43539016857288e-10):0.999999999756461,"
        "(C:2.43539016857288e-10,D:2.43539016857288e-10):0.999999999756461);"
    )
    assert math.isclose(
        report.final_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.accepted_move_count == 2
    assert report.evaluated_neighbor_count == 12
    assert report.total_branch_optimization_pass_count == 27
    assert report.total_branch_function_evaluation_count == 7627
    assert report.stopping_reason == "no-improving-neighbor"
    assert [row.event_kind for row in report.trace_rows] == [
        "start",
        "accepted-move",
        "accepted-move",
        "final",
    ]
    assert report.trace_rows[0].optimized_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "A|C",
        "A|B|C",
    ]
    assert report.trace_rows[1].pivot_branch_id == "A|C"
    assert report.trace_rows[1].sibling_clade_id == "B"
    assert report.trace_rows[1].exchanged_clade_id == "C"
    assert report.trace_rows[1].optimized_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "A|B",
        "A|B|C",
    ]
    assert report.trace_rows[2].pivot_branch_id == "A|B|C"
    assert report.trace_rows[2].sibling_clade_id == "D"
    assert report.trace_rows[2].exchanged_clade_id == "A|B"
    assert report.trace_rows[2].optimized_branch_clade_ids == [
        "A",
        "B",
        "C",
        "D",
        "A|B",
        "C|D",
    ]
    assert report.trace_rows[3].stopping_reason == "no-improving-neighbor"
    assert any(
        "search boundary" in warning
        for warning in report.trace_rows[0].boundary_warning_messages
    )
    assert any(
        "search boundary" in warning
        for warning in report.trace_rows[1].boundary_warning_messages
    )
    assert report.trace_rows[-1].boundary_warning_messages == []
    assert math.isclose(
        report.trace_rows[1].log_likelihood_delta or 0.0,
        17.499717347003227,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.trace_rows[2].log_likelihood_delta or 0.0,
        2.8075503046659165,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.trace_rows[-1].tree_after_newick == report.final_tree_newick


def test_likelihood_nni_search_local_branch_reoptimization_reports_neighborhood_scope() -> (
    None
):
    full_report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )
    local_report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        branch_reoptimization_policy="nni-local-affected-branches",
        upper_branch_length_bound=1.0,
    )

    assert local_report.branch_reoptimization_policy == "nni-local-affected-branches"
    assert local_report.trace_rows[0].branch_reoptimization_scope == "all-branches"
    assert local_report.trace_rows[0].optimized_branch_count == 6
    assert local_report.trace_rows[0].branch_reoptimization_converged is True
    accepted_rows = [
        row for row in local_report.trace_rows if row.event_kind == "accepted-move"
    ]
    assert accepted_rows
    assert all(
        row.branch_reoptimization_scope == "local-nni-neighborhood"
        for row in accepted_rows
    )
    assert all(
        row.optimized_branch_count < local_report.trace_rows[0].optimized_branch_count
        for row in accepted_rows
    )
    assert set(accepted_rows[0].optimized_branch_clade_ids) == {
        "A",
        "B",
        "C",
        "A|B",
        "A|B|C",
    }
    assert set(accepted_rows[-1].optimized_branch_clade_ids) == {
        "C",
        "D",
        "A|B",
        "C|D",
    }
    assert local_report.trace_rows[-1].branch_reoptimization_scope == "none"
    assert local_report.trace_rows[-1].branch_reoptimization_converged is None
    assert (
        local_report.trace_rows[1].optimized_branch_count
        < full_report.trace_rows[1].optimized_branch_count
    )
    assert (
        local_report.total_branch_function_evaluation_count
        < full_report.total_branch_function_evaluation_count
    )


def test_write_nucleotide_likelihood_nni_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_nni_artifacts(
        tmp_path / "likelihood-nni-run",
        report,
    )

    assert set(outputs) == {
        "best_tree_path",
        "candidate_table_path",
        "input_tree_path",
        "start_tree_path",
        "final_tree_path",
        "trace_path",
        "run_json_path",
    }
    assert outputs["best_tree_path"].name == "best_trees.nwk"
    assert outputs["candidate_table_path"].name == "candidate_table.tsv"
    assert (
        outputs["trace_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "event_index\tevent_kind\titeration\tmove_type\tcandidate_topology_fingerprint\tlog_likelihood_before\tlog_likelihood_after\tlog_likelihood_delta\taccepted_move\ttrace_reason\ttree_before_newick\ttree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\tstopping_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-nni-search"
    assert payload["model_name"] == "JC69"
    assert payload["branch_reoptimization_policy"] == "coordinate-branch-lengths"
    assert payload["substitution_parameter_policy"] == "fixed-from-model"
    assert payload["accepted_move_count"] == 2
    assert payload["evaluated_neighbor_count"] == 12
    assert payload["total_branch_optimization_pass_count"] == 27
    assert payload["total_branch_function_evaluation_count"] == 7627
    assert payload["stopping_reason"] == "no-improving-neighbor"
    assert payload["trace_rows"][0]["branch_reoptimization_scope"] == "all-branches"
    assert payload["trace_rows"][0]["optimized_branch_count"] == 6
    assert payload["trace_rows"][0]["optimized_branch_clade_ids"] == [
        "A",
        "B",
        "C",
        "D",
        "A|C",
        "A|B|C",
    ]
    assert payload["trace_rows"][0]["branch_reoptimization_converged"] is True
    assert payload["trace_rows"][0]["boundary_warning_messages"]
