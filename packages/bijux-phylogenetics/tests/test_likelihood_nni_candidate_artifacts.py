from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
    write_nucleotide_likelihood_nni_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_nucleotide_likelihood_nni_artifacts_materializes_candidate_table(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_nni_artifacts(
        tmp_path / "likelihood-nni-best-improvement-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "start_tree_path",
        "final_tree_path",
        "best_tree_path",
        "trace_path",
        "candidate_table_path",
        "run_json_path",
    }
    assert outputs["best_tree_path"].name == "best_trees.nwk"
    assert len(load_newick_tree_set(outputs["best_tree_path"])) == 1
    assert outputs["best_tree_path"].read_text(encoding="utf-8").strip()
    assert outputs["trace_path"].read_text(encoding="utf-8").startswith(
        "event_index\tevent_kind\titeration\tmove_type\tcandidate_topology_fingerprint\tlog_likelihood_before\tlog_likelihood_after\tlog_likelihood_delta\taccepted_move\ttrace_reason\ttree_before_newick\ttree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\tstopping_reason\n"
    )
    assert outputs["candidate_table_path"].read_text(encoding="utf-8").startswith(
        "iteration\tcandidate_order\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tcandidate_tree_newick\tlog_likelihood\tlog_likelihood_delta\timproving_move\tselected_best_move\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    selected_rows = [
        row for row in payload["candidate_rows"] if row["selected_best_move"]
    ]
    assert payload["trace_rows"][0]["move_type"] == "nni"
    assert payload["trace_rows"][0]["accepted_move"] is False
    assert payload["trace_rows"][0]["trace_reason"] == "search-start"
    assert payload["trace_rows"][1]["accepted_move"] is True
    assert payload["trace_rows"][1]["trace_reason"] == "accepted-improving-move"
    assert payload["trace_rows"][-1]["trace_reason"] == "no-improving-neighbor"
    assert payload["accepted_move_count"] == 2
    assert len(payload["candidate_rows"]) == 12
    assert payload["candidate_rows"][0]["iteration"] == 1
    assert payload["candidate_rows"][0]["candidate_order"] == 1
    assert len(selected_rows) == 2
    assert [(row["iteration"], row["candidate_order"]) for row in selected_rows] == [
        (1, 4),
        (2, 2),
    ]
    assert all(
        row["selected_best_move"] is False
        for row in payload["candidate_rows"][8:]
    )
