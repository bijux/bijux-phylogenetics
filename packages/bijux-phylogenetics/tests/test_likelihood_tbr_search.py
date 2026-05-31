from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick_tree_set
import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_tbr_from_alignment,
    write_nucleotide_likelihood_tbr_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_tbr_search_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_tbr_from_alignment
        is search_nucleotide_likelihood_tbr_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tbr_artifacts
        is write_nucleotide_likelihood_tbr_artifacts
    )


def test_likelihood_tbr_search_improves_and_reports_best_move() -> None:
    report = search_nucleotide_likelihood_tbr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-tbr-search"
    assert report.model_name == "JC69"
    assert report.accepted_move_count == 1
    assert report.evaluated_neighbor_count == 20
    assert report.branch_reoptimization_policy == "coordinate-branch-lengths"
    assert report.substitution_parameter_policy == "fixed-from-model"
    assert report.substitution_parameter_values == {}
    assert report.substitution_parameter_warnings == []
    assert report.stopping_reason == "no-improving-neighbor"
    assert math.isclose(
        report.start_log_likelihood,
        -8.174114933271198,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.final_log_likelihood,
        -7.742402024467651,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.final_log_likelihood > report.start_log_likelihood
    assert [row.event_kind for row in report.trace_rows] == [
        "start",
        "accepted-move",
        "final",
    ]
    assert report.trace_rows[1].cut_edge_id == "A|B|C|D"
    assert report.trace_rows[1].left_attachment_branch_id == "D"
    assert report.trace_rows[1].right_attachment_branch_id == "interface"
    assert report.trace_rows[1].optimized_branch_count == 8
    assert report.trace_rows[1].branch_reoptimization_scope == "all-branches"
    assert report.trace_rows[1].branch_reoptimization_converged is True
    assert report.trace_rows[2].stopping_reason == "no-improving-neighbor"


def test_write_nucleotide_likelihood_tbr_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_tbr_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_tbr_artifacts(
        tmp_path / "likelihood-tbr-run",
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
    assert (
        outputs["trace_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "event_index\tevent_kind\titeration\tmove_type\tcandidate_topology_fingerprint\tlog_likelihood_before\tlog_likelihood_after\tlog_likelihood_delta\taccepted_move\ttrace_reason\ttree_before_newick\ttree_after_newick\tcut_edge_id\tleft_attachment_branch_id\tright_attachment_branch_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\tstopping_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-tbr-search"
    assert payload["accepted_move_count"] == 1
    assert payload["evaluated_neighbor_count"] == 20
    assert payload["stopping_reason"] == "no-improving-neighbor"
    assert payload["trace_rows"][0]["move_type"] == "tbr"
    assert payload["trace_rows"][0]["accepted_move"] is False
    assert payload["trace_rows"][1]["accepted_move"] is True
    assert payload["trace_rows"][1]["trace_reason"] == "accepted-improving-move"
    assert payload["trace_rows"][-1]["trace_reason"] == "no-improving-neighbor"
    assert payload["trace_rows"][1]["cut_edge_id"] == "A|B|C|D"
    assert payload["trace_rows"][1]["left_attachment_branch_id"] == "D"
    assert payload["trace_rows"][1]["right_attachment_branch_id"] == "interface"
    assert payload["trace_rows"][1]["optimized_branch_count"] == 8
