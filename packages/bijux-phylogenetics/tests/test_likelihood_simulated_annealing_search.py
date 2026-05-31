from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_simulated_annealing,
    search_nucleotide_likelihood_simulated_annealing_from_alignment,
    write_nucleotide_likelihood_simulated_annealing_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_simulated_annealing_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_simulated_annealing_from_alignment
        is search_nucleotide_likelihood_simulated_annealing_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_simulated_annealing_artifacts
        is write_nucleotide_likelihood_simulated_annealing_artifacts
    )


def test_likelihood_simulated_annealing_accepts_worse_nni_move_by_temperature() -> None:
    report = search_nucleotide_likelihood_simulated_annealing(
        loads_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        proposal_move_family="nni",
        annealing_seed=1,
        iteration_count=4,
        initial_temperature=10.0,
        cooling_rate=0.8,
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-simulated-annealing-search"
    assert report.proposal_move_family == "nni"
    assert report.iteration_count_completed == 4
    assert report.accepted_move_count == 2
    assert report.rejected_move_count == 2
    assert report.accepted_worse_move_count == 1
    assert report.stopping_reason == "temperature-schedule-exhausted"
    assert math.isclose(
        report.start_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.final_log_likelihood,
        -36.94280000264263,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.best_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.best_tree_newick == report.start_tree_newick
    assert report.final_tree_newick != report.best_tree_newick

    first_row, second_row, third_row, fourth_row = report.trace_rows

    assert first_row.acceptance_decision == "accepted-worse-move"
    assert first_row.accepted_move is True
    assert first_row.best_tree_improved is False
    assert math.isclose(first_row.temperature, 10.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        first_row.log_likelihood_delta,
        -2.8075503046659165,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.acceptance_probability,
        0.7552133170757632,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.acceptance_uniform_draw,
        0.5692038748222122,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    assert second_row.acceptance_decision == "accepted-equal-move"
    assert second_row.accepted_move is True
    assert math.isclose(
        second_row.log_likelihood_delta,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        second_row.acceptance_probability,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    assert third_row.acceptance_decision == "rejected-worse-move"
    assert third_row.accepted_move is False
    assert math.isclose(
        third_row.log_likelihood_delta,
        -17.499717347003227,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        third_row.acceptance_probability,
        0.0649374464040797,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        third_row.acceptance_uniform_draw,
        0.7609624449125756,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    assert fourth_row.acceptance_decision == "rejected-worse-move"
    assert fourth_row.accepted_move is False
    assert math.isclose(
        fourth_row.temperature,
        5.12,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_likelihood_simulated_annealing_supports_rooted_spr_proposals() -> None:
    report = search_nucleotide_likelihood_simulated_annealing_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        proposal_move_family="spr",
        annealing_seed=3,
        iteration_count=3,
        initial_temperature=5.0,
        cooling_rate=0.7,
        upper_branch_length_bound=1.0,
    )

    assert report.proposal_move_family == "spr"
    assert report.accepted_move_count == 3
    assert report.rejected_move_count == 0
    assert report.accepted_worse_move_count == 1
    assert math.isclose(
        report.start_log_likelihood,
        -8.174114933271198,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.final_log_likelihood,
        -8.642039982465846,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.best_log_likelihood,
        -7.742402023371725,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [
        (row.pruned_clade_id, row.regraft_target_branch_id, row.acceptance_decision)
        for row in report.trace_rows
    ] == [
        ("B", "C", "accepted-improving-move"),
        ("A", "B|C", "accepted-improving-move"),
        ("E", "B", "accepted-worse-move"),
    ]


def test_write_likelihood_simulated_annealing_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_simulated_annealing(
        loads_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        proposal_move_family="nni",
        annealing_seed=1,
        iteration_count=4,
        initial_temperature=10.0,
        cooling_rate=0.8,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_simulated_annealing_artifacts(
        tmp_path / "likelihood-simulated-annealing-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "start_tree_path",
        "final_tree_path",
        "best_tree_path",
        "trace_table_path",
        "run_json_path",
    }
    assert outputs["trace_table_path"].read_text(encoding="utf-8").startswith(
        "iteration\ttemperature\tmove_family\tcurrent_log_likelihood_before\tproposed_log_likelihood\tlog_likelihood_delta\tacceptance_probability\tacceptance_uniform_draw\tacceptance_decision\taccepted_move\tbest_tree_improved\tcurrent_tree_before_newick\tproposed_tree_newick\tcurrent_tree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tpruned_clade_id\tregraft_target_branch_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-simulated-annealing-search"
    assert payload["proposal_move_family"] == "nni"
    assert payload["annealing_seed"] == 1
    assert payload["iteration_count_completed"] == 4
    assert payload["accepted_worse_move_count"] == 1
    assert payload["trace_rows"][0]["acceptance_decision"] == "accepted-worse-move"
    assert payload["trace_rows"][2]["accepted_move"] is False
