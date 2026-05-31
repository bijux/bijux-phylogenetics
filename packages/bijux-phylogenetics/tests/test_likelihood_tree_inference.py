from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.likelihood import (
    infer_nucleotide_likelihood_tree_from_alignment,
    write_nucleotide_likelihood_tree_inference_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_tree_inference_surface() -> None:
    assert (
        likelihood_api.infer_nucleotide_likelihood_tree_from_alignment
        is infer_nucleotide_likelihood_tree_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tree_inference_artifacts
        is write_nucleotide_likelihood_tree_inference_artifacts
    )


def test_native_tree_inference_runs_stepwise_model_selection_and_search() -> None:
    report = infer_nucleotide_likelihood_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-tree-inference"
    assert report.taxon_count == 4
    assert report.site_count == 12
    assert report.pattern_count == 2
    assert report.stepwise_addition_model_name == "JC69"
    assert report.model_selection_strategy == "model-selection"
    assert report.model_selection_criterion == "aic"
    assert report.search_method == "nni"
    assert report.start_tree_source_policy == "stepwise-addition-tree-plus-random-tree"
    assert report.random_start_tree_count == 2
    assert report.start_tree_seed == 17
    assert len(report.run_summaries) == 3
    assert report.selected_model_name in {"JC69", "K80", "F81", "HKY85", "GTR"}
    assert len(report.model_selection_report.rows) == 5
    assert report.best_run_source_label == "stepwise-addition-tree"
    assert report.best_search_report.algorithm == "nucleotide-likelihood-nni-search"
    assert report.best_search_report.final_tree_newick == report.best_final_tree_newick
    assert math.isclose(
        report.best_final_log_likelihood,
        -27.213282538844922,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [row.start_tree_source_label for row in report.run_summaries] == [
        "stepwise-addition-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row.best_run for row in report.run_summaries] == [True, False, False]


def test_native_tree_inference_fixed_model_writes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = infer_nucleotide_likelihood_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_tree_inference_artifacts(
        tmp_path / "likelihood-tree-inference",
        report,
    )

    assert report.model_selection_strategy == "fixed-model"
    assert report.model_selection_criterion is None
    assert report.selected_model_name == "JC69"
    assert len(report.model_selection_report.rows) == 1
    assert set(outputs) == {
        "stepwise_tree_path",
        "start_tree_path",
        "final_tree_path",
        "best_tree_set_path",
        "likelihood_table_path",
        "model_table_path",
        "search_trace_path",
        "run_json_path",
        "candidate_table_path",
    }
    assert outputs["likelihood_table_path"].read_text(encoding="utf-8").startswith(
        "start_tree_source_kind\tstart_tree_source_label\tstart_tree_generation_seed\tsearch_algorithm\tstart_log_likelihood\tfinal_log_likelihood\tfinal_likelihood_rank\tfinal_topology_fingerprint\tsearch_iteration_count\taccepted_move_count\tevaluated_neighbor_count\tbranch_reoptimization_policy\tsubstitution_parameter_policy\tstopping_reason\tbest_run\tstart_tree_newick\tfinal_tree_newick\n"
    )
    assert outputs["model_table_path"].read_text(encoding="utf-8").startswith(
        "model_name\tbase_model_name\trate_heterogeneity_model\tfit_succeeded\tparameter_count\tlog_likelihood\taic\taicc\tbic\tdelta_aic\takaike_weight\trank\tselected_by_aic\tselected_by_aicc\tselected_by_bic\tparameter_values\twarnings\n"
    )
    assert outputs["search_trace_path"].read_text(encoding="utf-8").startswith(
        "event_index\tevent_kind\titeration\tmove_type\tcandidate_topology_fingerprint\tlog_likelihood_before\tlog_likelihood_after\tlog_likelihood_delta\taccepted_move\ttrace_reason\ttree_before_newick\ttree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tbranch_reoptimization_policy\tbranch_reoptimization_scope\toptimized_branch_count\toptimized_branch_clade_ids\tbranch_reoptimization_converged\tbranch_optimization_pass_count\tbranch_function_evaluation_count\tboundary_warning_messages\tstopping_reason\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-tree-inference"
    assert payload["selected_model_name"] == "JC69"
    assert payload["best_run_source_label"] == "random-tree-seed-17"
    assert len(payload["run_summaries"]) == 3
    assert len(load_newick_tree_set(outputs["best_tree_set_path"])) >= 1
