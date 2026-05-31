from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    build_likelihood_stepwise_addition_tree_from_alignment,
    validate_likelihood_stepwise_addition_model,
)
from bijux_phylogenetics.phylo.topology import write_stepwise_addition_artifacts

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_gateway_exports_stepwise_addition_surface() -> None:
    assert (
        likelihood_api.build_likelihood_stepwise_addition_tree_from_alignment
        is build_likelihood_stepwise_addition_tree_from_alignment
    )
    assert (
        likelihood_api.validate_likelihood_stepwise_addition_model
        is validate_likelihood_stepwise_addition_model
    )


def test_build_likelihood_stepwise_addition_tree_reports_candidate_likelihoods_per_edge() -> (
    None
):
    tree, report = build_likelihood_stepwise_addition_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "greedy-stepwise-addition-tree"
    assert report.objective_name == "likelihood-jc69"
    assert report.objective_direction == "maximize"
    assert report.tip_count == 4
    assert report.strictly_bifurcating is True
    assert report.all_requested_taxa_present_once is True
    assert len(report.trace_rows) == 2

    first_step = report.trace_rows[0]
    assert first_step.taxon == "C"
    assert first_step.best_edge_id == "root"
    assert len(first_step.tested_edge_rows) == 3
    assert [row.branch_id for row in first_step.tested_edge_rows] == ["root", "A", "B"]
    assert math.isclose(first_step.best_score, -34.13524969213178, abs_tol=1e-12)
    assert first_step.best_score > first_step.tested_edge_rows[1].score

    second_step = report.trace_rows[1]
    assert second_step.taxon == "D"
    assert second_step.best_edge_id == "C"
    assert len(second_step.tested_edge_rows) == 5
    assert [row.branch_id for row in second_step.tested_edge_rows] == [
        "root",
        "A",
        "B",
        "C",
        "A|B",
    ]
    assert math.isclose(second_step.best_score, -34.13524969797671, abs_tol=1e-12)
    assert second_step.best_score > second_step.tested_edge_rows[0].score

    assert tree.to_newick() == (
        "((A:2.43539016857288e-10,B:2.43539016857288e-10):0.999999999756461,"
        "(C:2.43539016857288e-10,D:2.43539016857288e-10):0.999999999756461);"
    )
    assert math.isclose(report.final_score, -34.13524969797671, abs_tol=1e-12)


def test_build_likelihood_stepwise_addition_tree_rejects_insertion_order_taxa_mismatch() -> (
    None
):
    with pytest.raises(ValueError) as error_info:
        build_likelihood_stepwise_addition_tree_from_alignment(
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            insertion_order=["A", "B", "C", "X"],
        )

    assert (
        str(error_info.value)
        == "likelihood stepwise addition insertion_order must match alignment taxa exactly"
    )


def test_write_stepwise_addition_artifacts_materializes_likelihood_outputs(
    tmp_path: Path,
) -> None:
    _tree, report = build_likelihood_stepwise_addition_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        upper_branch_length_bound=1.0,
    )

    outputs = write_stepwise_addition_artifacts(
        tmp_path / "likelihood-stepwise-addition-run",
        report,
    )

    assert set(outputs) == {"tree_path", "trace_path", "run_json_path"}
    assert outputs["trace_path"].read_text(encoding="utf-8").startswith(
        "step_index\ttaxon\tinserted_taxa\ttested_edge_id\ttested_edge_descendant_taxa\ttested_edge_score\tbest_edge_id\tbest_edge_descendant_taxa\tbest_score\tselected\tcandidate_tree_newick\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["objective_name"] == "likelihood-jc69"
    assert payload["objective_direction"] == "maximize"
    assert payload["tip_count"] == 4
    assert payload["trace_rows"][1]["best_edge_id"] == "C"


def test_validate_likelihood_stepwise_addition_model_rejects_unsupported_model() -> None:
    with pytest.raises(ValueError) as error_info:
        validate_likelihood_stepwise_addition_model("k80")

    assert (
        str(error_info.value) == "likelihood stepwise addition model must be one of jc69"
    )
