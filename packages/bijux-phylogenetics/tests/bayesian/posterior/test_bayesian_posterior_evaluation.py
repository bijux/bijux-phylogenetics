from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.posterior_evaluation import (
    BayesianPosteriorDecompositionRow,
    BayesianPosteriorEvaluationReport,
    evaluate_bayesian_posterior,
    summarize_bayesian_posterior_decomposition,
    write_bayesian_posterior_decomposition_table,
    write_bayesian_posterior_evaluation_artifacts,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_evaluate_bayesian_posterior_combines_likelihood_and_prior_into_one_report() -> (
    None
):
    report = evaluate_bayesian_posterior(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"kappa": 2.0}
        ),
        update_prior_components=_build_prior_components,
        update_log_likelihood=lambda _state: -4.5,
    )

    assert report.total_log_prior == pytest.approx(-2.0)
    assert report.log_likelihood == pytest.approx(-4.5)
    assert report.hastings_independent_log_posterior == pytest.approx(-6.5)
    assert report.state.posterior_log_score == pytest.approx(-6.5)
    assert report.state.random_state is None
    assert report.decomposition_valid is True


def test_summarize_bayesian_posterior_decomposition_exposes_likelihood_and_prior_terms() -> (
    None
):
    state = build_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="tree-topology",
                family="uniform-rooted-labeled-bifurcating",
                log_prior=-1.25,
            ),
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.75,
                parameter_values={"rate": 1.0},
            ),
        ],
        log_likelihood=-4.5,
    )

    report = summarize_bayesian_posterior_decomposition(state)

    assert isinstance(report, BayesianPosteriorEvaluationReport)
    assert report.analysis_name == "generic-bayesian-state"
    assert report.posterior_term_source == "derived_from_likelihood_and_prior"
    assert report.likelihood_term_source == "evaluated"
    assert report.prior_term_source == "summed_prior_components"
    assert report.log_likelihood == pytest.approx(-4.5)
    assert report.total_log_prior == pytest.approx(-2.0)
    assert report.hastings_independent_log_posterior == pytest.approx(-6.5)
    assert report.decomposition_delta == pytest.approx(0.0)
    assert report.decomposition_valid is True
    assert [row.term_kind for row in report.rows] == ["likelihood", "prior", "prior"]
    assert report.rows[1] == BayesianPosteriorDecompositionRow(
        term_kind="prior",
        term_name="tree-topology",
        family="uniform-rooted-labeled-bifurcating",
        log_value=-1.25,
        parameter_values={},
    )
    assert report.rows[2].parameter_values == {"rate": 1.0}


def test_summarize_bayesian_posterior_decomposition_rejects_invalid_identity_tolerance() -> (
    None
):
    state = build_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.75,
            )
        ],
        log_likelihood=-4.5,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="identity_tolerance requires one nonnegative finite value",
    ):
        summarize_bayesian_posterior_decomposition(state, identity_tolerance=-1.0)


def test_write_bayesian_posterior_decomposition_table_writes_verifiable_identity_columns(
    tmp_path: Path,
) -> None:
    report = evaluate_bayesian_posterior(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_build_prior_components,
        update_log_likelihood=lambda _state: -4.5,
    )
    output_path = tmp_path / "posterior-decomposition.tsv"

    write_bayesian_posterior_decomposition_table(output_path, report)

    text = output_path.read_text(encoding="utf-8")
    assert text.startswith(
        "term_kind\tterm_name\tfamily\tlog_value\tparameter_values\tanalysis_name\ttopology_id\thastings_independent_log_posterior\tlog_likelihood\ttotal_log_prior\tdecomposition_delta\tdecomposition_valid\tposterior_term_source\tlikelihood_term_source\tprior_term_source\tidentity_tolerance\n"
    )
    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["term_kind"] == "likelihood"
    assert rows[0]["decomposition_valid"] == "true"


def test_write_bayesian_posterior_evaluation_artifacts_emits_summary_bundle(
    tmp_path: Path,
) -> None:
    report = evaluate_bayesian_posterior(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_build_prior_components,
        update_log_likelihood=lambda _state: -4.5,
    )

    outputs = write_bayesian_posterior_evaluation_artifacts(tmp_path, report)

    assert outputs["decomposition_table_path"].name == "posterior-decomposition.tsv"
    assert outputs["summary_path"].name == "posterior-evaluation.json"
    summary_payload = json.loads(outputs["summary_path"].read_text(encoding="utf-8"))
    assert summary_payload["decomposition_valid"] is True
    assert summary_payload["hastings_independent_log_posterior"] == pytest.approx(-6.5)
    assert summary_payload["state"]["posterior_log_score"] == pytest.approx(-6.5)


def _build_prior_components(
    _state,
):
    return [
        build_bayesian_prior_component_state(
            component_name="tree-topology",
            family="uniform-rooted-labeled-bifurcating",
            log_prior=-1.25,
        ),
        build_bayesian_prior_component_state(
            component_name="branch-lengths",
            family="exponential",
            log_prior=-0.75,
            parameter_values={"rate": 1.0},
        ),
    ]
