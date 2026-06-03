from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaRunReport,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.posterior_predictive_p_values import (
    summarize_posterior_predictive_p_values,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveAlignmentSimulationReport,
    PosteriorPredictiveObservedStatisticRow,
    PosteriorPredictiveReplicateStatisticRow,
    build_posterior_predictive_simulation_definition,
    simulate_fixed_topology_dna_posterior_predictive,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def test_posterior_predictive_p_values_summarize_alignment_tail_areas() -> None:
    report = PosteriorPredictiveAlignmentSimulationReport(
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=4,
            seed=0,
        ),
        model_name="HKY85",
        taxon_count=2,
        alignment_length=12,
        observed_statistic_rows=[
            PosteriorPredictiveObservedStatisticRow(
                statistic_name="gc-fraction",
                value=0.95,
            )
        ],
        replicate_statistic_rows=[
            PosteriorPredictiveReplicateStatisticRow(
                statistic_name="gc-fraction",
                replicate_index=index,
                posterior_sample_index=0,
                posterior_iteration_index=index,
                value=value,
            )
            for index, value in enumerate((0.10, 0.20, 0.30, 0.40))
        ],
        statistic_summary_rows=[],
        replicates=[],
    )

    p_value_report = summarize_posterior_predictive_p_values(report)

    assert p_value_report.report_kind == "alignment"
    assert p_value_report.model_name == "HKY85"
    assert p_value_report.statistic_count == 1
    assert (
        p_value_report.p_value_method
        == "posterior-predictive tail area with plus-one Monte Carlo smoothing"
    )

    row = p_value_report.statistic_rows[0]
    assert row.statistic_name == "gc-fraction"
    assert row.lower_tail_count == 4
    assert row.upper_tail_count == 0
    assert row.lower_tail_probability == pytest.approx(1.0, abs=1e-12)
    assert row.upper_tail_probability == pytest.approx(0.2, abs=1e-12)
    assert row.posterior_predictive_p_value == pytest.approx(0.4, abs=1e-12)
    assert row.replicate_mean == pytest.approx(0.25, abs=1e-12)
    assert row.replicate_median == pytest.approx(0.25, abs=1e-12)
    assert row.replicate_minimum == pytest.approx(0.10, abs=1e-12)
    assert row.replicate_maximum == pytest.approx(0.40, abs=1e-12)


def test_fixed_topology_dna_posterior_predictive_p_values_flag_extreme_gc_mismatch() -> (
    None
):
    observed_records = [
        AlignmentRecord(identifier="A", sequence="A" * 400),
        AlignmentRecord(identifier="B", sequence="A" * 400),
    ]
    simulation_report = simulate_fixed_topology_dna_posterior_predictive(
        run_report=_build_gc_heavy_fixed_topology_dna_run_report(),
        records=observed_records,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=39,
            seed=7,
        ),
    )

    p_value_report = summarize_posterior_predictive_p_values(simulation_report)

    gc_fraction_row = next(
        row
        for row in p_value_report.statistic_rows
        if row.statistic_name == "gc-fraction"
    )

    assert gc_fraction_row.observed_value == pytest.approx(0.0, abs=1e-12)
    assert gc_fraction_row.lower_tail_count == 0
    assert gc_fraction_row.lower_tail_probability == pytest.approx(0.025, abs=1e-12)
    assert gc_fraction_row.upper_tail_probability == pytest.approx(1.0, abs=1e-12)
    assert gc_fraction_row.posterior_predictive_p_value == pytest.approx(
        0.05,
        abs=1e-12,
    )
    assert gc_fraction_row.replicate_minimum > 0.85


def _build_gc_heavy_fixed_topology_dna_run_report() -> FixedTopologyDnaRunReport:
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=3.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.0
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
            ),
        ),
    )
    proposal_schedule = build_fixed_topology_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.2,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.2,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.2,
    )
    tree = PhyloTree.from_newick("(A:0.01,B:0.01);")
    tree.rooted = True
    sampled_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "HKY85"},
            scalar_parameters={"kappa": 2.0},
            vector_parameters={
                "base-frequencies": {
                    "A": 0.01,
                    "C": 0.49,
                    "G": 0.49,
                    "T": 0.01,
                }
            },
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.2,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:kappa",
                family="exponential",
                log_prior=-0.1,
            ),
            build_bayesian_prior_component_state(
                component_name="substitution:base-frequencies",
                family="dirichlet",
                log_prior=-0.1,
            ),
        ],
        log_likelihood=-1.0,
    )
    return FixedTopologyDnaRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        observation_policy="reject",
        chain_report=MetropolisHastingsRunReport(
            iteration_count=1,
            sample_every=1,
            seed=0,
            accepted_count=1,
            rejected_count=0,
            acceptance_rate=1.0,
            initial_state=sampled_state,
            final_state=sampled_state,
            sampled_states=[sampled_state],
            step_rows=[],
        ),
        posterior_rows=[],
    )
