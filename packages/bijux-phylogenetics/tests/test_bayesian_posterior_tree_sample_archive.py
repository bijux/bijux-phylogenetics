from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaRunReport,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    run_fixed_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    build_bayesian_posterior_tree_sample,
    build_metropolis_hastings_posterior_tree_sample_archive,
    load_bayesian_posterior_tree_sample_archive,
    write_bayesian_posterior_tree_sample_archive,
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
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_metropolis_hastings_posterior_tree_sample_archive_round_trips_metadata(
    tmp_path: Path,
) -> None:
    run_report = _build_fixed_topology_dna_run_report(seed=13)

    archive = build_metropolis_hastings_posterior_tree_sample_archive(
        chain_report=run_report.chain_report,
    )
    archive_path = write_bayesian_posterior_tree_sample_archive(
        tmp_path / "fixed-topology-dna.posterior-tree-samples.json",
        archive,
    )
    loaded_archive = load_bayesian_posterior_tree_sample_archive(archive_path)

    assert loaded_archive == archive
    assert loaded_archive.sample_count == len(run_report.chain_report.sampled_states)
    first_sample = loaded_archive.samples[0]
    reference_state = run_report.chain_report.sampled_states[0]
    assert first_sample.sample_index == 0
    assert first_sample.iteration_index == 0
    assert first_sample.model_id == "substitution-model=HKY85"
    assert first_sample.posterior_log_score == reference_state.posterior_log_score
    assert first_sample.tree.topology_id == reference_state.tree.topology_id
    assert first_sample.tree.branch_rows == reference_state.tree.branch_rows


def test_bayesian_posterior_tree_sample_accepts_explicit_model_id_override() -> None:
    state = build_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"kappa": 2.0},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-1.0,
                parameter_values={"rate": 4.0},
            )
        ],
        log_likelihood=-2.5,
    )

    sample = build_bayesian_posterior_tree_sample(
        sample_index=0,
        state=state,
        model_id="manual-hky-like",
    )

    assert sample.model_id == "manual-hky-like"
    assert sample.posterior_log_score == -3.5


def test_bayesian_posterior_tree_sample_rejects_missing_model_identity() -> None:
    state = build_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"kappa": 2.0},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-1.0,
                parameter_values={"rate": 4.0},
            )
        ],
        log_likelihood=-2.5,
    )

    with pytest.raises(PhylogeneticsError) as error_info:
        build_bayesian_posterior_tree_sample(
            sample_index=0,
            state=state,
        )

    assert error_info.value.code == "bayesian_posterior_tree_sample_model_id_missing"


def test_load_bayesian_posterior_tree_sample_archive_rejects_sample_count_mismatch(
    tmp_path: Path,
) -> None:
    run_report = _build_fixed_topology_dna_run_report(seed=17)
    archive = build_metropolis_hastings_posterior_tree_sample_archive(
        chain_report=run_report.chain_report,
    )
    archive_path = write_bayesian_posterior_tree_sample_archive(
        tmp_path / "fixed-topology-dna.posterior-tree-samples.json",
        archive,
    )
    payload = json.loads(archive_path.read_text(encoding="utf-8"))
    payload["sample_count"] = payload["sample_count"] + 1
    archive_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(PhylogeneticsError) as error_info:
        load_bayesian_posterior_tree_sample_archive(archive_path)

    assert (
        error_info.value.code == "bayesian_posterior_tree_sample_archive_count_mismatch"
    )


def _build_fixed_topology_dna_run_report(*, seed: int) -> FixedTopologyDnaRunReport:
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
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
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.35,
        base_frequency_move_weight=1.0,
        base_frequency_coordinate_standard_deviation=0.45,
    )
    return run_fixed_topology_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")),
        records=load_fasta_alignment(
            fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
        ),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=9,
        sample_every=1,
        seed=seed,
        observation_policy="reject",
    )
