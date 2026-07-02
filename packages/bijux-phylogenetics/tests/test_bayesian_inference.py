from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian import inference as bayesian_inference
from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_run_bayesian_inference_dispatches_fixed_topology_dna_runner(
    monkeypatch,
) -> None:
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
    tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    captured: dict[str, object] = {}
    sentinel_report = object()

    def fake_runner(**kwargs: object) -> object:
        captured.update(kwargs)
        return sentinel_report

    monkeypatch.setattr(
        bayesian_inference,
        "run_fixed_topology_dna_metropolis_hastings",
        fake_runner,
    )

    report = bayesian_inference.run_bayesian_inference(
        tree=tree,
        records=records,
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=15,
        sample_every=3,
        seed=11,
        observation_policy="reject",
    )

    assert report is sentinel_report
    assert captured == {
        "tree": tree,
        "records": records,
        "model_definition": model_definition,
        "proposal_schedule": proposal_schedule,
        "iteration_count": 15,
        "sample_every": 3,
        "seed": 11,
        "observation_policy": "reject",
    }
