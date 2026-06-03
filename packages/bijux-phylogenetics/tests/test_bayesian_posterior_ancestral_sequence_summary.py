from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    run_fixed_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
    run_joint_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    build_posterior_ancestral_sequence_definition,
    summarize_nucleotide_posterior_ancestral_sequences,
    write_posterior_ancestral_sequence_fasta,
    write_posterior_ancestral_state_probability_table,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_posterior_ancestral_sequence_summary_aggregates_fixed_topology_hky85_chain(
    tmp_path: Path,
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
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    run_report = run_fixed_topology_dna_metropolis_hastings(
        tree=load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")),
        records=records,
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=9,
        sample_every=1,
        seed=0,
    )

    summary = summarize_nucleotide_posterior_ancestral_sequences(
        run_report.chain_report.sampled_states,
        definition=build_posterior_ancestral_sequence_definition(
            records=records,
            posterior_probability_threshold=0.6,
            minimum_clade_posterior_support=0.5,
        ),
    )

    assert summary.sample_count == len(run_report.chain_report.sampled_states)
    assert summary.distinct_topology_count == 1
    assert summary.sampled_substitution_models == ["HKY85"]
    assert summary.sequence_records
    assert all(
        len(record.sequence) == len(records[0].sequence)
        for record in summary.sequence_records
    )
    _assert_site_probabilities_sum_to_clade_support(summary)

    fasta_path = write_posterior_ancestral_sequence_fasta(
        tmp_path / "posterior-ancestral-sequences.fasta",
        summary,
    )
    table_path = write_posterior_ancestral_state_probability_table(
        tmp_path / "posterior-ancestral-probabilities.tsv",
        summary,
    )

    fasta_text = fasta_path.read_text(encoding="utf-8")
    table_lines = table_path.read_text(encoding="utf-8").splitlines()
    assert fasta_text.startswith(">")
    assert "clade_id\trepresentative_node_id" in table_lines[0]
    assert len(table_lines) == len(summary.state_probability_rows) + 1


def test_posterior_ancestral_sequence_summary_tracks_topology_uncertainty_across_joint_chain() -> (
    None
):
    sequence_model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="K80",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            )
        ),
    )
    joint_model_definition = build_joint_topology_dna_model_definition(
        sequence_model_definition=sequence_model_definition,
        topology_prior=build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"]),
    )
    joint_proposal_schedule = build_joint_topology_dna_proposal_schedule(
        sequence_proposal_schedule=build_fixed_topology_dna_proposal_schedule(
            model_definition=sequence_model_definition,
            branch_length_move_weight=1.0,
            branch_length_log_scale_standard_deviation=0.25,
            kappa_move_weight=1.0,
            kappa_log_scale_standard_deviation=0.35,
        ),
        nni_move_weight=1.0,
    )
    records = _build_ambiguous_alignment_records()
    run_report = run_joint_topology_dna_metropolis_hastings(
        tree=_build_ambiguous_start_tree(),
        records=records,
        model_definition=joint_model_definition,
        proposal_schedule=joint_proposal_schedule,
        iteration_count=20,
        sample_every=1,
        seed=0,
    )

    summary = summarize_nucleotide_posterior_ancestral_sequences(
        run_report.chain_report.sampled_states,
        definition=build_posterior_ancestral_sequence_definition(
            records=records,
            posterior_probability_threshold=0.4,
            minimum_clade_posterior_support=0.4,
        ),
    )

    assert summary.distinct_topology_count > 1
    assert summary.sampled_substitution_models == ["K80"]
    assert any("multiple sampled topologies" in warning for warning in summary.warnings)
    assert any(
        row.clade_posterior_probability < 1.0
        for row in summary.state_probability_rows
        if len(row.descendant_taxa) < len(records)
    )
    assert any(
        row.clade_posterior_probability < 1.0
        for row in summary.sequence_records
        if len(row.descendant_taxa) < len(records)
    )
    _assert_site_probabilities_sum_to_clade_support(summary)


def _assert_site_probabilities_sum_to_clade_support(
    summary,
) -> None:
    grouped_probabilities: dict[tuple[str, int], list[float]] = {}
    clade_support_by_key: dict[tuple[str, int], float] = {}
    for row in summary.state_probability_rows:
        key = (row.clade_id, row.site_position)
        grouped_probabilities.setdefault(key, []).append(
            row.marginal_posterior_probability
        )
        clade_support_by_key.setdefault(key, row.clade_posterior_probability)
    for key, probabilities in grouped_probabilities.items():
        assert math.isclose(
            sum(probabilities),
            clade_support_by_key[key],
            rel_tol=0.0,
            abs_tol=1e-9,
        )


def _build_ambiguous_start_tree() -> PhyloTree:
    tree = PhyloTree.from_newick("(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);")
    tree.rooted = True
    return tree


def _build_ambiguous_alignment_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="A", sequence="A"),
        AlignmentRecord(identifier="B", sequence="C"),
        AlignmentRecord(identifier="C", sequence="G"),
        AlignmentRecord(identifier="D", sequence="T"),
    ]
