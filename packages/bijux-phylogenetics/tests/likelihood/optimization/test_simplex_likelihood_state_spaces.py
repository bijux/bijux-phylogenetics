from __future__ import annotations

from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    resolve_codon_state_space,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.discrete_observation_policies import (
    resolve_discrete_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    normalize_protein_likelihood_records,
    protein_leaf_likelihood_vector,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.simplex_search import (
    run_simplex_likelihood_search,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_simplex_search_keeps_protein_frequency_candidates_valid() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = normalize_protein_likelihood_records(
        load_fasta_alignment(
            fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
        ),
        model_name="protein simplex likelihood proof",
        gap_policy="reject",
        missing_policy="reject",
        ambiguity_policy="reject",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    candidate_vectors: list[dict[str, float]] = []
    initial_values = {
        state: 1.0 / len(PROTEIN_STATE_ORDER) for state in PROTEIN_STATE_ORDER
    }
    search_result = run_simplex_likelihood_search(
        initial_values_by_name=initial_values,
        expected_component_names=PROTEIN_STATE_ORDER,
        owner_name="protein simplex likelihood proof",
        evaluate=lambda candidate: _evaluate_protein_frequency_candidate(
            tree,
            compressed_patterns,
            candidate,
            candidate_vectors,
        ),
        improvement_tolerance=1e-6,
        max_coordinate_passes=1,
    )

    assert search_result.function_evaluation_count > 1
    assert search_result.optimization_pass_count >= 1
    assert numpy.isfinite(search_result.objective_value)
    assert candidate_vectors
    for candidate in candidate_vectors:
        _assert_valid_simplex_candidate(candidate)
    _assert_valid_simplex_candidate(
        search_result.optimized_parameterization.constrained_mapping()
    )


def test_simplex_search_keeps_codon_frequency_candidates_valid() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    state_space = resolve_codon_state_space()
    initial_values = {
        codon: 1.0 / len(state_space.state_order) for codon in state_space.state_order
    }
    candidate_vectors: list[dict[str, float]] = []
    search_result = run_simplex_likelihood_search(
        initial_values_by_name=initial_values,
        expected_component_names=state_space.state_order,
        owner_name="codon simplex likelihood proof",
        evaluate=lambda candidate: _evaluate_codon_frequency_candidate(
            tree,
            records,
            candidate,
            candidate_vectors,
        ),
        improvement_tolerance=1e-6,
        max_coordinate_passes=1,
    )

    assert search_result.function_evaluation_count > 1
    assert search_result.optimization_pass_count >= 1
    assert numpy.isfinite(search_result.objective_value)
    assert candidate_vectors
    for candidate in candidate_vectors:
        _assert_valid_simplex_candidate(candidate)
    _assert_valid_simplex_candidate(
        search_result.optimized_parameterization.constrained_mapping()
    )


def test_simplex_search_keeps_generic_state_frequency_candidates_valid() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0"),
        AlignmentRecord(identifier="B", sequence="2"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(
        records,
        site_tokenizer=lambda sequence: (sequence,),
    )
    state_order = ("0", "1", "2", "3")
    initial_values = dict.fromkeys(state_order, 0.25)
    candidate_vectors: list[dict[str, float]] = []
    search_result = run_simplex_likelihood_search(
        initial_values_by_name=initial_values,
        expected_component_names=state_order,
        owner_name="finite-state simplex likelihood proof",
        evaluate=lambda candidate: _evaluate_generic_state_frequency_candidate(
            tree,
            compressed_patterns,
            state_order,
            candidate,
            candidate_vectors,
        ),
        improvement_tolerance=1e-6,
        max_coordinate_passes=2,
    )

    assert search_result.function_evaluation_count > 1
    assert search_result.optimization_pass_count >= 1
    assert numpy.isfinite(search_result.objective_value)
    assert candidate_vectors
    for candidate in candidate_vectors:
        _assert_valid_simplex_candidate(candidate)
    _assert_valid_simplex_candidate(
        search_result.optimized_parameterization.constrained_mapping()
    )


def _evaluate_protein_frequency_candidate(
    tree: PhyloTree,
    compressed_patterns,
    candidate_frequencies: dict[str, float],
    candidate_vectors: list[dict[str, float]],
) -> tuple[dict[str, float], float]:
    candidate_vectors.append(dict(candidate_frequencies))
    _assert_valid_simplex_candidate(candidate_frequencies)
    surface = _build_protein_frequency_surface(candidate_frequencies)
    context = compile_finite_state_likelihood_context(tree, surface)
    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    return candidate_frequencies, log_likelihood


def _evaluate_codon_frequency_candidate(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    candidate_frequencies: dict[str, float],
    candidate_vectors: list[dict[str, float]],
) -> tuple[object, float]:
    candidate_vectors.append(dict(candidate_frequencies))
    _assert_valid_simplex_candidate(candidate_frequencies)
    report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        codon_frequencies=candidate_frequencies,
    )
    return report, report.log_likelihood


def _evaluate_generic_state_frequency_candidate(
    tree: PhyloTree,
    compressed_patterns,
    state_order: tuple[str, ...],
    candidate_frequencies: dict[str, float],
    candidate_vectors: list[dict[str, float]],
) -> tuple[dict[str, float], float]:
    candidate_vectors.append(dict(candidate_frequencies))
    _assert_valid_simplex_candidate(candidate_frequencies)
    surface = _build_generic_state_frequency_surface(
        state_order=state_order,
        stationary_frequencies=candidate_frequencies,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    return candidate_frequencies, log_likelihood


def _build_protein_frequency_surface(
    candidate_frequencies: dict[str, float],
) -> FiniteStateLikelihoodSurface:
    stationary_vector = numpy.asarray(
        [candidate_frequencies[state] for state in PROTEIN_STATE_ORDER],
        dtype=float,
    )
    rate_matrix = _build_equal_input_rate_matrix(
        state_labels=PROTEIN_STATE_ORDER,
        stationary_vector=stationary_vector,
    )
    return FiniteStateLikelihoodSurface(
        model_name="protein simplex likelihood proof",
        state_labels=PROTEIN_STATE_ORDER,
        root_prior=stationary_vector,
        root_prior_source="stationary",
        observation_policy="reject",
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="protein simplex likelihood proof",
                node_name=taxon_name,
                gap_policy="reject",
                missing_policy="reject",
                ambiguity_policy="reject",
            )
        ),
    )


def _build_generic_state_frequency_surface(
    *,
    state_order: tuple[str, ...],
    stationary_frequencies: dict[str, float],
) -> FiniteStateLikelihoodSurface:
    stationary_vector = numpy.asarray(
        [stationary_frequencies[state] for state in state_order],
        dtype=float,
    )
    rate_matrix = _build_equal_input_rate_matrix(
        state_labels=state_order,
        stationary_vector=stationary_vector,
    )
    return FiniteStateLikelihoodSurface(
        model_name="finite-state simplex likelihood proof",
        state_labels=state_order,
        root_prior=stationary_vector,
        root_prior_source="stationary",
        observation_policy="reject",
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_discrete_observation_leaf_vector(
                observed_state,
                model_name="finite-state simplex likelihood proof",
                node_name=taxon_name,
                state_order=list(state_order),
                observation_policy="reject",
            )
        ),
    )


def _build_equal_input_rate_matrix(
    *,
    state_labels: tuple[str, ...],
    stationary_vector: numpy.ndarray,
) -> numpy.ndarray:
    rate_matrix = numpy.tile(stationary_vector, (len(state_labels), 1))
    numpy.fill_diagonal(rate_matrix, 0.0)
    numpy.fill_diagonal(rate_matrix, -rate_matrix.sum(axis=1))
    return normalize_ctmc_rate_matrix_by_expected_substitution_rate(
        rate_matrix,
        stationary_vector,
        state_labels=state_labels,
    )


def _assert_valid_simplex_candidate(candidate: dict[str, float]) -> None:
    assert candidate
    assert min(candidate.values()) > 0.0
    assert pytest.approx(1.0, abs=1e-12, rel=0.0) == sum(candidate.values())
