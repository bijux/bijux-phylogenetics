from __future__ import annotations

from dataclasses import dataclass
import random

import numpy

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import DNA_STATE_ORDER
from bijux_phylogenetics.phylo.likelihood.f81 import f81_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import jc69_transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.k80 import k80_transition_probability_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

_UNIFORM_DNA_FREQUENCIES = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25,
}


@dataclass(frozen=True, slots=True)
class _GovernedNativeMaximumLikelihoodAccuracyCase:
    case_id: str
    model_name: str
    tree_newick: str
    sequence_length: int
    simulation_seed: int
    parameter_truth: dict[str, float]
    parameter_tolerances: dict[str, float]


_GOVERNED_NATIVE_MAXIMUM_LIKELIHOOD_ACCURACY_CASES: tuple[
    _GovernedNativeMaximumLikelihoodAccuracyCase,
    ...,
] = (
    _GovernedNativeMaximumLikelihoodAccuracyCase(
        case_id="k80-four-taxon-transition-bias",
        model_name="k80",
        tree_newick="((A:0.06,B:0.06):0.35,(C:0.06,D:0.06):0.33);",
        sequence_length=16,
        simulation_seed=203,
        parameter_truth={"kappa": 4.0},
        parameter_tolerances={"kappa": 2.0},
    ),
    _GovernedNativeMaximumLikelihoodAccuracyCase(
        case_id="f81-four-taxon-composition-bias",
        model_name="f81",
        tree_newick="((A:0.06,B:0.06):0.35,(C:0.06,D:0.06):0.33);",
        sequence_length=16,
        simulation_seed=307,
        parameter_truth={
            "base_frequency_a": 0.37,
            "base_frequency_c": 0.18,
            "base_frequency_g": 0.17,
            "base_frequency_t": 0.28,
        },
        parameter_tolerances={
            "base_frequency_a": 0.2,
            "base_frequency_c": 0.2,
            "base_frequency_g": 0.2,
            "base_frequency_t": 0.2,
        },
    ),
)


def governed_native_maximum_likelihood_accuracy_cases() -> tuple[
    _GovernedNativeMaximumLikelihoodAccuracyCase,
    ...,
]:
    return _GOVERNED_NATIVE_MAXIMUM_LIKELIHOOD_ACCURACY_CASES


def build_native_maximum_likelihood_accuracy_truth_tree(
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
) -> PhyloTree:
    return loads_newick(case.tree_newick)


def simulate_native_maximum_likelihood_accuracy_alignment(
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
) -> list[AlignmentRecord]:
    rng = random.Random(case.simulation_seed)  # nosec B311
    tree = build_native_maximum_likelihood_accuracy_truth_tree(case)
    root_sequence = [
        _sample_state(
            _root_state_probabilities(case),
            rng,
        )
        for _ in range(case.sequence_length)
    ]
    records_by_taxon: dict[str, str] = {}
    _propagate_alignment(
        tree.root,
        parent_sequence=root_sequence,
        case=case,
        rng=rng,
        records_by_taxon=records_by_taxon,
    )
    return [
        AlignmentRecord(identifier=taxon, sequence=sequence)
        for taxon, sequence in sorted(records_by_taxon.items())
    ]


def _propagate_alignment(
    node: TreeNode,
    *,
    parent_sequence: list[str],
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
    rng: random.Random,
    records_by_taxon: dict[str, str],
) -> None:
    if node.is_leaf():
        if node.name is None:
            raise ValueError("accuracy fixture leaf is missing a taxon name")
        records_by_taxon[node.name] = "".join(parent_sequence)
        return
    for child in node.children:
        branch_length = float(child.branch_length or 0.0)
        transition_matrix = _transition_matrix_for_case(
            case,
            branch_length=branch_length,
        )
        child_sequence = [
            _sample_state(
                transition_matrix[_state_index(parent_state)].tolist(),
                rng,
            )
            for parent_state in parent_sequence
        ]
        _propagate_alignment(
            child,
            parent_sequence=child_sequence,
            case=case,
            rng=rng,
            records_by_taxon=records_by_taxon,
        )


def _transition_matrix_for_case(
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
    *,
    branch_length: float,
) -> numpy.ndarray:
    if case.model_name == "jc69":
        return jc69_transition_probability_matrix(branch_length)
    if case.model_name == "k80":
        return k80_transition_probability_matrix(
            branch_length,
            kappa=case.parameter_truth["kappa"],
        )
    if case.model_name == "f81":
        return f81_transition_probability_matrix(
            branch_length,
            base_frequencies=_base_frequencies(case),
        )
    return hky85_transition_probability_matrix(
        branch_length,
        base_frequencies=_base_frequencies(case),
        kappa=case.parameter_truth["kappa"],
    )


def _root_state_probabilities(
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
) -> list[float]:
    base_frequencies = _base_frequencies(case)
    return [base_frequencies[state] for state in DNA_STATE_ORDER]


def _base_frequencies(
    case: _GovernedNativeMaximumLikelihoodAccuracyCase,
) -> dict[str, float]:
    if case.model_name in {"jc69", "k80"}:
        return dict(_UNIFORM_DNA_FREQUENCIES)
    return {
        "A": case.parameter_truth["base_frequency_a"],
        "C": case.parameter_truth["base_frequency_c"],
        "G": case.parameter_truth["base_frequency_g"],
        "T": case.parameter_truth["base_frequency_t"],
    }


def _state_index(state: str) -> int:
    return DNA_STATE_ORDER.index(state)


def _sample_state(probabilities: list[float], rng: random.Random) -> str:
    draw = rng.random()
    cumulative = 0.0
    for state, probability in zip(DNA_STATE_ORDER, probabilities, strict=True):
        cumulative += float(probability)
        if draw <= cumulative:
            return state
    return DNA_STATE_ORDER[-1]
