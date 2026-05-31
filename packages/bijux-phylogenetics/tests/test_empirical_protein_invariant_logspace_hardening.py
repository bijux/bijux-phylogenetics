from __future__ import annotations

from decimal import Decimal, getcontext
import math

import numpy

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    build_discrete_gamma_rate_categories,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_empirical_protein_tree_likelihood_with_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    evaluate_fixed_topology_protein_site_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def test_invariant_mixture_empirical_protein_deep_tree_avoids_variable_underflow_failure() -> (
    None
):
    tree = _deep_tree()
    records = _deep_tree_records()
    rate_matrix = _uniform_empirical_rate_matrix()
    root_prior = _uniform_root_prior()

    report = evaluate_empirical_protein_tree_likelihood_with_invariant_mixture(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        invariant_proportion=0.35,
        matrix_label="uniform-deep-tree",
    )

    fixed_rate_log_likelihood = _site_log_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="empirical protein matrix +I",
    )
    expected = math.log(1.0 - 0.35) + fixed_rate_log_likelihood

    assert math.isfinite(report.log_likelihood)
    assert math.isclose(report.log_likelihood, expected, rel_tol=0.0, abs_tol=1e-12)
    assert report.site_likelihoods[0].invariant_component_likelihood == 0.0
    assert report.site_likelihoods[0].mixture_likelihood == 0.0


def test_discrete_gamma_invariant_empirical_protein_deep_tree_matches_high_precision_mixture() -> (
    None
):
    tree = _deep_tree()
    records = _deep_tree_records()
    rate_matrix = _uniform_empirical_rate_matrix()
    root_prior = _uniform_root_prior()
    categories = build_discrete_gamma_rate_categories(alpha=0.8, category_count=4)

    report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            matrix_label="uniform-deep-tree",
        )
    )

    category_log_likelihoods = [
        _site_log_likelihood(
            tree,
            records,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            model_name="empirical protein matrix +G+I",
            rate_scale=category.rate,
        )
        for category in categories
    ]
    variable_component_log_likelihood = _high_precision_weighted_logsumexp(
        category_log_likelihoods,
        [category.weight for category in categories],
    )
    expected = math.log(1.0 - 0.35) + variable_component_log_likelihood

    assert math.isfinite(report.log_likelihood)
    assert math.isclose(report.log_likelihood, expected, rel_tol=0.0, abs_tol=1e-12)
    assert report.site_likelihoods[0].invariant_component_likelihood == 0.0
    assert report.site_likelihoods[0].mixture_likelihood == 0.0


def _deep_tree():
    tip_names = [f"T{index}" for index in range(256)]
    return loads_newick(_balanced_newick(tip_names) + ";")


def _balanced_newick(names: list[str]) -> str:
    if len(names) == 1:
        return names[0]
    midpoint = len(names) // 2
    return (
        f"({_balanced_newick(names[:midpoint])}:5.0,"
        f"{_balanced_newick(names[midpoint:])}:5.0)"
    )


def _deep_tree_records() -> list[AlignmentRecord]:
    state_order = _protein_state_order()
    return [
        AlignmentRecord(
            identifier=f"T{index}",
            sequence=state_order[index % len(state_order)],
        )
        for index in range(256)
    ]


def _uniform_empirical_rate_matrix() -> numpy.ndarray:
    rate_matrix = numpy.full((20, 20), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    return rate_matrix


def _uniform_root_prior() -> numpy.ndarray:
    return numpy.full(20, 1.0 / 20.0, dtype=float)


def _site_log_likelihood(
    tree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    model_name: str,
    rate_scale: float = 1.0,
) -> float:
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            rate_matrix,
            max(float(child.branch_length or 0.0), 0.0) * rate_scale,
        )
        for _parent, child in tree.iter_edges()
    }
    return evaluate_fixed_topology_protein_site_log_likelihood(
        tree,
        tuple(record.sequence for record in records),
        taxon_order=[record.identifier for record in records],
        model_name=model_name,
        root_prior=root_prior,
        transition_matrix_for_child=lambda child: transition_by_node_id[
            child.node_id or ""
        ],
    )


def _high_precision_weighted_logsumexp(
    log_values: list[float],
    weights: list[float],
) -> float:
    getcontext().prec = 90
    maximum = max(Decimal(str(value)) for value in log_values)
    total = Decimal("0")
    for log_value, weight in zip(log_values, weights, strict=True):
        total += Decimal(str(weight)) * (Decimal(str(log_value)) - maximum).exp()
    return float(maximum + total.ln())


def _protein_state_order() -> tuple[str, ...]:
    return (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )
