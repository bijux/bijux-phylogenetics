from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
    count_rooted_labeled_bifurcating_topologies,
    validate_tree_topology_prior_taxa,
)
from bijux_phylogenetics.runtime.errors import DuplicateTaxonError, PhylogeneticsError


@pytest.mark.parametrize(
    ("taxon_count", "expected_count"),
    [
        (2, 1),
        (3, 3),
        (4, 15),
        (5, 105),
    ],
)
def test_count_rooted_labeled_bifurcating_topologies_matches_known_values(
    taxon_count: int,
    expected_count: int,
) -> None:
    assert count_rooted_labeled_bifurcating_topologies(taxon_count) == expected_count


def test_build_uniform_rooted_tree_topology_prior_matches_analytical_log_probability() -> (
    None
):
    prior_model = build_uniform_rooted_tree_topology_prior(["D", "B", "A", "C"])

    assert prior_model.family == "uniform-rooted-labeled-bifurcating"
    assert prior_model.taxa == ["A", "B", "C", "D"]
    assert prior_model.topology_count == 15
    assert math.isclose(
        prior_model.log_topology_probability,
        -math.log(15.0),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_validate_tree_topology_prior_taxa_sorts_distinct_taxa() -> None:
    assert validate_tree_topology_prior_taxa(["Gamma", "Alpha", "Beta"]) == [
        "Alpha",
        "Beta",
        "Gamma",
    ]


def test_validate_tree_topology_prior_taxa_rejects_duplicates() -> None:
    with pytest.raises(
        DuplicateTaxonError,
        match="tree topology prior requires distinct taxa; duplicates: Alpha",
    ):
        validate_tree_topology_prior_taxa(["Alpha", "Beta", "Alpha"])


@pytest.mark.parametrize(
    ("taxa", "message"),
    [
        (["Alpha"], "at least two taxa"),
        (["Alpha", " "], "does not allow blank taxon labels"),
    ],
)
def test_validate_tree_topology_prior_taxa_rejects_invalid_taxon_sets(
    taxa: list[str],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        validate_tree_topology_prior_taxa(taxa)
