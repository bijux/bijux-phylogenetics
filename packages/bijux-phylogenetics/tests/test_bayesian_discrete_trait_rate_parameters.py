from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.discrete_trait_rate_parameters import (
    parameterize_discrete_trait_rate_rows,
    resolve_discrete_trait_rate_rows,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_phytools_comparative_fixture,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_equal_rates_parameterization_collapses_to_one_shared_rate() -> None:
    fit_report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    parameterization = parameterize_discrete_trait_rate_rows(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
    )

    assert parameterization.model == "equal-rates"
    assert parameterization.parameter_count == 1
    assert parameterization.parameter_values.keys() == {"shared-rate"}
    assert len(parameterization.groups[0].transition_pairs) == 6


@pytest.mark.slow
def test_symmetric_parameterization_groups_bidirectional_pairs_once() -> None:
    fit_report = fit_discrete_mk_model(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="region_state",
        taxon_column="taxon",
        model="symmetric",
    )

    parameterization = parameterize_discrete_trait_rate_rows(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
    )

    assert parameterization.parameter_count == 3
    assert set(parameterization.parameter_values) == {
        "north<->south",
        "north<->west",
        "south<->west",
    }
    assert all(len(group.transition_pairs) == 2 for group in parameterization.groups)


def test_ard_parameter_resolution_updates_directional_rows() -> None:
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )
    fit_report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
    )
    parameterization = parameterize_discrete_trait_rate_rows(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
    )

    resolved_rows = resolve_discrete_trait_rate_rows(
        model=fit_report.model,
        transition_rate_rows=fit_report.transition_rate_rows,
        parameter_values={"0->1": 0.25, "1->0": 0.75},
    )
    rate_lookup = {
        (row.source_state, row.target_state): row.rate
        for row in resolved_rows
        if row.transition_allowed
    }

    assert parameterization.parameter_count == 2
    assert set(parameterization.parameter_values) == {"0->1", "1->0"}
    assert math.isclose(rate_lookup[("0", "1")], 0.25, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(rate_lookup[("1", "0")], 0.75, rel_tol=0.0, abs_tol=1e-12)
