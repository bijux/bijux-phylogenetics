from __future__ import annotations

from .fixtures import PhytoolsRegistryFixtureCatalog
from .models import PhytoolsParityCase


def build_signal_cases(
    fixture_catalog: PhytoolsRegistryFixtureCatalog,
) -> list[PhytoolsParityCase]:
    """Build the governed live `phytools::phylosig` parity cases."""

    nonultrametric_signal_fixture = fixture_catalog.nonultrametric_signal_fixture
    weak_signal_fixture = fixture_catalog.weak_signal_fixture
    strong_signal_fixture = fixture_catalog.strong_signal_fixture
    return [
        PhytoolsParityCase(
            case_id="phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa",
            fixture_id=nonultrametric_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='lambda')",
            python_function_name="estimate_pagels_lambda",
            operation="phylogenetic-signal-lambda",
            input_fixtures=(
                nonultrametric_signal_fixture.tree_path,
                nonultrametric_signal_fixture.traits_path,
            ),
            tolerance=5e-4,
            trait_name=nonultrametric_signal_fixture.trait_name,
            taxon_column=nonultrametric_signal_fixture.taxon_column,
        ),
        PhytoolsParityCase(
            case_id="phylosig-lambda-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='lambda')",
            python_function_name="estimate_pagels_lambda",
            operation="phylogenetic-signal-lambda",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-3,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
        ),
        PhytoolsParityCase(
            case_id="phylosig-k-strong-signal-twenty-four-taxa",
            fixture_id=strong_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='K')",
            python_function_name="compute_phylogenetic_signal_test",
            operation="phylogenetic-signal-k",
            input_fixtures=(
                strong_signal_fixture.tree_path,
                strong_signal_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=strong_signal_fixture.trait_name,
            taxon_column=strong_signal_fixture.taxon_column,
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={
                "p_value": 0.03,
                "simulated_k_minimum": 0.01,
                "simulated_k_mean": 0.01,
            },
        ),
        PhytoolsParityCase(
            case_id="phylosig-k-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='K')",
            python_function_name="compute_phylogenetic_signal_test",
            operation="phylogenetic-signal-k",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={
                "p_value": 0.03,
                "simulated_k_minimum": 0.01,
                "simulated_k_mean": 0.01,
            },
        ),
    ]
