from __future__ import annotations

from bijux_phylogenetics.simulation import DiscreteHistoryRateRow

from .fixtures import PhytoolsRegistryFixtureCatalog
from .models import PhytoolsParityCase


def build_discrete_history_cases(
    fixture_catalog: PhytoolsRegistryFixtureCatalog,
) -> list[PhytoolsParityCase]:
    """Build the governed live `phytools` discrete-history parity cases."""

    simulation_tree_fixture = fixture_catalog.simulation_tree_fixture
    simulation_six_taxa_tree_fixture = fixture_catalog.simulation_six_taxa_tree_fixture
    binary_discrete_fixture = fixture_catalog.binary_discrete_fixture
    multistate_discrete_fixture = fixture_catalog.multistate_discrete_fixture
    binary_discrete_missing_fixture = fixture_catalog.binary_discrete_missing_fixture
    multistate_discrete_missing_fixture = (
        fixture_catalog.multistate_discrete_missing_fixture
    )
    return [
        PhytoolsParityCase(
            case_id="sim-history-binary-no-change-example-tree",
            fixture_id="example_tree_discrete_history_binary_no_change",
            function_name="phytools::sim.history",
            python_function_name="simulate_discrete_histories",
            operation="simulate-discrete-history",
            input_fixtures=(simulation_tree_fixture,),
            tolerance=1e-6,
            trait_name="simulated_state",
            simulation_states=("0", "1"),
            simulation_rate_rows=(
                DiscreteHistoryRateRow("0", "1", 0.01),
                DiscreteHistoryRateRow("1", "0", 0.01),
            ),
            simulation_root_state="0",
            simulation_replicate_count=128,
            simulation_seed=17,
            field_tolerances={
                "mean_total_transition_count": 0.5,
                "lower_95_total_transition_count": 1.0,
                "upper_95_total_transition_count": 1.0,
            },
            row_field_tolerances={
                "mean_value": 0.25,
                "lower_95_interval": 0.35,
                "upper_95_interval": 0.35,
                "presence_fraction": 0.2,
            },
        ),
        PhytoolsParityCase(
            case_id="sim-history-binary-high-rate-example-tree",
            fixture_id="example_tree_discrete_history_binary_high_rate",
            function_name="phytools::sim.history",
            python_function_name="simulate_discrete_histories",
            operation="simulate-discrete-history",
            input_fixtures=(simulation_tree_fixture,),
            tolerance=1e-6,
            trait_name="simulated_state",
            simulation_states=("0", "1"),
            simulation_rate_rows=(
                DiscreteHistoryRateRow("0", "1", 4.0),
                DiscreteHistoryRateRow("1", "0", 2.5),
            ),
            simulation_root_state="0",
            simulation_replicate_count=128,
            simulation_seed=17,
            field_tolerances={
                "mean_total_transition_count": 2.0,
                "lower_95_total_transition_count": 4.0,
                "upper_95_total_transition_count": 4.0,
            },
            row_field_tolerances={
                "mean_value": 1.0,
                "lower_95_interval": 2.0,
                "upper_95_interval": 2.0,
                "presence_fraction": 0.2,
            },
        ),
        PhytoolsParityCase(
            case_id="sim-history-multistate-no-change-six-taxa",
            fixture_id="example_tree_six_taxa_discrete_history_multistate_no_change",
            function_name="phytools::sim.history",
            python_function_name="simulate_discrete_histories",
            operation="simulate-discrete-history",
            input_fixtures=(simulation_six_taxa_tree_fixture,),
            tolerance=1e-6,
            trait_name="simulated_state",
            simulation_states=("red", "blue", "green"),
            simulation_rate_rows=(
                DiscreteHistoryRateRow("red", "blue", 0.01),
                DiscreteHistoryRateRow("red", "green", 0.005),
                DiscreteHistoryRateRow("blue", "red", 0.01),
                DiscreteHistoryRateRow("blue", "green", 0.01),
                DiscreteHistoryRateRow("green", "red", 0.005),
                DiscreteHistoryRateRow("green", "blue", 0.01),
            ),
            simulation_root_state="red",
            simulation_replicate_count=128,
            simulation_seed=17,
            field_tolerances={
                "mean_total_transition_count": 0.75,
                "lower_95_total_transition_count": 1.5,
                "upper_95_total_transition_count": 1.5,
            },
            row_field_tolerances={
                "mean_value": 0.15,
                "lower_95_interval": 1.25,
                "upper_95_interval": 1.6,
                "presence_fraction": 0.03,
            },
        ),
        PhytoolsParityCase(
            case_id="sim-history-multistate-high-rate-six-taxa",
            fixture_id="example_tree_six_taxa_discrete_history_multistate_high_rate",
            function_name="phytools::sim.history",
            python_function_name="simulate_discrete_histories",
            operation="simulate-discrete-history",
            input_fixtures=(simulation_six_taxa_tree_fixture,),
            tolerance=1e-6,
            trait_name="simulated_state",
            simulation_states=("red", "blue", "green"),
            simulation_rate_rows=(
                DiscreteHistoryRateRow("red", "blue", 6.0),
                DiscreteHistoryRateRow("red", "green", 3.5),
                DiscreteHistoryRateRow("blue", "red", 4.5),
                DiscreteHistoryRateRow("blue", "green", 5.5),
                DiscreteHistoryRateRow("green", "red", 3.0),
                DiscreteHistoryRateRow("green", "blue", 4.0),
            ),
            simulation_root_state="red",
            simulation_replicate_count=128,
            simulation_seed=17,
            field_tolerances={
                "mean_total_transition_count": 3.0,
                "lower_95_total_transition_count": 6.0,
                "upper_95_total_transition_count": 6.0,
            },
            row_field_tolerances={
                "mean_value": 1.5,
                "lower_95_interval": 3.0,
                "upper_95_interval": 3.0,
                "presence_fraction": 0.2,
            },
        ),
        PhytoolsParityCase(
            case_id="sim-history-binary-root-prior-example-tree",
            fixture_id="example_tree_discrete_history_binary_root_prior",
            function_name="phytools::sim.history",
            python_function_name="simulate_discrete_histories",
            operation="simulate-discrete-history",
            input_fixtures=(simulation_tree_fixture,),
            tolerance=1e-6,
            trait_name="simulated_state",
            simulation_states=("0", "1"),
            simulation_rate_rows=(
                DiscreteHistoryRateRow("0", "1", 0.5),
                DiscreteHistoryRateRow("1", "0", 1.0),
            ),
            simulation_root_state_probabilities={"0": 0.2, "1": 0.8},
            simulation_replicate_count=128,
            simulation_seed=17,
            field_tolerances={
                "mean_total_transition_count": 1.0,
                "lower_95_total_transition_count": 2.0,
                "upper_95_total_transition_count": 2.0,
            },
            row_field_tolerances={
                "mean_value": 0.35,
                "lower_95_interval": 0.5,
                "upper_95_interval": 0.5,
                "presence_fraction": 0.2,
            },
        ),
        PhytoolsParityCase(
            case_id="rerooting-er-binary-twenty-four-taxa",
            fixture_id=binary_discrete_fixture.fixture_id,
            function_name="phytools::rerootingMethod",
            python_function_name="reconstruct_discrete_ancestral_states",
            operation="discrete-ancestral-rerooting",
            input_fixtures=(
                binary_discrete_fixture.tree_path,
                binary_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_fixture.trait_name,
            taxon_column=binary_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"probability": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="rerooting-er-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::rerootingMethod",
            python_function_name="reconstruct_discrete_ancestral_states",
            operation="discrete-ancestral-rerooting",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"probability": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="rerooting-er-binary-missing-twenty-four-taxa",
            fixture_id=binary_discrete_missing_fixture.fixture_id,
            function_name="phytools::rerootingMethod",
            python_function_name="reconstruct_discrete_ancestral_states",
            operation="discrete-ancestral-rerooting",
            input_fixtures=(
                binary_discrete_missing_fixture.tree_path,
                binary_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_missing_fixture.trait_name,
            taxon_column=binary_discrete_missing_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"probability": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="rerooting-sym-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::rerootingMethod",
            python_function_name="reconstruct_discrete_ancestral_states",
            operation="discrete-ancestral-rerooting",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="symmetric",
            row_field_tolerances={"probability": 5e-5},
        ),
        PhytoolsParityCase(
            case_id="rerooting-sym-multistate-missing-twenty-four-taxa",
            fixture_id=multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::rerootingMethod",
            python_function_name="reconstruct_discrete_ancestral_states",
            operation="discrete-ancestral-rerooting",
            input_fixtures=(
                multistate_discrete_missing_fixture.tree_path,
                multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_missing_fixture.trait_name,
            taxon_column=multistate_discrete_missing_fixture.taxon_column,
            discrete_model="symmetric",
            row_field_tolerances={"probability": 5e-5},
        ),
    ]
