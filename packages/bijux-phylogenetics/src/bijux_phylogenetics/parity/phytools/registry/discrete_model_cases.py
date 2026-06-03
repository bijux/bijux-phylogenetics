from __future__ import annotations

from .fixtures import PhytoolsRegistryFixtureCatalog
from .models import PhytoolsParityCase


def build_discrete_model_cases(
    fixture_catalog: PhytoolsRegistryFixtureCatalog,
) -> list[PhytoolsParityCase]:
    """Build the governed live `phytools::fitMk` parity cases."""

    binary_discrete_fixture = fixture_catalog.binary_discrete_fixture
    multistate_discrete_fixture = fixture_catalog.multistate_discrete_fixture
    binary_discrete_missing_fixture = fixture_catalog.binary_discrete_missing_fixture
    multistate_discrete_missing_fixture = (
        fixture_catalog.multistate_discrete_missing_fixture
    )
    ard_binary_discrete_fixture = fixture_catalog.ard_binary_discrete_fixture
    ard_multistate_discrete_fixture = fixture_catalog.ard_multistate_discrete_fixture
    ard_binary_discrete_missing_fixture = (
        fixture_catalog.ard_binary_discrete_missing_fixture
    )
    ard_multistate_discrete_missing_fixture = (
        fixture_catalog.ard_multistate_discrete_missing_fixture
    )
    return [
        PhytoolsParityCase(
            case_id="fitmk-er-binary-twenty-four-taxa",
            fixture_id=binary_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                binary_discrete_fixture.tree_path,
                binary_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_fixture.trait_name,
            taxon_column=binary_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-binary-missing-twenty-four-taxa",
            fixture_id=binary_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                binary_discrete_missing_fixture.tree_path,
                binary_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_missing_fixture.trait_name,
            taxon_column=binary_discrete_missing_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-multistate-missing-twenty-four-taxa",
            fixture_id=multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_missing_fixture.tree_path,
                multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_missing_fixture.trait_name,
            taxon_column=multistate_discrete_missing_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-sym-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='SYM')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="symmetric",
            field_tolerances={
                "log_likelihood": 2e-4,
                "aic": 2e-4,
                "aicc": 2e-4,
            },
            row_field_tolerances={"rate": 1e-4},
        ),
        PhytoolsParityCase(
            case_id="fitmk-sym-multistate-missing-twenty-four-taxa",
            fixture_id=multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='SYM')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_missing_fixture.tree_path,
                multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_missing_fixture.trait_name,
            taxon_column=multistate_discrete_missing_fixture.taxon_column,
            discrete_model="symmetric",
            field_tolerances={
                "log_likelihood": 2e-4,
                "aic": 2e-4,
                "aicc": 2e-4,
            },
            row_field_tolerances={"rate": 1e-4},
        ),
        PhytoolsParityCase(
            case_id="fitmk-ard-binary-twenty-four-taxa",
            fixture_id=ard_binary_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ARD')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                ard_binary_discrete_fixture.tree_path,
                ard_binary_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=ard_binary_discrete_fixture.trait_name,
            taxon_column=ard_binary_discrete_fixture.taxon_column,
            discrete_model="all-rates-different",
            field_tolerances={
                "log_likelihood": 1e-3,
                "aic": 1e-3,
                "aicc": 1e-3,
            },
            row_field_tolerances={"rate": 1e-3},
        ),
        PhytoolsParityCase(
            case_id="fitmk-ard-multistate-twenty-four-taxa",
            fixture_id=ard_multistate_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ARD')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                ard_multistate_discrete_fixture.tree_path,
                ard_multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=ard_multistate_discrete_fixture.trait_name,
            taxon_column=ard_multistate_discrete_fixture.taxon_column,
            discrete_model="all-rates-different",
            field_tolerances={
                "log_likelihood": 1e-3,
                "aic": 1e-3,
                "aicc": 1e-3,
            },
            row_field_tolerances={"rate": 1e-3},
            compare_rows=False,
        ),
        PhytoolsParityCase(
            case_id="fitmk-ard-binary-missing-twenty-four-taxa",
            fixture_id=ard_binary_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ARD')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                ard_binary_discrete_missing_fixture.tree_path,
                ard_binary_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=ard_binary_discrete_missing_fixture.trait_name,
            taxon_column=ard_binary_discrete_missing_fixture.taxon_column,
            discrete_model="all-rates-different",
            field_tolerances={
                "log_likelihood": 1e-3,
                "aic": 1e-3,
                "aicc": 1e-3,
            },
            row_field_tolerances={"rate": 1e-3},
        ),
        PhytoolsParityCase(
            case_id="fitmk-ard-multistate-missing-twenty-four-taxa",
            fixture_id=ard_multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ARD')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                ard_multistate_discrete_missing_fixture.tree_path,
                ard_multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=ard_multistate_discrete_missing_fixture.trait_name,
            taxon_column=ard_multistate_discrete_missing_fixture.taxon_column,
            discrete_model="all-rates-different",
            field_tolerances={
                "log_likelihood": 1e-3,
                "aic": 1e-3,
                "aicc": 1e-3,
            },
            row_field_tolerances={"rate": 1e-3},
            compare_rows=False,
        ),
    ]
