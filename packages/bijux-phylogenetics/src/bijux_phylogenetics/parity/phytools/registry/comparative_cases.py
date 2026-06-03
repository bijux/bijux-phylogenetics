from __future__ import annotations

from .fixtures import PhytoolsRegistryFixtureCatalog
from .models import PhytoolsParityCase


def build_comparative_cases(
    fixture_catalog: PhytoolsRegistryFixtureCatalog,
) -> list[PhytoolsParityCase]:
    """Build the governed live `phytools` comparative parity cases."""

    pgls_continuous_fixture = fixture_catalog.pgls_continuous_fixture
    pgls_categorical_fixture = fixture_catalog.pgls_categorical_fixture
    pgls_interaction_fixture = fixture_catalog.pgls_interaction_fixture
    phyl_resid_brownian_fixture = fixture_catalog.phyl_resid_brownian_fixture
    phyl_resid_lambda_fixture = fixture_catalog.phyl_resid_lambda_fixture
    phyl_resid_lambda_missing_fixture = (
        fixture_catalog.phyl_resid_lambda_missing_fixture
    )
    phyl_anova_fixture = fixture_catalog.phyl_anova_fixture
    phyl_anova_missing_fixture = fixture_catalog.phyl_anova_missing_fixture
    return [
        PhytoolsParityCase(
            case_id="pgls-sey-brownian-continuous-four-taxa",
            fixture_id=pgls_continuous_fixture.fixture_id,
            function_name="phytools::pgls.SEy",
            python_function_name="run_pgls",
            operation="comparative-pgls-brownian",
            input_fixtures=(
                pgls_continuous_fixture.tree_path,
                pgls_continuous_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=pgls_continuous_fixture.trait_name,
            taxon_column=pgls_continuous_fixture.taxon_column,
            comparative_formula="response ~ predictor_one",
            comparative_lambda_value=1.0,
            row_field_tolerances={"value": 1e-6},
        ),
        PhytoolsParityCase(
            case_id="pgls-sey-brownian-categorical-eight-taxa",
            fixture_id=pgls_categorical_fixture.fixture_id,
            function_name="phytools::pgls.SEy",
            python_function_name="run_pgls",
            operation="comparative-pgls-brownian",
            input_fixtures=(
                pgls_categorical_fixture.tree_path,
                pgls_categorical_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=pgls_categorical_fixture.trait_name,
            taxon_column=pgls_categorical_fixture.taxon_column,
            comparative_formula="response ~ habitat + diet",
            comparative_lambda_value=1.0,
            row_field_tolerances={"value": 1e-6},
        ),
        PhytoolsParityCase(
            case_id="pgls-sey-brownian-interaction-eight-taxa",
            fixture_id=pgls_interaction_fixture.fixture_id,
            function_name="phytools::pgls.SEy",
            python_function_name="run_pgls",
            operation="comparative-pgls-brownian",
            input_fixtures=(
                pgls_interaction_fixture.tree_path,
                pgls_interaction_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=pgls_interaction_fixture.trait_name,
            taxon_column=pgls_interaction_fixture.taxon_column,
            comparative_formula="response ~ habitat * diet",
            comparative_lambda_value=1.0,
            row_field_tolerances={"value": 1e-6},
        ),
        PhytoolsParityCase(
            case_id="phyl-resid-bm-allometry-six-taxa",
            fixture_id=phyl_resid_brownian_fixture.fixture_id,
            function_name="phytools::phyl.resid(method='BM')",
            python_function_name="summarize_phylogenetic_residuals",
            operation="phylogenetic-residuals",
            input_fixtures=(
                phyl_resid_brownian_fixture.tree_path,
                phyl_resid_brownian_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=phyl_resid_brownian_fixture.trait_name,
            taxon_column=phyl_resid_brownian_fixture.taxon_column,
            comparative_predictors=("body_mass",),
            comparative_lambda_value=1.0,
            row_field_tolerances={"value": 1e-6},
        ),
        PhytoolsParityCase(
            case_id="phyl-resid-lambda-allometry-six-taxa",
            fixture_id=phyl_resid_lambda_fixture.fixture_id,
            function_name="phytools::phyl.resid(method='lambda')",
            python_function_name="summarize_phylogenetic_residuals",
            operation="phylogenetic-residuals",
            input_fixtures=(
                phyl_resid_lambda_fixture.tree_path,
                phyl_resid_lambda_fixture.traits_path,
            ),
            tolerance=5e-4,
            trait_name=phyl_resid_lambda_fixture.trait_name,
            taxon_column=phyl_resid_lambda_fixture.taxon_column,
            comparative_predictors=("body_mass",),
            field_tolerances={"lambda_value": 5e-4, "log_likelihood": 5e-4},
            row_field_tolerances={"value": 5e-4},
        ),
        PhytoolsParityCase(
            case_id="phyl-resid-lambda-missing-six-taxa",
            fixture_id=phyl_resid_lambda_missing_fixture.fixture_id,
            function_name="phytools::phyl.resid(method='lambda')",
            python_function_name="summarize_phylogenetic_residuals",
            operation="phylogenetic-residuals",
            input_fixtures=(
                phyl_resid_lambda_missing_fixture.tree_path,
                phyl_resid_lambda_missing_fixture.traits_path,
            ),
            tolerance=5e-4,
            trait_name=phyl_resid_lambda_missing_fixture.trait_name,
            taxon_column=phyl_resid_lambda_missing_fixture.taxon_column,
            comparative_predictors=("body_mass",),
            field_tolerances={"lambda_value": 5e-4, "log_likelihood": 5e-4},
            row_field_tolerances={"value": 5e-4},
        ),
        PhytoolsParityCase(
            case_id="phyl-anova-group-effect-six-taxa",
            fixture_id=phyl_anova_fixture.fixture_id,
            function_name="phytools::phylANOVA",
            python_function_name="summarize_phylogenetic_anova",
            operation="phylogenetic-anova",
            input_fixtures=(
                phyl_anova_fixture.tree_path,
                phyl_anova_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=phyl_anova_fixture.trait_name,
            taxon_column=phyl_anova_fixture.taxon_column,
            comparative_predictors=("habitat",),
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={"p_value": 0.15},
            row_field_tolerances={
                "uncorrected_p_value": 0.15,
                "adjusted_p_value": 0.15,
            },
        ),
        PhytoolsParityCase(
            case_id="phyl-anova-group-effect-missing-six-taxa",
            fixture_id=phyl_anova_missing_fixture.fixture_id,
            function_name="phytools::phylANOVA",
            python_function_name="summarize_phylogenetic_anova",
            operation="phylogenetic-anova",
            input_fixtures=(
                phyl_anova_missing_fixture.tree_path,
                phyl_anova_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=phyl_anova_missing_fixture.trait_name,
            taxon_column=phyl_anova_missing_fixture.taxon_column,
            comparative_predictors=("habitat",),
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={"p_value": 0.15},
            row_field_tolerances={
                "uncorrected_p_value": 0.15,
                "adjusted_p_value": 0.15,
            },
        ),
    ]
