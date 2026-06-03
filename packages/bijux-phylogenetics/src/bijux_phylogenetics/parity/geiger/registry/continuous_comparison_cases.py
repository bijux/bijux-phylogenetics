from __future__ import annotations

from .fixtures import GeigerRegistryFixtureCatalog
from .models import FITCONTINUOUS_REFERENCE_MODEL_ORDER, GeigerParityCase

_FITCONTINUOUS_MODEL_COMPARISON_FIELDS = (
    "taxon_count",
    "trait_name",
    "model_name",
    "selected_model",
    "comparable_model_count",
    "noncomparable_model_count",
    "runner_up_model",
    "runner_up_aicc_delta",
    "warning_count",
    "optimizer_settings",
)
_FITCONTINUOUS_MODEL_COMPARISON_OPTIMIZER = {
    "reference_control_policy": "fitcontinuous-ranked-model-comparison",
    "candidate_models": list(FITCONTINUOUS_REFERENCE_MODEL_ORDER),
    "bijux_comparison_surface": "compare_fitcontinuous_model_ranking",
}


def _comparison_case(
    *,
    case_id: str,
    fixture_id: str,
    input_fixtures: tuple[object, object],
    tolerance: float,
    trait_name: str,
    taxon_column: str,
    kappa_bounds: tuple[float, float],
) -> GeigerParityCase:
    return GeigerParityCase(
        case_id=case_id,
        fixture_id=fixture_id,
        function_name="geiger::fitContinuous(model comparison)",
        python_function_name="compare_fitcontinuous_model_ranking",
        operation="compare-fitcontinuous-models",
        model_name="model-comparison",
        python_mode="fitcontinuous-model-ranking",
        input_fixtures=input_fixtures,
        tolerance=tolerance,
        trait_name=trait_name,
        taxon_column=taxon_column,
        optimizer_settings=_FITCONTINUOUS_MODEL_COMPARISON_OPTIMIZER,
        candidate_model_names=FITCONTINUOUS_REFERENCE_MODEL_ORDER,
        lambda_bounds=(0.0, 1.0),
        kappa_bounds=kappa_bounds,
        delta_bounds=(0.0, 3.0),
        ou_bounds=(0.0, 10.0),
        early_burst_bounds=(0.0, 50.0),
        comparison_fields=_FITCONTINUOUS_MODEL_COMPARISON_FIELDS,
        field_tolerances={"runner_up_aicc_delta": 1e-2},
    )


def build_continuous_comparison_cases(
    fixture_catalog: GeigerRegistryFixtureCatalog,
) -> list[GeigerParityCase]:
    """Build the governed live `geiger` model-comparison parity cases."""

    return [
        _comparison_case(
            case_id="fitcontinuous-model-comparison-brownian-review",
            fixture_id=fixture_catalog.comparison_brownian_fixture.fixture_id,
            input_fixtures=(
                fixture_catalog.comparison_brownian_fixture.tree_path,
                fixture_catalog.comparison_brownian_fixture.traits_path,
            ),
            tolerance=0.2,
            trait_name=fixture_catalog.comparison_brownian_fixture.trait_name,
            taxon_column=fixture_catalog.comparison_brownian_fixture.taxon_column,
            kappa_bounds=(0.0, 1.0),
        ),
        _comparison_case(
            case_id="fitcontinuous-model-comparison-ou-review",
            fixture_id=fixture_catalog.comparison_ou_fixture.fixture_id,
            input_fixtures=(
                fixture_catalog.comparison_ou_fixture.tree_path,
                fixture_catalog.comparison_ou_fixture.traits_path,
            ),
            tolerance=0.2,
            trait_name=fixture_catalog.comparison_ou_fixture.trait_name,
            taxon_column=fixture_catalog.comparison_ou_fixture.taxon_column,
            kappa_bounds=(0.0, 1.0),
        ),
        _comparison_case(
            case_id="fitcontinuous-model-comparison-early-burst-review",
            fixture_id=fixture_catalog.comparison_early_burst_fixture.fixture_id,
            input_fixtures=(
                fixture_catalog.comparison_early_burst_fixture.tree_path,
                fixture_catalog.comparison_early_burst_fixture.traits_path,
            ),
            tolerance=0.2,
            trait_name=fixture_catalog.comparison_early_burst_fixture.trait_name,
            taxon_column=fixture_catalog.comparison_early_burst_fixture.taxon_column,
            kappa_bounds=(0.0, 1.0),
        ),
        _comparison_case(
            case_id="fitcontinuous-model-comparison-white-review",
            fixture_id=fixture_catalog.comparison_white_fixture.fixture_id,
            input_fixtures=(
                fixture_catalog.comparison_white_fixture.tree_path,
                fixture_catalog.comparison_white_fixture.traits_path,
            ),
            tolerance=0.2,
            trait_name=fixture_catalog.comparison_white_fixture.trait_name,
            taxon_column=fixture_catalog.comparison_white_fixture.taxon_column,
            kappa_bounds=(0.0, 3.0),
        ),
    ]
