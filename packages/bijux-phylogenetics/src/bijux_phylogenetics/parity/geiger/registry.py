from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.fixtures import get_shared_geiger_continuous_fixture


@dataclass(frozen=True, slots=True)
class GeigerParityCase:
    """One governed live `geiger` parity case."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    model_name: str
    python_mode: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    trait_name: str
    taxon_column: str | None = None
    optimizer_settings: dict[str, object] | None = None
    comparison_fields: tuple[str, ...] = ()
    ou_bounds: tuple[float, float] | None = None
    early_burst_bounds: tuple[float, float] | None = None
    field_tolerances: dict[str, float] | None = None


def _package_root() -> Path:
    return Path(__file__).resolve().parents[4]


def list_geiger_parity_cases() -> list[GeigerParityCase]:
    """Return the governed live `geiger` parity cases."""
    package_root = _package_root()
    tests_root = package_root / "tests" / "fixtures"
    brownian_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )
    brownian_missing_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_missing_values_twenty_four_taxa"
    )
    ou_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_ou_known_truth_twenty_four_taxa"
    )
    early_burst_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
    )
    return [
        GeigerParityCase(
            case_id="fitcontinuous-bm-example-tree",
            fixture_id="example_tree_comparative_response",
            function_name="geiger::fitContinuous(model='BM')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="BM",
            python_mode="brownian",
            input_fixtures=(
                tests_root / "trees" / "example_tree.nwk",
                tests_root / "metadata" / "example_traits_comparative.tsv",
            ),
            tolerance=1e-6,
            trait_name="response",
            taxon_column="taxon",
            optimizer_settings={
                "reference_control_policy": "fitcontinuous-default",
                "bijux_optimizer_name": "closed-form-profile-solution",
                "bijux_parameter_search": "none",
            },
            comparison_fields=(
                "taxon_count",
                "trait_name",
                "model_name",
                "excluded_taxon_count",
                "missing_value_policy",
                "standard_error_policy",
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
                "aicc",
            ),
            field_tolerances={"aicc": 1e-5},
        ),
        GeigerParityCase(
            case_id="fitcontinuous-bm-brownian-sigma-recovery",
            fixture_id=brownian_fixture.fixture_id,
            function_name="geiger::fitContinuous(model='BM')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="BM",
            python_mode="brownian",
            input_fixtures=(brownian_fixture.tree_path, brownian_fixture.traits_path),
            tolerance=0.2,
            trait_name=brownian_fixture.trait_name,
            taxon_column=brownian_fixture.taxon_column,
            optimizer_settings={
                "reference_control_policy": "fitcontinuous-default",
                "bijux_optimizer_name": "closed-form-profile-solution",
                "bijux_parameter_search": "none",
            },
            comparison_fields=(
                "taxon_count",
                "trait_name",
                "model_name",
                "excluded_taxon_count",
                "missing_value_policy",
                "standard_error_policy",
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
                "aicc",
            ),
            field_tolerances={"aicc": 1e-5},
        ),
        GeigerParityCase(
            case_id="fitcontinuous-bm-missing-values-review",
            fixture_id=brownian_missing_fixture.fixture_id,
            function_name="geiger::fitContinuous(model='BM')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="BM",
            python_mode="brownian",
            input_fixtures=(
                brownian_missing_fixture.tree_path,
                brownian_missing_fixture.traits_path,
            ),
            tolerance=0.2,
            trait_name=brownian_missing_fixture.trait_name,
            taxon_column=brownian_missing_fixture.taxon_column,
            optimizer_settings={
                "reference_control_policy": "fitcontinuous-default",
                "bijux_optimizer_name": "closed-form-profile-solution",
                "bijux_parameter_search": "none",
            },
            comparison_fields=(
                "taxon_count",
                "trait_name",
                "model_name",
                "excluded_taxon_count",
                "excluded_taxa",
                "missing_value_taxa",
                "non_numeric_taxa",
                "missing_from_traits",
                "missing_value_policy",
                "standard_error_policy",
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
                "aicc",
            ),
            field_tolerances={"aicc": 1e-5},
        ),
        GeigerParityCase(
            case_id="fitcontinuous-ou-ou-parameter-recovery",
            fixture_id=ou_fixture.fixture_id,
            function_name="geiger::fitContinuous(model='OU')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="OU",
            python_mode="ornstein-uhlenbeck",
            input_fixtures=(ou_fixture.tree_path, ou_fixture.traits_path),
            tolerance=0.75,
            trait_name=ou_fixture.trait_name,
            taxon_column=ou_fixture.taxon_column,
            optimizer_settings={
                "reference_control_policy": "fitcontinuous-default",
                "bijux_optimizer_name": "governed-two-stage-grid-search",
                "bijux_coarse_grid_point_count": 81,
                "bijux_fine_grid_point_count": 81,
                "bijux_parameter_bounds": {"lower": 0.0, "upper": 10.0},
            },
            comparison_fields=(
                "taxon_count",
                "trait_name",
                "model_name",
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
                "parameter_name",
                "parameter_value",
            ),
            ou_bounds=(0.0, 10.0),
        ),
        GeigerParityCase(
            case_id="fitcontinuous-eb-early-burst-rate-recovery",
            fixture_id=early_burst_fixture.fixture_id,
            function_name="geiger::fitContinuous(model='EB')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="EB",
            python_mode="early-burst",
            input_fixtures=(
                early_burst_fixture.tree_path,
                early_burst_fixture.traits_path,
            ),
            tolerance=0.5,
            trait_name=early_burst_fixture.trait_name,
            taxon_column=early_burst_fixture.taxon_column,
            optimizer_settings={
                "reference_control_policy": "fitcontinuous-default",
                "bijux_optimizer_name": "governed-two-stage-grid-search",
                "bijux_coarse_grid_point_count": 81,
                "bijux_fine_grid_point_count": 81,
                "bijux_parameter_bounds": {"lower": 0.0, "upper": 10.0},
            },
            comparison_fields=(
                "taxon_count",
                "trait_name",
                "model_name",
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
                "parameter_name",
                "parameter_value",
            ),
            early_burst_bounds=(0.0, 10.0),
        ),
    ]
