from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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


def _continuous_mode_recovery_root() -> Path:
    return (
        _package_root()
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "datasets"
        / "simulation"
        / "continuous_mode_recovery_panel"
    )


def list_geiger_parity_cases() -> list[GeigerParityCase]:
    """Return the governed live `geiger` parity cases."""
    package_root = _package_root()
    tests_root = package_root / "tests" / "fixtures"
    recovery_root = _continuous_mode_recovery_root()
    recovery_tree_path = recovery_root / "reference-tree.nwk"
    recovery_traits_root = recovery_root / "expected" / "simulated-traits"
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
            fixture_id="continuous_mode_recovery_panel:brownian-sigma-recovery",
            function_name="geiger::fitContinuous(model='BM')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="BM",
            python_mode="brownian",
            input_fixtures=(
                recovery_tree_path,
                recovery_traits_root / "brownian-sigma-recovery.tsv",
            ),
            tolerance=0.2,
            trait_name="value",
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
                "root_state",
                "rate",
                "log_likelihood",
                "aic",
            ),
        ),
        GeigerParityCase(
            case_id="fitcontinuous-ou-ou-parameter-recovery",
            fixture_id="continuous_mode_recovery_panel:ou-parameter-recovery",
            function_name="geiger::fitContinuous(model='OU')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="OU",
            python_mode="ornstein-uhlenbeck",
            input_fixtures=(
                recovery_tree_path,
                recovery_traits_root / "ou-parameter-recovery.tsv",
            ),
            tolerance=0.75,
            trait_name="value",
            taxon_column="taxon",
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
            fixture_id="continuous_mode_recovery_panel:early-burst-rate-recovery",
            function_name="geiger::fitContinuous(model='EB')",
            python_function_name="fit_continuous_evolutionary_mode",
            operation="fit-continuous",
            model_name="EB",
            python_mode="early-burst",
            input_fixtures=(
                recovery_tree_path,
                recovery_traits_root / "early-burst-rate-recovery.tsv",
            ),
            tolerance=0.5,
            trait_name="value",
            taxon_column="taxon",
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
