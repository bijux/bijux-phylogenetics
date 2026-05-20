from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures.geiger_continuous import (
    SharedGeigerContinuousFixture,
    get_shared_geiger_continuous_fixture,
)
from bijux_phylogenetics.datasets.shared_fixtures.geiger_discrete import (
    SharedGeigerDiscreteFixture,
    get_shared_geiger_discrete_fixture,
)


@dataclass(frozen=True, slots=True)
class GeigerRegistryFixtureCatalog:
    """Resolved fixture catalog for the governed `geiger` parity registry."""

    tests_root: Path
    brownian_fixture: SharedGeigerContinuousFixture
    brownian_missing_fixture: SharedGeigerContinuousFixture
    white_strong_fixture: SharedGeigerContinuousFixture
    white_weak_fixture: SharedGeigerContinuousFixture
    white_missing_fixture: SharedGeigerContinuousFixture
    lambda_strong_fixture: SharedGeigerContinuousFixture
    lambda_weak_fixture: SharedGeigerContinuousFixture
    lambda_missing_fixture: SharedGeigerContinuousFixture
    kappa_strong_fixture: SharedGeigerContinuousFixture
    kappa_weak_fixture: SharedGeigerContinuousFixture
    kappa_missing_fixture: SharedGeigerContinuousFixture
    delta_strong_fixture: SharedGeigerContinuousFixture
    delta_weak_fixture: SharedGeigerContinuousFixture
    delta_missing_fixture: SharedGeigerContinuousFixture
    ou_missing_fixture: SharedGeigerContinuousFixture
    ou_fixture: SharedGeigerContinuousFixture
    ou_lower_boundary_fixture: SharedGeigerContinuousFixture
    early_burst_fixture: SharedGeigerContinuousFixture
    early_burst_boundary_fixture: SharedGeigerContinuousFixture
    comparison_brownian_fixture: SharedGeigerContinuousFixture
    comparison_ou_fixture: SharedGeigerContinuousFixture
    comparison_early_burst_fixture: SharedGeigerContinuousFixture
    comparison_white_fixture: SharedGeigerContinuousFixture
    discrete_er_binary_fixture: SharedGeigerDiscreteFixture
    discrete_transform_weak_fixture: SharedGeigerDiscreteFixture
    discrete_lambda_missing_fixture: SharedGeigerDiscreteFixture
    discrete_kappa_strong_fixture: SharedGeigerDiscreteFixture
    discrete_kappa_weak_fixture: SharedGeigerDiscreteFixture
    discrete_kappa_missing_fixture: SharedGeigerDiscreteFixture
    discrete_delta_boundary_fixture: SharedGeigerDiscreteFixture
    discrete_early_burst_early_fixture: SharedGeigerDiscreteFixture
    discrete_early_burst_weak_fixture: SharedGeigerDiscreteFixture
    discrete_early_burst_late_fixture: SharedGeigerDiscreteFixture
    discrete_early_burst_missing_fixture: SharedGeigerDiscreteFixture
    discrete_er_missing_fixture: SharedGeigerDiscreteFixture
    discrete_er_mismatch_fixture: SharedGeigerDiscreteFixture
    discrete_sym_three_state_fixture: SharedGeigerDiscreteFixture
    discrete_sym_four_state_fixture: SharedGeigerDiscreteFixture
    discrete_sym_missing_fixture: SharedGeigerDiscreteFixture
    discrete_ard_binary_fixture: SharedGeigerDiscreteFixture
    discrete_ard_four_state_fixture: SharedGeigerDiscreteFixture
    discrete_ard_missing_fixture: SharedGeigerDiscreteFixture


def _package_root() -> Path:
    return Path(__file__).resolve().parents[5]


def build_geiger_registry_fixture_catalog() -> GeigerRegistryFixtureCatalog:
    """Resolve the shared fixtures consumed by the governed `geiger` parity cases."""

    package_root = _package_root()
    tests_root = package_root / "tests" / "fixtures"
    return GeigerRegistryFixtureCatalog(
        tests_root=tests_root,
        brownian_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        brownian_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        white_strong_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        white_weak_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_white_noise_twenty_four_taxa"
        ),
        white_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        lambda_strong_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        lambda_weak_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_white_noise_twenty_four_taxa"
        ),
        lambda_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        kappa_strong_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        kappa_weak_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_white_noise_twenty_four_taxa"
        ),
        kappa_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        delta_strong_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        delta_weak_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_white_noise_twenty_four_taxa"
        ),
        delta_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        ou_missing_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_missing_values_twenty_four_taxa"
        ),
        ou_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_ou_known_truth_twenty_four_taxa"
        ),
        ou_lower_boundary_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_nonultrametric_control_twenty_four_taxa"
        ),
        early_burst_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
        ),
        early_burst_boundary_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        comparison_brownian_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_brownian_signal_twenty_four_taxa"
        ),
        comparison_ou_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_ou_known_truth_twenty_four_taxa"
        ),
        comparison_early_burst_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
        ),
        comparison_white_fixture=get_shared_geiger_continuous_fixture(
            "geiger_continuous_white_noise_twenty_four_taxa"
        ),
        discrete_er_binary_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_er_binary_twenty_four_taxa"
        ),
        discrete_transform_weak_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_transform_weak_signal_twenty_four_taxa"
        ),
        discrete_lambda_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_lambda_missing_binary_twenty_four_taxa"
        ),
        discrete_kappa_strong_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_kappa_branch_sensitive_twenty_four_taxa"
        ),
        discrete_kappa_weak_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_kappa_weak_signal_twenty_four_taxa"
        ),
        discrete_kappa_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_kappa_missing_three_state_twenty_four_taxa"
        ),
        discrete_delta_boundary_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_delta_late_change_binary_twenty_four_taxa"
        ),
        discrete_early_burst_early_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_early_burst_early_change_twenty_four_taxa"
        ),
        discrete_early_burst_weak_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_early_burst_weak_signal_twenty_four_taxa"
        ),
        discrete_early_burst_late_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_early_burst_late_change_twenty_four_taxa"
        ),
        discrete_early_burst_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_early_burst_missing_binary_twenty_four_taxa"
        ),
        discrete_er_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_missing_three_state_twenty_four_taxa"
        ),
        discrete_er_mismatch_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_mismatch_four_state_twenty_four_taxa"
        ),
        discrete_sym_three_state_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_sym_three_state_twenty_four_taxa"
        ),
        discrete_sym_four_state_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_sym_four_state_twenty_four_taxa"
        ),
        discrete_sym_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_missing_three_state_twenty_four_taxa"
        ),
        discrete_ard_binary_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_ard_binary_twenty_four_taxa"
        ),
        discrete_ard_four_state_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_ard_four_state_twenty_four_taxa"
        ),
        discrete_ard_missing_fixture=get_shared_geiger_discrete_fixture(
            "geiger_discrete_missing_three_state_twenty_four_taxa"
        ),
    )
