from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures import (
    SharedPhytoolsComparativeFixture,
    get_shared_phytools_comparative_fixture,
)


@dataclass(frozen=True, slots=True)
class PhytoolsRegistryFixtureCatalog:
    """Resolved fixture catalog for the governed `phytools` parity registry."""

    simulation_tree_fixture: Path
    simulation_six_taxa_tree_fixture: Path
    strong_signal_fixture: SharedPhytoolsComparativeFixture
    nonultrametric_signal_fixture: SharedPhytoolsComparativeFixture
    weak_signal_fixture: SharedPhytoolsComparativeFixture
    missing_signal_fixture: SharedPhytoolsComparativeFixture
    pgls_continuous_fixture: SharedPhytoolsComparativeFixture
    pgls_categorical_fixture: SharedPhytoolsComparativeFixture
    pgls_interaction_fixture: SharedPhytoolsComparativeFixture
    phyl_resid_brownian_fixture: SharedPhytoolsComparativeFixture
    phyl_resid_lambda_fixture: SharedPhytoolsComparativeFixture
    phyl_resid_lambda_missing_fixture: SharedPhytoolsComparativeFixture
    phyl_anova_fixture: SharedPhytoolsComparativeFixture
    phyl_anova_missing_fixture: SharedPhytoolsComparativeFixture
    binary_discrete_fixture: SharedPhytoolsComparativeFixture
    multistate_discrete_fixture: SharedPhytoolsComparativeFixture
    binary_discrete_missing_fixture: SharedPhytoolsComparativeFixture
    multistate_discrete_missing_fixture: SharedPhytoolsComparativeFixture
    ard_binary_discrete_fixture: SharedPhytoolsComparativeFixture
    ard_multistate_discrete_fixture: SharedPhytoolsComparativeFixture
    ard_binary_discrete_missing_fixture: SharedPhytoolsComparativeFixture
    ard_multistate_discrete_missing_fixture: SharedPhytoolsComparativeFixture


def _package_root() -> Path:
    return Path(__file__).resolve().parents[5]


def build_phytools_registry_fixture_catalog() -> PhytoolsRegistryFixtureCatalog:
    """Resolve the shared fixtures consumed by the governed `phytools` parity cases."""

    tests_root = _package_root() / "tests" / "fixtures" / "trees"
    return PhytoolsRegistryFixtureCatalog(
        simulation_tree_fixture=tests_root / "example_tree.nwk",
        simulation_six_taxa_tree_fixture=tests_root / "example_tree_six_taxa.nwk",
        strong_signal_fixture=get_shared_phytools_comparative_fixture(
            "phytools_continuous_strong_signal_twenty_four_taxa"
        ),
        nonultrametric_signal_fixture=get_shared_phytools_comparative_fixture(
            "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
        ),
        weak_signal_fixture=get_shared_phytools_comparative_fixture(
            "phytools_continuous_weak_signal_twenty_four_taxa"
        ),
        missing_signal_fixture=get_shared_phytools_comparative_fixture(
            "phytools_continuous_missing_values_twenty_four_taxa"
        ),
        pgls_continuous_fixture=get_shared_phytools_comparative_fixture(
            "phytools_pgls_brownian_continuous_four_taxa"
        ),
        pgls_categorical_fixture=get_shared_phytools_comparative_fixture(
            "phytools_pgls_brownian_categorical_eight_taxa"
        ),
        pgls_interaction_fixture=get_shared_phytools_comparative_fixture(
            "phytools_pgls_brownian_interaction_eight_taxa"
        ),
        phyl_resid_brownian_fixture=get_shared_phytools_comparative_fixture(
            "phytools_phyl_resid_bm_allometry_six_taxa"
        ),
        phyl_resid_lambda_fixture=get_shared_phytools_comparative_fixture(
            "phytools_phyl_resid_lambda_allometry_six_taxa"
        ),
        phyl_resid_lambda_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_phyl_resid_lambda_missing_six_taxa"
        ),
        phyl_anova_fixture=get_shared_phytools_comparative_fixture(
            "phytools_phyl_anova_group_effect_six_taxa"
        ),
        phyl_anova_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_phyl_anova_group_effect_missing_six_taxa"
        ),
        binary_discrete_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_binary_twenty_four_taxa"
        ),
        multistate_discrete_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_multistate_twenty_four_taxa"
        ),
        binary_discrete_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_binary_missing_twenty_four_taxa"
        ),
        multistate_discrete_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_multistate_missing_twenty_four_taxa"
        ),
        ard_binary_discrete_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_ard_binary_twenty_four_taxa"
        ),
        ard_multistate_discrete_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_ard_multistate_twenty_four_taxa"
        ),
        ard_binary_discrete_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_ard_binary_missing_twenty_four_taxa"
        ),
        ard_multistate_discrete_missing_fixture=get_shared_phytools_comparative_fixture(
            "phytools_discrete_ard_multistate_missing_twenty_four_taxa"
        ),
    )
