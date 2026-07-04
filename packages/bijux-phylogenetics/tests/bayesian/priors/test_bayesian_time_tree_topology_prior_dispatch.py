from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentEpoch,
    build_constant_population_coalescent_tree_prior,
    build_crown_conditioned_birth_death_tree_prior,
    build_crown_conditioned_yule_tree_prior,
    build_skyline_coalescent_tree_prior,
    evaluate_birth_death_tree_log_prior,
    evaluate_constant_population_coalescent_tree_log_prior,
    evaluate_skyline_coalescent_tree_log_prior,
    evaluate_yule_tree_log_prior,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    TIME_TREE_TOPOLOGY_PRIOR_FAMILIES,
    evaluate_time_tree_topology_log_prior,
    require_time_tree_topology_prior_model,
    supports_time_tree_topology_prior_model,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("prior_model", "expected_log_prior"),
    [
        (
            build_crown_conditioned_yule_tree_prior(speciation_rate=0.5),
            evaluate_yule_tree_log_prior,
        ),
        (
            build_crown_conditioned_birth_death_tree_prior(
                speciation_rate=1.2,
                extinction_rate=0.4,
                sampling_fraction=0.75,
            ),
            evaluate_birth_death_tree_log_prior,
        ),
        (
            build_constant_population_coalescent_tree_prior(
                effective_population_size=2.0
            ),
            evaluate_constant_population_coalescent_tree_log_prior,
        ),
        (
            build_skyline_coalescent_tree_prior(
                [
                    SkylineCoalescentEpoch(
                        younger_boundary_age=0.0,
                        older_boundary_age=1.5,
                        effective_population_size=2.0,
                    ),
                    SkylineCoalescentEpoch(
                        younger_boundary_age=1.5,
                        older_boundary_age=None,
                        effective_population_size=4.0,
                    ),
                ]
            ),
            evaluate_skyline_coalescent_tree_log_prior,
        ),
    ],
)
def test_shared_time_tree_topology_prior_dispatch_matches_family_evaluators(
    prior_model: object,
    expected_log_prior,
) -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))

    shared_report = evaluate_time_tree_topology_log_prior(tree, prior_model)
    specific_report = expected_log_prior(tree, prior_model)

    assert shared_report.family == specific_report.family
    assert shared_report.log_prior == specific_report.log_prior


def test_time_tree_topology_prior_support_guard_accepts_supported_models() -> None:
    assert TIME_TREE_TOPOLOGY_PRIOR_FAMILIES == (
        "constant-population-coalescent",
        "crown-conditioned-birth-death",
        "crown-conditioned-yule",
        "skyline-coalescent",
    )
    assert supports_time_tree_topology_prior_model(
        build_crown_conditioned_yule_tree_prior(speciation_rate=0.5)
    )
    assert supports_time_tree_topology_prior_model(
        build_crown_conditioned_birth_death_tree_prior(
            speciation_rate=1.0,
            extinction_rate=0.2,
            sampling_fraction=0.8,
        )
    )
    assert supports_time_tree_topology_prior_model(
        build_constant_population_coalescent_tree_prior(effective_population_size=3.0)
    )
    assert supports_time_tree_topology_prior_model(
        build_skyline_coalescent_tree_prior([SkylineCoalescentEpoch(0.0, None, 3.0)])
    )
    assert not supports_time_tree_topology_prior_model(object())


def test_time_tree_topology_prior_dispatch_rejects_unsupported_model() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))

    with pytest.raises(PhylogeneticsError) as error_info:
        evaluate_time_tree_topology_log_prior(tree, object())

    assert error_info.value.code == "time_tree_topology_prior_type_invalid"


def test_time_tree_topology_prior_requirement_uses_caller_message() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        require_time_tree_topology_prior_model(
            object(),
            owner_name="strict-clock topology surface",
            code="strict_clock_topology_prior_invalid",
        )

    assert str(error_info.value) == (
        "strict-clock topology surface requires one supported rooted time-tree prior model"
    )
    assert error_info.value.code == "strict_clock_topology_prior_invalid"
