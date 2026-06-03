from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state,
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
    strip_partition_model_parameter_state,
)
from bijux_phylogenetics.bayesian.state import build_bayesian_model_parameter_state
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_partition_model_state_round_trip_preserves_linkage_and_partition_values() -> (
    None
):
    partition_models = (
        build_partition_substitution_model_definition(
            partition_name="gene_alpha",
            model_name="HKY85+G",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_beta",
            model_name="HKY85+G",
        ),
    )
    linkage_plan = build_partition_parameter_linkage_plan(
        partition_names=("gene_alpha", "gene_beta"),
        linkage_policies={
            "kappa": "linked",
            "base-frequencies": "unlinked",
            "gamma-alpha": "linked",
        },
    )

    model_parameter_state = build_partition_model_parameter_state(
        partition_models=partition_models,
        linkage_plan=linkage_plan,
        partition_parameter_states=(
            _partition_state(
                partition_name="gene_alpha",
                kappa=2.0,
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            _partition_state(
                partition_name="gene_beta",
                kappa=2.0,
                base_frequencies=(0.25, 0.15, 0.25, 0.35),
                gamma_alpha=0.7,
            ),
        ),
        preserved_categorical_parameters={"substitution-model": "partitioned-dna"},
        preserved_scalar_parameters={"clock-rate": 0.5},
    )

    resolved_linkage_plan = (
        resolve_partition_parameter_linkage_plan_from_model_parameters(
            model_parameters=model_parameter_state,
            partition_names=("gene_alpha", "gene_beta"),
        )
    )
    resolved_states = resolve_partition_parameter_states_from_model_parameters(
        model_parameters=model_parameter_state,
        partition_models=partition_models,
        linkage_plan=resolved_linkage_plan,
    )

    assert (
        model_parameter_state.categorical_parameters["substitution-model"]
        == "partitioned-dna"
    )
    assert model_parameter_state.scalar_parameters["clock-rate"] == 0.5
    assert resolved_linkage_plan.groups_for_target("kappa") == {
        "gene_alpha": "kappa-shared",
        "gene_beta": "kappa-shared",
    }
    assert resolved_linkage_plan.groups_for_target("base-frequencies") == {
        "gene_alpha": "gene_alpha",
        "gene_beta": "gene_beta",
    }
    assert resolved_states[0].kappa == 2.0
    assert resolved_states[1].kappa == 2.0
    assert resolved_states[0].base_frequencies == {
        "A": 0.3,
        "C": 0.2,
        "G": 0.1,
        "T": 0.4,
    }
    assert resolved_states[1].base_frequencies == {
        "A": 0.25,
        "C": 0.15,
        "G": 0.25,
        "T": 0.35,
    }


def test_strip_partition_model_parameter_state_removes_reserved_keys_only() -> None:
    model_parameter_state = build_bayesian_model_parameter_state(
        categorical_parameters={
            "substitution-model": "partitioned-dna",
            "partition-linkage:kappa:gene_alpha": "kappa-shared",
        },
        scalar_parameters={
            "clock-rate": 0.5,
            "partition-parameter:kappa:kappa-shared": 2.0,
        },
        vector_parameters={
            "partition-parameter:base-frequencies:gene_alpha": {
                "A": 0.3,
                "C": 0.2,
                "G": 0.1,
                "T": 0.4,
            },
            "exchangeabilities": {
                "AC": 1.0,
                "AG": 2.0,
                "AT": 1.0,
                "CG": 3.0,
                "CT": 2.0,
                "GT": 1.0,
            },
        },
    )

    stripped = strip_partition_model_parameter_state(model_parameter_state)

    assert stripped.categorical_parameters == {"substitution-model": "partitioned-dna"}
    assert stripped.scalar_parameters == {"clock-rate": 0.5}
    assert stripped.vector_parameters == {
        "exchangeabilities": {
            "AC": 1.0,
            "AG": 2.0,
            "AT": 1.0,
            "CG": 3.0,
            "CT": 2.0,
            "GT": 1.0,
        }
    }


def test_partition_model_state_encoding_uses_partition_names_not_state_order() -> None:
    partition_models = (
        build_partition_substitution_model_definition(
            partition_name="gene_alpha",
            model_name="HKY85",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_beta",
            model_name="HKY85",
        ),
    )
    linkage_plan = build_partition_parameter_linkage_plan(
        partition_names=("gene_alpha", "gene_beta"),
        linkage_policies={
            "kappa": "unlinked",
            "base-frequencies": "unlinked",
        },
    )

    model_parameter_state = build_partition_model_parameter_state(
        partition_models=partition_models,
        linkage_plan=linkage_plan,
        partition_parameter_states=(
            _partition_state(
                partition_name="gene_beta",
                kappa=3.0,
                base_frequencies=(0.25, 0.15, 0.25, 0.35),
            ),
            _partition_state(
                partition_name="gene_alpha",
                kappa=2.0,
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
            ),
        ),
    )

    resolved_states = resolve_partition_parameter_states_from_model_parameters(
        model_parameters=model_parameter_state,
        partition_models=partition_models,
        linkage_plan=linkage_plan,
    )

    assert [state.partition_name for state in resolved_states] == [
        "gene_alpha",
        "gene_beta",
    ]
    assert resolved_states[0].kappa == 2.0
    assert resolved_states[1].kappa == 3.0


def test_partition_model_state_encoding_rejects_duplicate_partition_state_entries() -> (
    None
):
    partition_models = (
        build_partition_substitution_model_definition(
            partition_name="gene_alpha",
            model_name="HKY85",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_beta",
            model_name="HKY85",
        ),
    )
    linkage_plan = build_partition_parameter_linkage_plan(
        partition_names=("gene_alpha", "gene_beta"),
    )

    with pytest.raises(
        PhylogeneticsError,
        match="received more than one state for the same partition",
    ):
        build_partition_model_parameter_state(
            partition_models=partition_models,
            linkage_plan=linkage_plan,
            partition_parameter_states=(
                _partition_state(
                    partition_name="gene_alpha",
                    kappa=2.0,
                    base_frequencies=(0.3, 0.2, 0.1, 0.4),
                ),
                _partition_state(
                    partition_name="gene_alpha",
                    kappa=2.0,
                    base_frequencies=(0.3, 0.2, 0.1, 0.4),
                ),
            ),
        )


def _partition_state(
    *,
    partition_name: str,
    kappa: float | None = None,
    base_frequencies: tuple[float, float, float, float] | None = None,
    gamma_alpha: float | None = None,
):
    from bijux_phylogenetics.bayesian.partition_model_priors import (
        PartitionSubstitutionParameterState,
    )

    return PartitionSubstitutionParameterState(
        partition_name=partition_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        gamma_alpha=gamma_alpha,
    )
