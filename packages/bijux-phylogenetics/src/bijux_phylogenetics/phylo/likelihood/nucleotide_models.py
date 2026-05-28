from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
    estimate_empirical_dna_base_frequencies,
    normalize_dna_exchangeabilities_by_anchor,
    validate_dna_base_frequencies,
    validate_positive_kappa,
)
from bijux_phylogenetics.phylo.likelihood.f81 import (
    f81_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    gtr_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    k80_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_root_priors import (
    resolve_nucleotide_root_prior,
)

_SELECTED_NUCLEOTIDE_MODELS = frozenset({"jc69", "k80", "f81", "hky85", "gtr"})


@dataclass(frozen=True, slots=True)
class SelectedNucleotideLikelihoodSpecification:
    """Resolved model parameters for one fixed-topology nucleotide likelihood."""

    model_name: str
    root_prior: numpy.ndarray
    root_prior_source: str
    parameter_values: dict[str, float]
    transition_matrix_for_branch_length: Callable[[float], numpy.ndarray]


def validate_selected_nucleotide_likelihood_model(model_name: str) -> str:
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _SELECTED_NUCLEOTIDE_MODELS:
        raise ValueError(
            "selected nucleotide likelihood model must be one of "
            + ", ".join(sorted(_SELECTED_NUCLEOTIDE_MODELS))
        )
    return normalized_model_name


def resolve_selected_nucleotide_likelihood_specification(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
) -> SelectedNucleotideLikelihoodSpecification:
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    if normalized_model_name == "jc69":
        _reject_irrelevant_parameter(owner_name, "kappa", kappa)
        _reject_irrelevant_parameter(owner_name, "base_frequencies", base_frequencies)
        _reject_irrelevant_parameter(owner_name, "exchangeabilities", exchangeabilities)
        resolved_root_prior = resolve_nucleotide_root_prior(
            records,
            owner_name=owner_name,
            default_policy="equal",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            stationary_frequencies=UNIFORM_DNA_ROOT_PRIOR,
        )
        return SelectedNucleotideLikelihoodSpecification(
            model_name="JC69",
            root_prior=resolved_root_prior.root_prior,
            root_prior_source=resolved_root_prior.root_prior_source,
            parameter_values={},
            transition_matrix_for_branch_length=jc69_transition_probability_matrix,
        )
    if normalized_model_name == "k80":
        _reject_irrelevant_parameter(owner_name, "base_frequencies", base_frequencies)
        _reject_irrelevant_parameter(owner_name, "exchangeabilities", exchangeabilities)
        if kappa is None:
            raise ValueError(f"{owner_name} requires 'kappa'")
        validated_kappa = validate_positive_kappa(kappa, model_name="K80")
        resolved_root_prior = resolve_nucleotide_root_prior(
            records,
            owner_name=owner_name,
            default_policy="equal",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            stationary_frequencies=UNIFORM_DNA_ROOT_PRIOR,
        )
        return SelectedNucleotideLikelihoodSpecification(
            model_name="K80",
            root_prior=resolved_root_prior.root_prior,
            root_prior_source=resolved_root_prior.root_prior_source,
            parameter_values={"kappa": validated_kappa},
            transition_matrix_for_branch_length=lambda branch_length: (
                k80_transition_probability_matrix(
                    branch_length,
                    kappa=validated_kappa,
                )
            ),
        )
    if normalized_model_name == "f81":
        _reject_irrelevant_parameter(owner_name, "kappa", kappa)
        _reject_irrelevant_parameter(owner_name, "exchangeabilities", exchangeabilities)
        stationary = _resolve_dna_base_frequencies(
            records,
            base_frequencies=base_frequencies,
            model_name="F81",
        )
        resolved_root_prior = resolve_nucleotide_root_prior(
            records,
            owner_name=owner_name,
            default_policy="stationary",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            stationary_frequencies=stationary,
        )
        return SelectedNucleotideLikelihoodSpecification(
            model_name="F81",
            root_prior=resolved_root_prior.root_prior,
            root_prior_source=resolved_root_prior.root_prior_source,
            parameter_values=_base_frequency_parameter_values(stationary),
            transition_matrix_for_branch_length=lambda branch_length: (
                f81_transition_probability_matrix(
                    branch_length,
                    base_frequencies=stationary,
                )
            ),
        )
    if normalized_model_name == "hky85":
        _reject_irrelevant_parameter(owner_name, "exchangeabilities", exchangeabilities)
        if kappa is None:
            raise ValueError(f"{owner_name} requires 'kappa'")
        stationary = _resolve_dna_base_frequencies(
            records,
            base_frequencies=base_frequencies,
            model_name="HKY85",
        )
        validated_kappa = validate_positive_kappa(kappa, model_name="HKY85")
        resolved_root_prior = resolve_nucleotide_root_prior(
            records,
            owner_name=owner_name,
            default_policy="stationary",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            stationary_frequencies=stationary,
        )
        return SelectedNucleotideLikelihoodSpecification(
            model_name="HKY85",
            root_prior=resolved_root_prior.root_prior,
            root_prior_source=resolved_root_prior.root_prior_source,
            parameter_values={
                **_base_frequency_parameter_values(stationary),
                "kappa": validated_kappa,
            },
            transition_matrix_for_branch_length=lambda branch_length: (
                hky85_transition_probability_matrix(
                    branch_length,
                    base_frequencies=stationary,
                    kappa=validated_kappa,
                )
            ),
        )

    _reject_irrelevant_parameter(owner_name, "kappa", kappa)
    if exchangeabilities is None:
        raise ValueError(f"{owner_name} requires 'exchangeabilities'")
    stationary = _resolve_dna_base_frequencies(
        records,
        base_frequencies=base_frequencies,
        model_name="GTR",
    )
    normalized_exchangeabilities = normalize_dna_exchangeabilities_by_anchor(
        exchangeabilities,
        model_name="GTR",
    )
    resolved_root_prior = resolve_nucleotide_root_prior(
        records,
        owner_name=owner_name,
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary,
    )
    return SelectedNucleotideLikelihoodSpecification(
        model_name="GTR",
        root_prior=resolved_root_prior.root_prior,
        root_prior_source=resolved_root_prior.root_prior_source,
        parameter_values={
            **_base_frequency_parameter_values(stationary),
            "exchangeability_ac": float(normalized_exchangeabilities[0]),
            "exchangeability_ag": float(normalized_exchangeabilities[1]),
            "exchangeability_at": float(normalized_exchangeabilities[2]),
            "exchangeability_cg": float(normalized_exchangeabilities[3]),
            "exchangeability_ct": float(normalized_exchangeabilities[4]),
            "exchangeability_gt": float(normalized_exchangeabilities[5]),
        },
        transition_matrix_for_branch_length=lambda branch_length: (
            gtr_transition_probability_matrix(
                branch_length,
                exchangeabilities=normalized_exchangeabilities,
                base_frequencies=stationary,
            )
        ),
    )


def _resolve_dna_base_frequencies(
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    model_name: str,
) -> numpy.ndarray:
    if base_frequencies is None:
        return estimate_empirical_dna_base_frequencies(records)
    return validate_dna_base_frequencies(base_frequencies, model_name=model_name)


def _base_frequency_parameter_values(
    stationary: numpy.ndarray,
) -> dict[str, float]:
    return {
        "base_frequency_a": float(stationary[0]),
        "base_frequency_c": float(stationary[1]),
        "base_frequency_g": float(stationary[2]),
        "base_frequency_t": float(stationary[3]),
    }


def _reject_irrelevant_parameter(
    owner_name: str,
    parameter_name: str,
    value: object,
) -> None:
    if value is not None:
        raise ValueError(
            f"{owner_name} does not accept '{parameter_name}' because that model does not use it"
        )
