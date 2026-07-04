from __future__ import annotations

from collections.abc import Mapping, Sequence
import math

from bijux_phylogenetics.ancestral.discrete import DiscreteTransitionRateRow
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
    evaluate_discrete_trait_rate_log_prior,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_codon_frequency_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_gamma_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import resolve_codon_state_space
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
)

from .contracts import (
    SubstitutionParameterPriorCoverageEntry,
    SubstitutionParameterPriorCoverageReport,
)

_LIKELIHOOD_TOLERANCE = 1e-12
_DNA_BASE_COMPONENT_NAMES = ("A", "C", "G", "T")


def build_substitution_parameter_prior_coverage() -> (
    SubstitutionParameterPriorCoverageReport
):
    """Prove that each supported substitution-prior target changes prior scores."""

    entries = [
        _build_kappa_coverage_entry(),
        _build_exchangeability_coverage_entry(),
        _build_base_frequency_coverage_entry(),
        _build_gamma_alpha_coverage_entry(),
        _build_invariant_proportion_coverage_entry(),
        _build_codon_frequency_coverage_entry(),
        _build_mk_rate_coverage_entry(),
    ]
    passed_target_count = sum(1 for entry in entries if entry.passed)
    return SubstitutionParameterPriorCoverageReport(
        goal_id=76,
        passed=passed_target_count == len(entries),
        target_count=len(entries),
        passed_target_count=passed_target_count,
        entries=entries,
        limitations=[
            "the coverage report isolates one governed score comparison per supported prior target; it does not attempt whole-model MCMC calibration or posterior predictive validation",
            "the substitution-parameter cases use explicit parameter values rather than stochastic proposals so score sensitivity is attributable to one configured prior target at a time",
            "the Mk-rate case proves the owned discrete-trait rate-prior surface on a deterministic transition matrix rather than on a fitted empirical dataset",
        ],
    )


def _build_kappa_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    return _build_bundle_target_coverage_entry(
        target_name="kappa",
        prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=0.5
            )
        ),
        baseline_parameters={"kappa": 0.8},
        comparison_parameters={"kappa": 1.8},
        covered_names=["value"],
        notes=[
            "uses an exponential prior so a higher transition-transversion bias changes the scalar log density",
        ],
    )


def _build_exchangeability_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    return _build_bundle_target_coverage_entry(
        target_name="exchangeabilities",
        prior_bundle=build_substitution_parameter_prior_bundle(
            exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=DNA_EXCHANGEABILITY_LABELS,
                concentration_parameters={
                    "AC": 2.0,
                    "AG": 5.0,
                    "AT": 3.0,
                    "CG": 4.0,
                    "CT": 6.0,
                    "GT": 2.5,
                },
            )
        ),
        baseline_parameters={
            "exchangeabilities": {
                "AC": 0.10,
                "AG": 0.28,
                "AT": 0.12,
                "CG": 0.18,
                "CT": 0.22,
                "GT": 0.10,
            }
        },
        comparison_parameters={
            "exchangeabilities": {
                "AC": 0.18,
                "AG": 0.14,
                "AT": 0.10,
                "CG": 0.12,
                "CT": 0.30,
                "GT": 0.16,
            }
        },
        covered_names=DNA_EXCHANGEABILITY_LABELS,
        notes=[
            "uses a non-uniform Dirichlet prior so the score changes with the realized reversible-rate composition",
        ],
    )


def _build_base_frequency_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    return _build_bundle_target_coverage_entry(
        target_name="base-frequencies",
        prior_bundle=build_substitution_parameter_prior_bundle(
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=_DNA_BASE_COMPONENT_NAMES,
                concentration_parameters={
                    "A": 5.0,
                    "C": 2.0,
                    "G": 3.5,
                    "T": 4.0,
                },
            )
        ),
        baseline_parameters={
            "base_frequencies": {"A": 0.35, "C": 0.15, "G": 0.20, "T": 0.30}
        },
        comparison_parameters={
            "base_frequencies": {"A": 0.20, "C": 0.28, "G": 0.22, "T": 0.30}
        },
        covered_names=_DNA_BASE_COMPONENT_NAMES,
        notes=[
            "uses composition-biased Dirichlet weights so the base-frequency score is not flat over the simplex",
        ],
    )


def _build_gamma_alpha_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    return _build_bundle_target_coverage_entry(
        target_name="gamma-alpha",
        prior_bundle=build_substitution_parameter_prior_bundle(
            gamma_alpha_prior=build_gamma_positive_substitution_parameter_prior(
                shape=2.5,
                scale=0.75,
            )
        ),
        baseline_parameters={"gamma_alpha": 0.8},
        comparison_parameters={"gamma_alpha": 2.2},
        covered_names=["value"],
        notes=[
            "uses a gamma prior so alternative discrete-gamma concentration values produce different prior support",
        ],
    )


def _build_invariant_proportion_coverage_entry() -> (
    SubstitutionParameterPriorCoverageEntry
):
    return _build_bundle_target_coverage_entry(
        target_name="invariant-proportion",
        prior_bundle=build_substitution_parameter_prior_bundle(
            invariant_proportion_prior=build_beta_probability_substitution_parameter_prior(
                alpha=2.0,
                beta=5.0,
            )
        ),
        baseline_parameters={"invariant_proportion": 0.10},
        comparison_parameters={"invariant_proportion": 0.35},
        covered_names=["value"],
        notes=[
            "uses a beta prior so the invariant-mass contribution changes away from the low-proportion mode",
        ],
    )


def _build_codon_frequency_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    codon_component_names = resolve_codon_state_space().state_order
    return _build_bundle_target_coverage_entry(
        target_name="codon-frequencies",
        prior_bundle=build_substitution_parameter_prior_bundle(
            codon_frequency_prior=build_dirichlet_codon_frequency_substitution_parameter_prior(
                concentration_parameters={
                    component_name: 1.0 + ((index % 5) * 0.5)
                    for index, component_name in enumerate(codon_component_names)
                }
            )
        ),
        baseline_parameters={
            "codon_frequencies": _build_codon_frequency_weights(
                codon_component_names,
                leading_weights=(8.0, 5.0, 3.0),
                trailing_weight=1.0,
            )
        },
        comparison_parameters={
            "codon_frequencies": _build_codon_frequency_weights(
                codon_component_names,
                leading_weights=(2.0, 6.0, 9.0),
                trailing_weight=1.0,
            )
        },
        covered_names=codon_component_names,
        notes=[
            "covers the full 61-state codon simplex instead of collapsing codon composition back to nucleotide frequencies",
        ],
    )


def _build_mk_rate_coverage_entry() -> SubstitutionParameterPriorCoverageEntry:
    prior_model = build_exponential_discrete_trait_rate_prior(rate=0.8)
    baseline_rows = [
        DiscreteTransitionRateRow(
            source_state="0",
            target_state="1",
            transition_allowed=True,
            step_distance=1,
            rate=0.4,
        ),
        DiscreteTransitionRateRow(
            source_state="1",
            target_state="0",
            transition_allowed=True,
            step_distance=1,
            rate=0.7,
        ),
    ]
    comparison_rows = [
        DiscreteTransitionRateRow(
            source_state="0",
            target_state="1",
            transition_allowed=True,
            step_distance=1,
            rate=1.2,
        ),
        DiscreteTransitionRateRow(
            source_state="1",
            target_state="0",
            transition_allowed=True,
            step_distance=1,
            rate=0.9,
        ),
    ]
    baseline_report = evaluate_discrete_trait_rate_log_prior(
        model="all-rates-different",
        transition_rate_rows=baseline_rows,
        prior_model=prior_model,
    )
    comparison_report = evaluate_discrete_trait_rate_log_prior(
        model="all-rates-different",
        transition_rate_rows=comparison_rows,
        prior_model=prior_model,
    )
    changed = not math.isclose(
        baseline_report.total_log_prior,
        comparison_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=_LIKELIHOOD_TOLERANCE,
    )
    covered_names = [row.parameter_name for row in baseline_report.rows]
    passed = (
        changed
        and baseline_report.parameter_count == 2
        and comparison_report.parameter_count == 2
        and covered_names == [row.parameter_name for row in comparison_report.rows]
    )
    return SubstitutionParameterPriorCoverageEntry(
        target_name="mk-rates",
        surface_name="discrete-trait-rate-prior",
        prior_family=baseline_report.family,
        contribution_count=baseline_report.parameter_count,
        covered_names=covered_names,
        baseline_log_prior=baseline_report.total_log_prior,
        comparison_log_prior=comparison_report.total_log_prior,
        log_prior_delta=float(
            format(
                comparison_report.total_log_prior - baseline_report.total_log_prior,
                ".15g",
            )
        ),
        changed=changed,
        passed=passed,
        notes=[
            "uses the owned all-rates-different Mk parameterization so both directional rates carry explicit prior mass",
        ],
    )


def _build_bundle_target_coverage_entry(
    *,
    target_name: str,
    prior_bundle: object,
    baseline_parameters: Mapping[str, object],
    comparison_parameters: Mapping[str, object],
    covered_names: Sequence[str],
    notes: list[str],
) -> SubstitutionParameterPriorCoverageEntry:
    baseline_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        **baseline_parameters,
    )
    comparison_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        **comparison_parameters,
    )
    changed = not math.isclose(
        baseline_report.total_log_prior,
        comparison_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=_LIKELIHOOD_TOLERANCE,
    )
    passed = (
        changed
        and baseline_report.prior_count == 1
        and comparison_report.prior_count == 1
        and len(baseline_report.rows) == 1
        and len(comparison_report.rows) == 1
        and baseline_report.rows[0].target_name == target_name
        and comparison_report.rows[0].target_name == target_name
    )
    return SubstitutionParameterPriorCoverageEntry(
        target_name=target_name,
        surface_name="substitution-parameter-bundle",
        prior_family=baseline_report.rows[0].family,
        contribution_count=len(baseline_report.rows),
        covered_names=list(covered_names),
        baseline_log_prior=baseline_report.total_log_prior,
        comparison_log_prior=comparison_report.total_log_prior,
        log_prior_delta=float(
            format(
                comparison_report.total_log_prior - baseline_report.total_log_prior,
                ".15g",
            )
        ),
        changed=changed,
        passed=passed,
        notes=notes,
    )


def _build_codon_frequency_weights(
    component_names: Sequence[str],
    *,
    leading_weights: tuple[float, float, float],
    trailing_weight: float,
) -> dict[str, float]:
    weights = dict.fromkeys(component_names, trailing_weight)
    for component_name, weight in zip(
        component_names[: len(leading_weights)],
        leading_weights,
        strict=True,
    ):
        weights[component_name] = weight
    return weights
