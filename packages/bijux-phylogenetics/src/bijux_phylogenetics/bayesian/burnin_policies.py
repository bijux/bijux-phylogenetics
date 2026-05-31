from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsRunReport,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.posterior_sets.diagnostics import (
    standardized_mean_shift,
)
from bijux_phylogenetics.bayesian.state import BayesianPhylogeneticState
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES = (
    "none",
    "fixed-count",
    "fixed-fraction",
    "diagnostic-suggested",
)


@dataclass(frozen=True, slots=True)
class MetropolisHastingsBurninPolicy:
    """One validated burn-in policy over one sampled Metropolis-Hastings trace."""

    policy_name: str
    discarded_sample_count: int | None = None
    discarded_fraction: float | None = None
    mean_shift_threshold: float | None = None
    minimum_retained_sample_count: int | None = None


@dataclass(frozen=True, slots=True)
class BurninSampleRow:
    """One sampled state row retained or discarded by one burn-in policy."""

    sample_index: int
    iteration_index: int
    posterior_log_score: float
    state: BayesianPhylogeneticState


@dataclass(frozen=True, slots=True)
class MetropolisHastingsBurninReport:
    """One applied burn-in policy over one native Metropolis-Hastings trace."""

    policy: MetropolisHastingsBurninPolicy
    sample_every: int
    total_sample_count: int
    discarded_sample_count: int
    retained_sample_count: int
    discarded_rows: list[BurninSampleRow]
    retained_rows: list[BurninSampleRow]
    diagnostic_report: MetropolisHastingsBurninDiagnosticReport | None = None


@dataclass(frozen=True, slots=True)
class MetropolisHastingsBurninDiagnosticCandidate:
    """One candidate retained tail evaluated by the diagnostic burn-in heuristic."""

    discarded_sample_count: int
    retained_sample_count: int
    maximum_mean_shift: float
    parameter_mean_shifts: dict[str, float]
    acceptable: bool


@dataclass(frozen=True, slots=True)
class MetropolisHastingsBurninDiagnosticReport:
    """One diagnostic suggestion ledger over one native sampled trace."""

    mean_shift_threshold: float
    minimum_retained_sample_count: int
    selected_discarded_sample_count: int
    stabilized_tail_found: bool
    candidate_rows: list[MetropolisHastingsBurninDiagnosticCandidate]


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainBurninReport:
    """One named applied burn-in policy over one independent chain."""

    chain_name: str
    burnin_report: MetropolisHastingsBurninReport


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsBurninReport:
    """One applied burn-in policy over one independent-chain run collection."""

    policy: MetropolisHastingsBurninPolicy
    chain_reports: list[IndependentMetropolisHastingsChainBurninReport]
    minimum_retained_sample_count: int
    maximum_retained_sample_count: int


def build_metropolis_hastings_burnin_policy(
    *,
    policy_name: str,
    discarded_sample_count: int | None = None,
    discarded_fraction: float | None = None,
    mean_shift_threshold: float | None = None,
    minimum_retained_sample_count: int | None = None,
) -> MetropolisHastingsBurninPolicy:
    """Build one validated burn-in policy over one sampled trace."""
    validated_policy_name = _validate_policy_name(policy_name)
    validated_discarded_sample_count = _validate_optional_nonnegative_integer(
        value=discarded_sample_count,
        field_name="discarded_sample_count",
    )
    validated_discarded_fraction = _validate_optional_fraction(
        value=discarded_fraction,
        field_name="discarded_fraction",
    )
    validated_mean_shift_threshold = _validate_optional_positive_float(
        value=mean_shift_threshold,
        field_name="mean_shift_threshold",
    )
    validated_minimum_retained_sample_count = _validate_optional_positive_integer(
        value=minimum_retained_sample_count,
        field_name="minimum_retained_sample_count",
    )
    if validated_policy_name == "none":
        _reject_unexpected_fixed_count_fields(
            discarded_sample_count=validated_discarded_sample_count,
            discarded_fraction=validated_discarded_fraction,
            mean_shift_threshold=validated_mean_shift_threshold,
            minimum_retained_sample_count=validated_minimum_retained_sample_count,
            policy_name=validated_policy_name,
        )
    elif validated_policy_name == "fixed-count":
        if validated_discarded_sample_count is None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-count' requires one discarded_sample_count",
                code="metropolis_hastings_burnin_policy_discarded_sample_count_missing",
            )
        if validated_discarded_fraction is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-count' does not accept discarded_fraction",
                code="metropolis_hastings_burnin_policy_discarded_fraction_unexpected",
            )
        if validated_mean_shift_threshold is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-count' does not accept mean_shift_threshold",
                code="metropolis_hastings_burnin_policy_mean_shift_threshold_unexpected",
            )
        if validated_minimum_retained_sample_count is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-count' does not accept minimum_retained_sample_count",
                code="metropolis_hastings_burnin_policy_minimum_retained_sample_count_unexpected",
            )
    elif validated_policy_name == "fixed-fraction":
        if validated_discarded_fraction is None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-fraction' requires one discarded_fraction",
                code="metropolis_hastings_burnin_policy_discarded_fraction_missing",
            )
        if validated_discarded_sample_count is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-fraction' does not accept discarded_sample_count",
                code="metropolis_hastings_burnin_policy_discarded_sample_count_unexpected",
            )
        if validated_mean_shift_threshold is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-fraction' does not accept mean_shift_threshold",
                code="metropolis_hastings_burnin_policy_mean_shift_threshold_unexpected",
            )
        if validated_minimum_retained_sample_count is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'fixed-fraction' does not accept minimum_retained_sample_count",
                code="metropolis_hastings_burnin_policy_minimum_retained_sample_count_unexpected",
            )
    else:
        if validated_discarded_sample_count is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'diagnostic-suggested' does not accept discarded_sample_count",
                code="metropolis_hastings_burnin_policy_discarded_sample_count_unexpected",
            )
        if validated_discarded_fraction is not None:
            raise PhylogeneticsError(
                "metropolis-hastings burn-in policy 'diagnostic-suggested' does not accept discarded_fraction",
                code="metropolis_hastings_burnin_policy_discarded_fraction_unexpected",
            )
        if validated_mean_shift_threshold is None:
            validated_mean_shift_threshold = 0.5
        if validated_minimum_retained_sample_count is None:
            validated_minimum_retained_sample_count = 4
    return MetropolisHastingsBurninPolicy(
        policy_name=validated_policy_name,
        discarded_sample_count=validated_discarded_sample_count,
        discarded_fraction=validated_discarded_fraction,
        mean_shift_threshold=validated_mean_shift_threshold,
        minimum_retained_sample_count=validated_minimum_retained_sample_count,
    )


def apply_metropolis_hastings_burnin_policy(
    *,
    chain_report: MetropolisHastingsRunReport,
    policy: MetropolisHastingsBurninPolicy,
) -> MetropolisHastingsBurninReport:
    """Apply one validated burn-in policy over one native Metropolis-Hastings run."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings burn-in policy application requires one MetropolisHastingsRunReport",
            code="metropolis_hastings_burnin_chain_report_type_invalid",
        )
    if not isinstance(policy, MetropolisHastingsBurninPolicy):
        raise PhylogeneticsError(
            "metropolis-hastings burn-in policy application requires one MetropolisHastingsBurninPolicy",
            code="metropolis_hastings_burnin_policy_type_invalid",
        )
    sample_rows = _build_sample_rows(chain_report)
    diagnostic_report: MetropolisHastingsBurninDiagnosticReport | None = None
    if policy.policy_name == "none":
        discarded_sample_count = 0
    elif policy.policy_name == "fixed-count":
        discarded_sample_count = policy.discarded_sample_count or 0
    elif policy.policy_name == "fixed-fraction":
        discarded_sample_count = int(
            len(sample_rows) * (policy.discarded_fraction or 0.0)
        )
    else:
        diagnostic_report = diagnose_metropolis_hastings_burnin(
            chain_report=chain_report,
            mean_shift_threshold=policy.mean_shift_threshold or 0.5,
            minimum_retained_sample_count=policy.minimum_retained_sample_count or 4,
        )
        discarded_sample_count = diagnostic_report.selected_discarded_sample_count
    _validate_retained_sample_count(
        total_sample_count=len(sample_rows),
        discarded_sample_count=discarded_sample_count,
    )
    discarded_rows = sample_rows[:discarded_sample_count]
    retained_rows = sample_rows[discarded_sample_count:]
    return MetropolisHastingsBurninReport(
        policy=policy,
        sample_every=chain_report.sample_every,
        total_sample_count=len(sample_rows),
        discarded_sample_count=discarded_sample_count,
        retained_sample_count=len(retained_rows),
        discarded_rows=discarded_rows,
        retained_rows=retained_rows,
        diagnostic_report=diagnostic_report,
    )


def diagnose_metropolis_hastings_burnin(
    *,
    chain_report: MetropolisHastingsRunReport,
    mean_shift_threshold: float = 0.5,
    minimum_retained_sample_count: int = 4,
) -> MetropolisHastingsBurninDiagnosticReport:
    """Diagnose one burn-in suggestion from scalar trace stabilization."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings burn-in diagnostic requires one MetropolisHastingsRunReport",
            code="metropolis_hastings_burnin_diagnostic_chain_report_type_invalid",
        )
    validated_mean_shift_threshold = _validate_optional_positive_float(
        value=mean_shift_threshold,
        field_name="mean_shift_threshold",
    )
    validated_minimum_retained_sample_count = _validate_optional_positive_integer(
        value=minimum_retained_sample_count,
        field_name="minimum_retained_sample_count",
    )
    if validated_mean_shift_threshold is None:
        raise PhylogeneticsError(
            "metropolis-hastings burn-in diagnostic requires one positive mean_shift_threshold",
            code="metropolis_hastings_burnin_diagnostic_mean_shift_threshold_missing",
        )
    if validated_minimum_retained_sample_count is None:
        raise PhylogeneticsError(
            "metropolis-hastings burn-in diagnostic requires one positive minimum_retained_sample_count",
            code="metropolis_hastings_burnin_diagnostic_minimum_retained_sample_count_missing",
        )
    sample_rows = _build_sample_rows(chain_report)
    if validated_minimum_retained_sample_count > len(sample_rows):
        raise PhylogeneticsError(
            "metropolis-hastings burn-in diagnostic requires minimum_retained_sample_count to fit within the sampled trace length",
            code="metropolis_hastings_burnin_diagnostic_minimum_retained_sample_count_too_large",
            details={
                "minimum_retained_sample_count": validated_minimum_retained_sample_count,
                "total_sample_count": len(sample_rows),
            },
        )
    diagnostic_series = _build_diagnostic_series(chain_report)
    candidate_rows: list[MetropolisHastingsBurninDiagnosticCandidate] = []
    maximum_discard_count = len(sample_rows) - validated_minimum_retained_sample_count
    for discarded_sample_count in range(maximum_discard_count + 1):
        parameter_mean_shifts = {
            parameter_name: standardized_mean_shift(values[discarded_sample_count:])
            for parameter_name, values in diagnostic_series.items()
        }
        maximum_mean_shift = max(parameter_mean_shifts.values(), default=0.0)
        candidate_rows.append(
            MetropolisHastingsBurninDiagnosticCandidate(
                discarded_sample_count=discarded_sample_count,
                retained_sample_count=len(sample_rows) - discarded_sample_count,
                maximum_mean_shift=maximum_mean_shift,
                parameter_mean_shifts=parameter_mean_shifts,
                acceptable=maximum_mean_shift <= validated_mean_shift_threshold,
            )
        )
    acceptable_candidates = [
        candidate_row for candidate_row in candidate_rows if candidate_row.acceptable
    ]
    if acceptable_candidates:
        selected_candidate = acceptable_candidates[0]
        stabilized_tail_found = True
    else:
        selected_candidate = min(
            candidate_rows,
            key=lambda candidate_row: (
                candidate_row.maximum_mean_shift,
                candidate_row.discarded_sample_count,
            ),
        )
        stabilized_tail_found = False
    return MetropolisHastingsBurninDiagnosticReport(
        mean_shift_threshold=validated_mean_shift_threshold,
        minimum_retained_sample_count=validated_minimum_retained_sample_count,
        selected_discarded_sample_count=selected_candidate.discarded_sample_count,
        stabilized_tail_found=stabilized_tail_found,
        candidate_rows=candidate_rows,
    )


def apply_independent_metropolis_hastings_burnin_policy(
    *,
    run_report: IndependentMetropolisHastingsRunReport,
    policy: MetropolisHastingsBurninPolicy,
) -> IndependentMetropolisHastingsBurninReport:
    """Apply one burn-in policy over every chain in one independent-chain run."""
    if not isinstance(run_report, IndependentMetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings burn-in policy application requires one IndependentMetropolisHastingsRunReport",
            code="independent_metropolis_hastings_burnin_run_report_type_invalid",
        )
    if not isinstance(policy, MetropolisHastingsBurninPolicy):
        raise PhylogeneticsError(
            "independent metropolis-hastings burn-in policy application requires one MetropolisHastingsBurninPolicy",
            code="independent_metropolis_hastings_burnin_policy_type_invalid",
        )
    chain_reports = [
        IndependentMetropolisHastingsChainBurninReport(
            chain_name=chain_report.chain_name,
            burnin_report=apply_metropolis_hastings_burnin_policy(
                chain_report=chain_report.chain_report,
                policy=policy,
            ),
        )
        for chain_report in run_report.chain_reports
    ]
    retained_sample_counts = [
        chain_report.burnin_report.retained_sample_count
        for chain_report in chain_reports
    ]
    return IndependentMetropolisHastingsBurninReport(
        policy=policy,
        chain_reports=chain_reports,
        minimum_retained_sample_count=min(retained_sample_counts),
        maximum_retained_sample_count=max(retained_sample_counts),
    )


def _reject_unexpected_fixed_count_fields(
    *,
    discarded_sample_count: int | None,
    discarded_fraction: float | None,
    mean_shift_threshold: float | None,
    minimum_retained_sample_count: int | None,
    policy_name: str,
) -> None:
    if discarded_sample_count is not None:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy '{policy_name}' does not accept discarded_sample_count",
            code="metropolis_hastings_burnin_policy_discarded_sample_count_unexpected",
        )
    if discarded_fraction is not None:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy '{policy_name}' does not accept discarded_fraction",
            code="metropolis_hastings_burnin_policy_discarded_fraction_unexpected",
        )
    if mean_shift_threshold is not None:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy '{policy_name}' does not accept mean_shift_threshold",
            code="metropolis_hastings_burnin_policy_mean_shift_threshold_unexpected",
        )
    if minimum_retained_sample_count is not None:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy '{policy_name}' does not accept minimum_retained_sample_count",
            code="metropolis_hastings_burnin_policy_minimum_retained_sample_count_unexpected",
        )


def _validate_policy_name(value: str) -> str:
    validated_value = _validate_nonblank_string(value, field_name="policy_name")
    if validated_value not in METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES:
        raise PhylogeneticsError(
            "metropolis-hastings burn-in policy name is unsupported",
            code="metropolis_hastings_burnin_policy_name_invalid",
            details={
                "policy_name": validated_value,
                "allowed_policy_names": list(METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES),
            },
        )
    return validated_value


def _build_sample_rows(
    chain_report: MetropolisHastingsRunReport,
) -> list[BurninSampleRow]:
    sample_rows = [
        BurninSampleRow(
            sample_index=sample_index,
            iteration_index=sample_index * chain_report.sample_every,
            posterior_log_score=sampled_state.posterior_log_score,
            state=sampled_state,
        )
        for sample_index, sampled_state in enumerate(chain_report.sampled_states)
    ]
    if not sample_rows:
        raise PhylogeneticsError(
            "metropolis-hastings burn-in policy application requires at least one sampled state",
            code="metropolis_hastings_burnin_sample_rows_empty",
        )
    return sample_rows


def _validate_retained_sample_count(
    *,
    total_sample_count: int,
    discarded_sample_count: int,
) -> None:
    if discarded_sample_count >= total_sample_count:
        raise PhylogeneticsError(
            "metropolis-hastings burn-in policy must retain at least one sampled state",
            code="metropolis_hastings_burnin_retained_sample_count_empty",
            details={
                "total_sample_count": total_sample_count,
                "discarded_sample_count": discarded_sample_count,
            },
        )


def _build_diagnostic_series(
    chain_report: MetropolisHastingsRunReport,
) -> dict[str, list[float]]:
    scalar_parameter_names = sorted(
        {
            parameter_name
            for sampled_state in chain_report.sampled_states
            for parameter_name in sampled_state.model_parameters.scalar_parameters
        }
    )
    diagnostic_series = {
        parameter_name: [
            sampled_state.model_parameters.scalar_parameters[parameter_name]
            for sampled_state in chain_report.sampled_states
            if parameter_name in sampled_state.model_parameters.scalar_parameters
        ]
        for parameter_name in scalar_parameter_names
    }
    diagnostic_series["posterior-log-score"] = [
        sampled_state.posterior_log_score
        for sampled_state in chain_report.sampled_states
    ]
    return diagnostic_series


def _validate_nonblank_string(value: str, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be one string",
            code="metropolis_hastings_burnin_policy_field_type_invalid",
            details={"field_name": field_name},
        )
    validated_value = value.strip()
    if not validated_value:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be nonblank",
            code="metropolis_hastings_burnin_policy_field_blank",
            details={"field_name": field_name},
        )
    return validated_value


def _validate_optional_nonnegative_integer(
    *,
    value: int | None,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool):
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be one integer when provided",
            code="metropolis_hastings_burnin_policy_integer_type_invalid",
            details={"field_name": field_name},
        )
    if value < 0:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be nonnegative",
            code="metropolis_hastings_burnin_policy_integer_negative",
            details={"field_name": field_name, "value": value},
        )
    return value


def _validate_optional_positive_integer(
    *,
    value: int | None,
    field_name: str,
) -> int | None:
    if value is None:
        return None
    validated_value = _validate_optional_nonnegative_integer(
        value=value,
        field_name=field_name,
    )
    if validated_value is None or validated_value > 0:
        return validated_value
    raise PhylogeneticsError(
        f"metropolis-hastings burn-in policy requires '{field_name}' to be positive when provided",
        code="metropolis_hastings_burnin_policy_integer_not_positive",
        details={"field_name": field_name, "value": value},
    )


def _validate_optional_fraction(
    *,
    value: float | None,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be one numeric fraction when provided",
            code="metropolis_hastings_burnin_policy_fraction_type_invalid",
            details={"field_name": field_name},
        )
    validated_value = float(value)
    if not 0.0 <= validated_value < 1.0:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to lie in [0, 1)",
            code="metropolis_hastings_burnin_policy_fraction_out_of_range",
            details={"field_name": field_name, "value": validated_value},
        )
    return validated_value


def _validate_optional_positive_float(
    *,
    value: float | None,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be one finite float when provided",
            code="metropolis_hastings_burnin_policy_float_type_invalid",
            details={"field_name": field_name},
        )
    validated_value = float(value)
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"metropolis-hastings burn-in policy requires '{field_name}' to be positive when provided",
            code="metropolis_hastings_burnin_policy_float_not_positive",
            details={"field_name": field_name, "value": validated_value},
        )
    return validated_value
