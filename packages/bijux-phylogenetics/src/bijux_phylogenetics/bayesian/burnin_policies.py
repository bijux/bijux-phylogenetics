from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
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
    if policy.policy_name == "none":
        discarded_sample_count = 0
    elif policy.policy_name == "fixed-count":
        discarded_sample_count = policy.discarded_sample_count or 0
    elif policy.policy_name == "fixed-fraction":
        discarded_sample_count = int(
            len(sample_rows) * (policy.discarded_fraction or 0.0)
        )
    else:
        raise PhylogeneticsError(
            "diagnostic-suggested burn-in policy requires diagnostic suggestion support",
            code="metropolis_hastings_burnin_policy_diagnostic_support_missing",
        )
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


def _build_sample_rows(chain_report: MetropolisHastingsRunReport) -> list[BurninSampleRow]:
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
