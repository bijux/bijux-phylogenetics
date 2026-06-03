from __future__ import annotations

from dataclasses import dataclass
from statistics import mean, median

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .posterior_predictive_simulation import (
    PosteriorPredictiveAlignmentSimulationReport,
    PosteriorPredictiveContinuousTraitSimulationReport,
    PosteriorPredictiveDiscreteTraitSimulationReport,
    PosteriorPredictiveObservedStatisticRow,
    PosteriorPredictiveReplicateStatisticRow,
)

_POSTERIOR_PREDICTIVE_P_VALUE_METHOD = (
    "posterior-predictive tail area with plus-one Monte Carlo smoothing"
)


@dataclass(frozen=True, slots=True)
class PosteriorPredictivePValueRow:
    """One posterior-predictive tail-area summary for one observed statistic."""

    statistic_name: str
    observed_value: float
    replicate_count: int
    lower_tail_count: int
    upper_tail_count: int
    lower_tail_probability: float
    upper_tail_probability: float
    posterior_predictive_p_value: float
    replicate_mean: float
    replicate_median: float
    replicate_minimum: float
    replicate_maximum: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictivePValueReport:
    """Posterior-predictive p-value summaries across one simulation report."""

    report_kind: str
    model_name: str
    statistic_count: int
    p_value_method: str
    statistic_rows: list[PosteriorPredictivePValueRow]
    warnings: list[str]


def summarize_posterior_predictive_p_values(
    simulation_report: (
        PosteriorPredictiveAlignmentSimulationReport
        | PosteriorPredictiveDiscreteTraitSimulationReport
        | PosteriorPredictiveContinuousTraitSimulationReport
    ),
) -> PosteriorPredictivePValueReport:
    """Compare observed statistics against posterior predictive replicates."""
    report_kind = _resolve_report_kind(simulation_report)
    observed_statistic_rows = _validate_observed_statistic_rows(
        simulation_report.observed_statistic_rows,
        owner_name="posterior predictive p-value summary",
    )
    replicate_values_by_statistic = _group_replicate_values_by_statistic(
        simulation_report.replicate_statistic_rows,
        owner_name="posterior predictive p-value summary",
    )
    statistic_rows: list[PosteriorPredictivePValueRow] = []
    for observed_row in observed_statistic_rows:
        replicate_values = replicate_values_by_statistic.get(
            observed_row.statistic_name
        )
        if replicate_values is None:
            raise PhylogeneticsError(
                "posterior predictive p-value summary requires replicate values for every observed statistic",
                code="posterior_predictive_p_value_replicates_missing",
                details={"statistic_name": observed_row.statistic_name},
            )
        lower_tail_count = sum(
            value <= observed_row.value for value in replicate_values
        )
        upper_tail_count = sum(
            value >= observed_row.value for value in replicate_values
        )
        replicate_count = len(replicate_values)
        lower_tail_probability = _round_float(
            (lower_tail_count + 1) / (replicate_count + 1)
        )
        upper_tail_probability = _round_float(
            (upper_tail_count + 1) / (replicate_count + 1)
        )
        statistic_rows.append(
            PosteriorPredictivePValueRow(
                statistic_name=observed_row.statistic_name,
                observed_value=observed_row.value,
                replicate_count=replicate_count,
                lower_tail_count=lower_tail_count,
                upper_tail_count=upper_tail_count,
                lower_tail_probability=lower_tail_probability,
                upper_tail_probability=upper_tail_probability,
                posterior_predictive_p_value=_round_float(
                    min(1.0, 2.0 * min(lower_tail_probability, upper_tail_probability))
                ),
                replicate_mean=_round_float(mean(replicate_values)),
                replicate_median=_round_float(median(replicate_values)),
                replicate_minimum=_round_float(min(replicate_values)),
                replicate_maximum=_round_float(max(replicate_values)),
            )
        )
    return PosteriorPredictivePValueReport(
        report_kind=report_kind,
        model_name=simulation_report.model_name,
        statistic_count=len(statistic_rows),
        p_value_method=_POSTERIOR_PREDICTIVE_P_VALUE_METHOD,
        statistic_rows=statistic_rows,
        warnings=[],
    )


def _resolve_report_kind(
    simulation_report: object,
) -> str:
    if isinstance(simulation_report, PosteriorPredictiveAlignmentSimulationReport):
        return "alignment"
    if isinstance(simulation_report, PosteriorPredictiveDiscreteTraitSimulationReport):
        return "discrete-trait"
    if isinstance(
        simulation_report, PosteriorPredictiveContinuousTraitSimulationReport
    ):
        return "continuous-trait"
    raise PhylogeneticsError(
        "posterior predictive p-value summary requires one native posterior predictive simulation report",
        code="posterior_predictive_p_value_report_type_invalid",
        details={"report_type": type(simulation_report).__name__},
    )


def _validate_observed_statistic_rows(
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow],
    *,
    owner_name: str,
) -> list[PosteriorPredictiveObservedStatisticRow]:
    statistic_rows_by_name: dict[str, PosteriorPredictiveObservedStatisticRow] = {}
    for row in observed_statistic_rows:
        if row.statistic_name in statistic_rows_by_name:
            raise PhylogeneticsError(
                f"{owner_name} requires each observed statistic name to appear once",
                code="posterior_predictive_p_value_observed_statistic_duplicated",
                details={"statistic_name": row.statistic_name},
            )
        statistic_rows_by_name[row.statistic_name] = row
    if not statistic_rows_by_name:
        raise PhylogeneticsError(
            f"{owner_name} requires at least one observed statistic row",
            code="posterior_predictive_p_value_observed_statistics_empty",
        )
    return list(statistic_rows_by_name.values())


def _group_replicate_values_by_statistic(
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow],
    *,
    owner_name: str,
) -> dict[str, list[float]]:
    if not replicate_statistic_rows:
        raise PhylogeneticsError(
            f"{owner_name} requires at least one replicate statistic row",
            code="posterior_predictive_p_value_replicate_statistics_empty",
        )
    values_by_statistic: dict[str, list[float]] = {}
    for row in replicate_statistic_rows:
        values_by_statistic.setdefault(row.statistic_name, []).append(row.value)
    return values_by_statistic


def _round_float(value: float) -> float:
    return float(format(value, ".15g"))


__all__ = [
    "PosteriorPredictivePValueReport",
    "PosteriorPredictivePValueRow",
    "summarize_posterior_predictive_p_values",
]
