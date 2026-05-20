from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

DEFAULT_BURNIN_FRACTIONS: tuple[float, ...] = (0.05, 0.1, 0.25, 0.5)
MAJORITY_CLADE_THRESHOLD = 0.5


@dataclass(frozen=True, slots=True)
class BurninParameterEstimate:
    burnin_fraction: float
    mean: float
    median: float
    hpd_95_lower: float
    hpd_95_upper: float
    effective_sample_size: float


@dataclass(slots=True)
class BurninSensitivityParameterShift:
    parameter: str
    estimates: list[BurninParameterEstimate]
    minimum_mean: float
    maximum_mean: float
    mean_range: float
    minimum_median: float
    maximum_median: float
    median_range: float
    common_hpd_95_lower: float | None
    common_hpd_95_upper: float | None
    minimum_effective_sample_size: float
    maximum_effective_sample_size: float
    unstable: bool
    instability_reasons: list[str]


@dataclass(frozen=True, slots=True)
class BurninCladeProbabilityEstimate:
    burnin_fraction: float
    posterior_probability: float


@dataclass(slots=True)
class BurninSensitivityCladeShift:
    clade: str
    estimates: list[BurninCladeProbabilityEstimate]
    minimum_posterior_probability: float
    maximum_posterior_probability: float
    posterior_probability_range: float
    crosses_majority_threshold: bool
    unstable: bool
    instability_reasons: list[str]


def normalize_burnin_fractions(
    burnin_fractions: tuple[float, ...] | list[float],
) -> tuple[float, ...]:
    """Validate, deduplicate, and sort requested burn-in fractions."""
    if not burnin_fractions:
        raise ValueError("burnin_fractions must contain at least one value")
    normalized = tuple(
        sorted(dict.fromkeys(float(value) for value in burnin_fractions))
    )
    for value in normalized:
        if not 0.0 <= value < 1.0:
            raise ValueError(f"burnin_fraction must be between 0 and 1, got {value}")
    return normalized


def summarize_burnin_parameter_shifts(
    parameter_summaries_by_fraction: dict[float, list[Any]],
) -> list[BurninSensitivityParameterShift]:
    """Compare posterior parameter summaries across burn-in fractions."""
    ordered_fractions = sorted(parameter_summaries_by_fraction)
    parameter_names = sorted(
        {
            str(summary.parameter)
            for summaries in parameter_summaries_by_fraction.values()
            for summary in summaries
        }
    )
    shifts: list[BurninSensitivityParameterShift] = []
    for parameter in parameter_names:
        estimates = [
            BurninParameterEstimate(
                burnin_fraction=fraction,
                mean=float(summary.mean),
                median=float(summary.median),
                hpd_95_lower=float(summary.hpd_95_lower),
                hpd_95_upper=float(summary.hpd_95_upper),
                effective_sample_size=float(summary.effective_sample_size),
            )
            for fraction in ordered_fractions
            for summary in parameter_summaries_by_fraction[fraction]
            if str(summary.parameter) == parameter
        ]
        if not estimates:
            continue
        hpd_lower = max(estimate.hpd_95_lower for estimate in estimates)
        hpd_upper = min(estimate.hpd_95_upper for estimate in estimates)
        instability_reasons: list[str] = []
        if hpd_lower > hpd_upper:
            instability_reasons.append(
                "95% HPD intervals do not share a common overlap"
            )
        shifts.append(
            BurninSensitivityParameterShift(
                parameter=parameter,
                estimates=estimates,
                minimum_mean=min(estimate.mean for estimate in estimates),
                maximum_mean=max(estimate.mean for estimate in estimates),
                mean_range=round(
                    max(estimate.mean for estimate in estimates)
                    - min(estimate.mean for estimate in estimates),
                    6,
                ),
                minimum_median=min(estimate.median for estimate in estimates),
                maximum_median=max(estimate.median for estimate in estimates),
                median_range=round(
                    max(estimate.median for estimate in estimates)
                    - min(estimate.median for estimate in estimates),
                    6,
                ),
                common_hpd_95_lower=None
                if hpd_lower > hpd_upper
                else round(hpd_lower, 6),
                common_hpd_95_upper=None
                if hpd_lower > hpd_upper
                else round(hpd_upper, 6),
                minimum_effective_sample_size=min(
                    estimate.effective_sample_size for estimate in estimates
                ),
                maximum_effective_sample_size=max(
                    estimate.effective_sample_size for estimate in estimates
                ),
                unstable=bool(instability_reasons),
                instability_reasons=instability_reasons,
            )
        )
    return shifts


def summarize_burnin_clade_shifts(
    clade_frequencies_by_fraction: dict[float, list[Any]],
) -> list[BurninSensitivityCladeShift]:
    """Compare posterior clade probabilities across burn-in fractions."""
    ordered_fractions = sorted(clade_frequencies_by_fraction)
    clade_names = sorted(
        {
            str(row.clade)
            for rows in clade_frequencies_by_fraction.values()
            for row in rows
        }
    )
    shifts: list[BurninSensitivityCladeShift] = []
    for clade in clade_names:
        estimates = [
            BurninCladeProbabilityEstimate(
                burnin_fraction=fraction,
                posterior_probability=float(
                    next(
                        (
                            row.frequency
                            for row in clade_frequencies_by_fraction[fraction]
                            if str(row.clade) == clade
                        ),
                        0.0,
                    )
                ),
            )
            for fraction in ordered_fractions
        ]
        minimum_probability = min(
            estimate.posterior_probability for estimate in estimates
        )
        maximum_probability = max(
            estimate.posterior_probability for estimate in estimates
        )
        crosses_majority_threshold = (
            minimum_probability < MAJORITY_CLADE_THRESHOLD <= maximum_probability
        )
        instability_reasons: list[str] = []
        if crosses_majority_threshold:
            instability_reasons.append(
                "posterior clade probability crosses the majority-rule threshold"
            )
        shifts.append(
            BurninSensitivityCladeShift(
                clade=clade,
                estimates=estimates,
                minimum_posterior_probability=minimum_probability,
                maximum_posterior_probability=maximum_probability,
                posterior_probability_range=round(
                    maximum_probability - minimum_probability, 6
                ),
                crosses_majority_threshold=crosses_majority_threshold,
                unstable=bool(instability_reasons),
                instability_reasons=instability_reasons,
            )
        )
    return shifts


def write_burnin_parameter_shift_table(
    path: Path,
    shifts: list[BurninSensitivityParameterShift],
) -> Path:
    """Write one row per posterior parameter burn-in comparison."""
    return write_taxon_rows(
        path,
        columns=[
            "parameter",
            "minimum_mean",
            "maximum_mean",
            "mean_range",
            "minimum_median",
            "maximum_median",
            "median_range",
            "common_hpd_95_lower",
            "common_hpd_95_upper",
            "minimum_effective_sample_size",
            "maximum_effective_sample_size",
            "unstable",
            "instability_reasons",
            "median_by_fraction",
            "hpd_95_by_fraction",
        ],
        rows=[
            {
                "parameter": shift.parameter,
                "minimum_mean": format(shift.minimum_mean, ".15g"),
                "maximum_mean": format(shift.maximum_mean, ".15g"),
                "mean_range": format(shift.mean_range, ".15g"),
                "minimum_median": format(shift.minimum_median, ".15g"),
                "maximum_median": format(shift.maximum_median, ".15g"),
                "median_range": format(shift.median_range, ".15g"),
                "common_hpd_95_lower": ""
                if shift.common_hpd_95_lower is None
                else format(shift.common_hpd_95_lower, ".15g"),
                "common_hpd_95_upper": ""
                if shift.common_hpd_95_upper is None
                else format(shift.common_hpd_95_upper, ".15g"),
                "minimum_effective_sample_size": format(
                    shift.minimum_effective_sample_size, ".15g"
                ),
                "maximum_effective_sample_size": format(
                    shift.maximum_effective_sample_size, ".15g"
                ),
                "unstable": str(shift.unstable).lower(),
                "instability_reasons": "; ".join(shift.instability_reasons),
                "median_by_fraction": "; ".join(
                    f"{format(estimate.burnin_fraction, '.15g')}={format(estimate.median, '.15g')}"
                    for estimate in shift.estimates
                ),
                "hpd_95_by_fraction": "; ".join(
                    f"{format(estimate.burnin_fraction, '.15g')}=[{format(estimate.hpd_95_lower, '.15g')}, {format(estimate.hpd_95_upper, '.15g')}]"
                    for estimate in shift.estimates
                ),
            }
            for shift in shifts
        ],
    )


def write_burnin_clade_shift_table(
    path: Path,
    shifts: list[BurninSensitivityCladeShift],
) -> Path:
    """Write one row per posterior clade burn-in comparison."""
    return write_taxon_rows(
        path,
        columns=[
            "clade",
            "minimum_posterior_probability",
            "maximum_posterior_probability",
            "posterior_probability_range",
            "crosses_majority_threshold",
            "unstable",
            "instability_reasons",
            "probability_by_fraction",
        ],
        rows=[
            {
                "clade": shift.clade,
                "minimum_posterior_probability": format(
                    shift.minimum_posterior_probability, ".15g"
                ),
                "maximum_posterior_probability": format(
                    shift.maximum_posterior_probability, ".15g"
                ),
                "posterior_probability_range": format(
                    shift.posterior_probability_range, ".15g"
                ),
                "crosses_majority_threshold": str(
                    shift.crosses_majority_threshold
                ).lower(),
                "unstable": str(shift.unstable).lower(),
                "instability_reasons": "; ".join(shift.instability_reasons),
                "probability_by_fraction": "; ".join(
                    f"{format(estimate.burnin_fraction, '.15g')}={format(estimate.posterior_probability, '.15g')}"
                    for estimate in shift.estimates
                ),
            }
            for shift in shifts
        ],
    )
