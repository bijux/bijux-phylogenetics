from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import TYPE_CHECKING

from bijux_phylogenetics.comparative.evolutionary_modes import (
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
)

from .registry import GeigerParityCase, list_geiger_parity_cases

if TYPE_CHECKING:
    from .runner import GeigerParityObservation, GeigerParityReport


@dataclass(frozen=True, slots=True)
class GeigerLikelihoodPolicyRow:
    """One governed record of a geiger-versus-Bijux likelihood-constant policy."""

    case_id: str
    function_name: str
    model_name: str
    status: str
    reference_likelihood_constant_policy: str
    bijux_likelihood_constant_policy: str
    case_level_raw_log_likelihood_comparable: bool
    raw_log_likelihood_match_within_tolerance: bool | None
    reference_aic_matches_raw_log_likelihood: bool | None
    bijux_aic_matches_raw_log_likelihood: bool | None
    relative_aic_comparable: bool
    ranking_permitted: bool
    ranking_guard_outcome: str
    policy_evidence: str


def build_geiger_likelihood_policy_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerLikelihoodPolicyRow]:
    """Build the governed geiger likelihood-constant policy rows."""

    case_by_id = {case.case_id: case for case in list_geiger_parity_cases()}
    return [
        _build_likelihood_policy_row(observation, case_by_id[observation.case_id])
        for observation in observations
    ]


def write_geiger_likelihood_policy_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write the governed geiger likelihood-constant policy table as TSV."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "model_name",
                "status",
                "reference_likelihood_constant_policy",
                "bijux_likelihood_constant_policy",
                "case_level_raw_log_likelihood_comparable",
                "raw_log_likelihood_match_within_tolerance",
                "reference_aic_matches_raw_log_likelihood",
                "bijux_aic_matches_raw_log_likelihood",
                "relative_aic_comparable",
                "ranking_permitted",
                "ranking_guard_outcome",
                "policy_evidence",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.likelihood_policy_rows:
            writer.writerow(asdict(row))
    return path


def _build_likelihood_policy_row(
    observation: GeigerParityObservation,
    case: GeigerParityCase,
) -> GeigerLikelihoodPolicyRow:
    reference_policy = _likelihood_policy_value(
        observation.reference_summary,
        CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    )
    bijux_policy = _likelihood_policy_value(
        observation.bijux_summary,
        CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    )
    if case.operation not in {"fit-continuous", "compare-fitcontinuous-models"}:
        return GeigerLikelihoodPolicyRow(
            case_id=observation.case_id,
            function_name=observation.function_name,
            model_name=observation.model_name,
            status=observation.status,
            reference_likelihood_constant_policy=reference_policy,
            bijux_likelihood_constant_policy=bijux_policy,
            case_level_raw_log_likelihood_comparable=False,
            raw_log_likelihood_match_within_tolerance=None,
            reference_aic_matches_raw_log_likelihood=None,
            bijux_aic_matches_raw_log_likelihood=None,
            relative_aic_comparable=False,
            ranking_permitted=False,
            ranking_guard_outcome="outside-current-fitcontinuous-likelihood-policy-scope",
            policy_evidence=(
                "This governed likelihood-constant registry tranche only covers the "
                "owned fitContinuous parity surface; discrete fitDiscrete lanes remain "
                "outside the current analytical likelihood-policy scope."
            ),
        )
    if case.operation == "compare-fitcontinuous-models":
        reference_aic_matches = _all_model_comparison_rows_match_aic_formula(
            observation.reference_rows
        )
        bijux_aic_matches = _all_model_comparison_rows_match_aic_formula(
            observation.bijux_rows
        )
        ranking_permitted = _ranking_permitted(observation.bijux_summary)
        ranking_guard_outcome = (
            "shared-fitcontinuous-policy-ranking-permitted"
            if ranking_permitted
            else "ranking-blocked-by-likelihood-policy-guard"
        )
        return GeigerLikelihoodPolicyRow(
            case_id=observation.case_id,
            function_name=observation.function_name,
            model_name=observation.model_name,
            status=observation.status,
            reference_likelihood_constant_policy=reference_policy,
            bijux_likelihood_constant_policy=bijux_policy,
            case_level_raw_log_likelihood_comparable=False,
            raw_log_likelihood_match_within_tolerance=None,
            reference_aic_matches_raw_log_likelihood=reference_aic_matches,
            bijux_aic_matches_raw_log_likelihood=bijux_aic_matches,
            relative_aic_comparable=True,
            ranking_permitted=ranking_permitted,
            ranking_guard_outcome=ranking_guard_outcome,
            policy_evidence=(
                "Model-comparison cases do not expose one case-level log likelihood; "
                "instead, each candidate row must satisfy AIC = 2k - 2logLik and the "
                "owned fitcontinuous ranking guard must keep all ranked models on one "
                f"shared policy: {FITCONTINUOUS_MODEL_RANKING_POLICY}"
            ),
        )
    return GeigerLikelihoodPolicyRow(
        case_id=observation.case_id,
        function_name=observation.function_name,
        model_name=observation.model_name,
        status=observation.status,
        reference_likelihood_constant_policy=reference_policy,
        bijux_likelihood_constant_policy=bijux_policy,
        case_level_raw_log_likelihood_comparable=True,
        raw_log_likelihood_match_within_tolerance=_same_numeric_value(
            _optional_float(
                None
                if observation.reference_summary is None
                else observation.reference_summary.get("log_likelihood")
            ),
            _optional_float(
                None
                if observation.bijux_summary is None
                else observation.bijux_summary.get("log_likelihood")
            ),
            tolerance=observation.tolerance,
        ),
        reference_aic_matches_raw_log_likelihood=_aic_matches_summary_log_likelihood(
            observation.reference_summary,
            case=case,
            tolerance=observation.tolerance,
        ),
        bijux_aic_matches_raw_log_likelihood=_aic_matches_summary_log_likelihood(
            observation.bijux_summary,
            case=case,
            tolerance=observation.tolerance,
        ),
        relative_aic_comparable=True,
        ranking_permitted=False,
        ranking_guard_outcome="single-fit-case-no-cross-case-ranking",
        policy_evidence=(
            "Single fitContinuous cases keep the owned and reference summaries on the "
            f"same direct policy: {CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY}"
        ),
    )


def _likelihood_policy_value(
    summary: dict[str, object] | None,
    default: str,
) -> str:
    if summary is None:
        return default
    value = summary.get("likelihood_constant_policy")
    return value if isinstance(value, str) and value else default


def _aic_matches_summary_log_likelihood(
    summary: dict[str, object] | None,
    *,
    case: GeigerParityCase,
    tolerance: float,
) -> bool | None:
    if summary is None:
        return None
    aic = _optional_float(summary.get("aic"))
    log_likelihood = _optional_float(summary.get("log_likelihood"))
    if aic is None or log_likelihood is None:
        return None
    expected_aic = (2.0 * _parameter_count(case)) - (2.0 * log_likelihood)
    return math.isclose(aic, expected_aic, rel_tol=tolerance, abs_tol=tolerance)


def _all_model_comparison_rows_match_aic_formula(
    rows: list[dict[str, object]] | None,
) -> bool | None:
    if rows is None:
        return None
    comparable_rows = [row for row in rows if row.get("comparable") is True]
    if not comparable_rows:
        return None
    return all(
        _model_comparison_row_matches_aic_formula(row) for row in comparable_rows
    )


def _model_comparison_row_matches_aic_formula(row: dict[str, object]) -> bool:
    parameter_count = _optional_float(row.get("parameter_count"))
    log_likelihood = _optional_float(row.get("log_likelihood"))
    aic = _optional_float(row.get("aic"))
    if parameter_count is None or log_likelihood is None or aic is None:
        return False
    expected_aic = (2.0 * parameter_count) - (2.0 * log_likelihood)
    return math.isclose(aic, expected_aic, rel_tol=1e-9, abs_tol=1e-9)


def _ranking_permitted(summary: dict[str, object] | None) -> bool:
    if summary is None:
        return False
    blocked = summary.get("noncomparable_likelihood_models")
    if isinstance(blocked, list):
        return len(blocked) == 0
    return True


def _parameter_count(case: GeigerParityCase) -> int:
    if case.operation == "fit-continuous":
        return 2 if case.model_name in {"BM", "white"} else 3
    raise ValueError(
        "likelihood policy summary AIC checks only apply to single fitContinuous cases"
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        numeric = float(value)
        return None if math.isnan(numeric) else numeric
    return None


def _same_numeric_value(
    left: float | None,
    right: float | None,
    *,
    tolerance: float,
) -> bool | None:
    if left is None or right is None:
        return None
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)
