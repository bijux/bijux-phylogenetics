from __future__ import annotations

import math

from .contracts import ComparativeModelComparisonRow


def compute_aic(log_likelihood: float, *, parameter_count: int) -> float:
    """Compute AIC from one maximized log-likelihood surface."""
    return (2.0 * parameter_count) - (2.0 * log_likelihood)


def compute_aicc(aic: float, *, sample_size: int, parameter_count: int) -> float:
    """Compute AICc, or infinity when the correction is undefined."""
    denominator = sample_size - parameter_count - 1
    if denominator <= 0:
        return math.inf
    return aic + ((2.0 * parameter_count * (parameter_count + 1)) / denominator)


def rank_model_comparison_rows(
    rows: list[ComparativeModelComparisonRow],
    *,
    delta_threshold: float,
) -> tuple[str | None, list[str]]:
    """Rank comparable model rows by AICc and annotate Akaike-weight support."""
    comparable_rows = [
        row
        for row in rows
        if row.comparable and math.isfinite(row.aic) and math.isfinite(row.aicc)
    ]
    if not comparable_rows:
        return None, []
    shared_policy = comparable_rows[0].likelihood_constant_policy
    if shared_policy is None:
        blocked_models = [row.model for row in comparable_rows]
        for row in comparable_rows:
            row.comparable = False
            row.comparability_note = "likelihood constant policy is missing, so AIC and AICc ranking is blocked"
        _clear_unranked_rows(rows)
        return None, blocked_models
    policies = {row.likelihood_constant_policy for row in comparable_rows}
    if len(policies) > 1:
        blocked_models = [row.model for row in comparable_rows]
        for row in comparable_rows:
            row.comparable = False
            row.comparability_note = (
                "likelihood constant policy differs across candidate models, so AIC and "
                "AICc ranking is blocked for the full comparison surface"
            )
        _clear_unranked_rows(rows)
        return None, blocked_models

    best_aic = min(row.aic for row in comparable_rows)
    best_aicc = min(row.aicc for row in comparable_rows)
    ranked_rows = sorted(
        comparable_rows,
        key=lambda row: (row.aicc, row.aic, row.model),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row.rank = rank
        row.delta_aic = row.aic - best_aic
        row.delta_aicc = row.aicc - best_aicc
        row.selected = math.isclose(
            row.aicc,
            best_aicc,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    raw_weights = [math.exp(-0.5 * row.delta_aicc) for row in ranked_rows]
    weight_total = sum(raw_weights)
    for row, raw_weight in zip(ranked_rows, raw_weights, strict=True):
        row.akaike_weight = raw_weight / weight_total if weight_total else 0.0
        row.within_delta_aic_threshold = row.delta_aic <= delta_threshold
        row.within_delta_aicc_threshold = row.delta_aicc <= delta_threshold
    _clear_unranked_rows([row for row in rows if row not in comparable_rows])
    rows.sort(
        key=lambda row: (
            row.rank is None,
            math.inf if row.rank is None else row.rank,
            row.model,
        )
    )
    return shared_policy, []


def _clear_unranked_rows(rows: list[ComparativeModelComparisonRow]) -> None:
    for row in rows:
        row.rank = None
        row.delta_aic = math.inf
        row.delta_aicc = math.inf
        row.selected = False
        row.akaike_weight = None
        row.within_delta_aic_threshold = None
        row.within_delta_aicc_threshold = None
