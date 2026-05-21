from __future__ import annotations

import math

from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
    compute_aic,
    compute_aicc,
    rank_model_comparison_rows,
)


def test_compute_aicc_returns_infinity_when_correction_is_undefined() -> None:
    aic = compute_aic(-10.0, parameter_count=4)

    assert math.isinf(compute_aicc(aic, sample_size=5, parameter_count=4))


def test_rank_model_comparison_rows_blocks_incompatible_likelihood_policies() -> None:
    rows = [
        ComparativeModelComparisonRow(
            model="brownian",
            parameter_count=2,
            log_likelihood=-5.0,
            aic=14.0,
            aicc=15.0,
            likelihood_constant_policy="policy-a",
        ),
        ComparativeModelComparisonRow(
            model="white-noise",
            parameter_count=2,
            log_likelihood=-4.0,
            aic=12.0,
            aicc=13.0,
            likelihood_constant_policy="policy-b",
        ),
    ]

    shared_policy, blocked_models = rank_model_comparison_rows(
        rows,
        delta_threshold=2.0,
    )

    assert shared_policy is None
    assert blocked_models == ["brownian", "white-noise"]
    assert all(row.comparable is False for row in rows)
    assert all(row.rank is None for row in rows)
