from __future__ import annotations

import math

from bijux_phylogenetics.comparative._math import (
    student_t_cdf,
    student_t_quantile,
    student_t_two_sided_p_value,
)


def test_student_t_cdf_matches_cauchy_reference_value() -> None:
    assert math.isclose(student_t_cdf(1.0, 1.0), 0.75, abs_tol=1e-12)


def test_student_t_quantile_matches_known_reference_value() -> None:
    assert math.isclose(
        student_t_quantile(0.975, 10.0),
        2.2281388519649385,
        rel_tol=1e-11,
        abs_tol=1e-11,
    )


def test_student_t_two_sided_p_value_inverts_reference_quantile() -> None:
    statistic = student_t_quantile(0.975, 10.0)
    assert math.isclose(
        student_t_two_sided_p_value(statistic, 10.0),
        0.05,
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
