from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    build_invariant_proportion_boundary_warnings,
    build_rate_heterogeneity_alpha_boundary_warnings,
)


def test_invariant_boundary_warning_builder_marks_lower_and_upper_hits() -> None:
    lower = build_invariant_proportion_boundary_warnings(
        invariant_proportion=0.0,
        lower_invariant_proportion_bound=0.0,
        upper_invariant_proportion_bound=0.95,
    )
    upper = build_invariant_proportion_boundary_warnings(
        invariant_proportion=0.95,
        lower_invariant_proportion_bound=0.0,
        upper_invariant_proportion_bound=0.95,
    )

    assert [warning.warning_kind for warning in lower] == ["invariant-boundary"]
    assert [warning.boundary_side for warning in lower] == ["lower"]
    assert [warning.message for warning in lower] == [
        "invariant proportion hit lower search boundary"
    ]
    assert [warning.warning_kind for warning in upper] == ["invariant-boundary"]
    assert [warning.boundary_side for warning in upper] == ["upper"]
    assert [warning.message for warning in upper] == [
        "invariant proportion hit upper search boundary"
    ]


def test_alpha_boundary_warning_builder_distinguishes_flat_upper_surface() -> None:
    lower = build_rate_heterogeneity_alpha_boundary_warnings(
        alpha=0.05,
        lower_alpha_bound=0.05,
        upper_alpha_bound=20.0,
    )
    upper = build_rate_heterogeneity_alpha_boundary_warnings(
        alpha=20.0,
        lower_alpha_bound=0.05,
        upper_alpha_bound=20.0,
    )

    assert [warning.warning_kind for warning in lower] == ["alpha-boundary"]
    assert [warning.boundary_side for warning in lower] == ["lower"]
    assert [warning.message for warning in lower] == ["alpha hit lower search boundary"]
    assert [warning.warning_kind for warning in upper] == ["flat-alpha"]
    assert [warning.boundary_side for warning in upper] == ["upper"]
    assert [warning.message for warning in upper] == ["alpha hit upper search boundary"]
