from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.discrete_evolution import (
    compare_discrete_state_models,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    render_discrete_state_evolution_report,
    render_tree_with_geographic_states,
    run_discrete_state_transition_model,
    validate_discrete_state_coding,
    write_node_state_probability_table,
)


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_validate_discrete_state_coding_reports_unsupported_and_delimited_states() -> None:
    report = validate_discrete_state_coding(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography_invalid.tsv"),
        trait="region",
        allowed_states=["north", "south", "island"],
    )

    assert report.valid is False
    assert report.observed_states == ["island", "north"]
    assert {issue.code for issue in report.issues} == {"unsupported-state-delimiter", "unsupported-state-label"}


def test_detect_state_imbalance_problems_flags_single_state_and_dominance() -> None:
    report = detect_state_imbalance_problems(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography_single_state.tsv"),
        trait="region",
    )

    assert report.state_counts == {"north": 4}
    assert {warning.code for warning in report.warnings} == {"single-state-dataset", "dominant-state-skew"}


def test_run_discrete_state_transition_model_returns_transition_matrix_and_events() -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )

    assert report.transition_model.model == "equal-rates"
    assert report.transition_model.state_order == ["island", "north", "south"]
    assert report.transition_summary.branch_count == 6
    assert report.transition_summary.transition_count >= 1
    assert all(sum(row.target_rates.values()) > 0.99 for row in report.transition_model.transition_matrix)


def test_estimate_ancestral_geographic_states_and_compare_models_return_node_differences() -> None:
    report = estimate_ancestral_geographic_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="all-rates-different",
    )
    comparison = compare_discrete_state_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
    )

    assert report.estimates[0].state_probabilities
    assert comparison.better_model in {"equal-rates", "all-rates-different"}
    assert len(comparison.rows) == 2
    assert len(comparison.node_differences) == len(report.estimates)


def test_write_node_probability_table_and_render_report_outputs_files(tmp_path: Path) -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    table_path = tmp_path / "node-probabilities.tsv"
    svg_path = tmp_path / "geography.svg"
    html_path = tmp_path / "discrete-report.html"

    write_node_state_probability_table(table_path, report)
    render_result = render_tree_with_geographic_states(
        fixture("example_tree.nwk"),
        report,
        out_path=svg_path,
    )
    report_result = render_discrete_state_evolution_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_geography.tsv"),
        trait="region",
        out_path=html_path,
        model="equal-rates",
        compare_model="all-rates-different",
    )

    assert "state_probabilities" in table_path.read_text(encoding="utf-8")
    assert render_result.rendered_internal_annotation_count == 3
    assert "discrete-state-evolution" in html_path.read_text(encoding="utf-8")
    assert report_result.report_kind == "discrete-state-evolution"
