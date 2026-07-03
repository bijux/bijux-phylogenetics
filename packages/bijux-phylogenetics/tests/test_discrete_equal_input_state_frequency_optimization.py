from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    optimize_discrete_equal_input_state_frequencies,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_equal_input_state_frequency_optimization_keeps_simplex_valid() -> (
    None
):
    report = optimize_discrete_equal_input_state_frequencies(
        load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk")),
        [
            AlignmentRecord(identifier="A", sequence="0"),
            AlignmentRecord(identifier="B", sequence="2"),
        ],
        state_order=("0", "1", "2", "3"),
        initial_state_frequencies={
            "0": 0.25,
            "1": 0.25,
            "2": 0.25,
            "3": 0.25,
        },
        model_name="finite-state equal-input",
        max_coordinate_passes=4,
    )

    optimized_values = {
        row.state_label: row.optimized_value for row in report.state_rows
    }

    assert report.model_name == "finite-state equal-input"
    assert report.state_count == 4
    assert report.state_frequency_source == "optimized"
    assert report.fit_assessment.classification in {"failed", "weak"}
    assert report.function_evaluation_count > 1
    assert report.optimization_pass_count >= 1
    assert report.optimized_log_likelihood >= report.initial_log_likelihood
    assert math.isclose(
        sum(optimized_values.values()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert min(optimized_values.values()) > 0.0
    assert optimized_values["0"] > optimized_values["1"]
    assert optimized_values["2"] > optimized_values["3"]
    if not report.converged:
        assert (
            "finite-state equal-input stationary-frequency search did not converge"
            in report.warnings
        )
    assert report.parameter_diagnostics


def test_discrete_equal_input_state_frequency_optimization_emits_structured_floor_warnings() -> (
    None
):
    report = optimize_discrete_equal_input_state_frequencies(
        load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk")),
        [
            AlignmentRecord(identifier="A", sequence="0"),
            AlignmentRecord(identifier="B", sequence="0"),
        ],
        state_order=("0", "1", "2", "3"),
        initial_state_frequencies={
            "0": 0.25,
            "1": 0.25,
            "2": 0.25,
            "3": 0.25,
        },
        model_name="finite-state equal-input",
        minimum_state_frequency=0.2,
        max_coordinate_passes=6,
    )

    assert report.boundary_warning_rows
    assert report.fit_assessment.classification == "overparameterized"
    assert all(
        warning.warning_kind == "frequency-boundary"
        for warning in report.boundary_warning_rows
    )
    assert all(
        warning.boundary_side == "lower" for warning in report.boundary_warning_rows
    )
    assert all(
        warning.message in report.warnings for warning in report.boundary_warning_rows
    )
    assert {warning.affected_parameter for warning in report.boundary_warning_rows} == {
        "state_frequency[1]",
        "state_frequency[2]",
        "state_frequency[3]",
    }
