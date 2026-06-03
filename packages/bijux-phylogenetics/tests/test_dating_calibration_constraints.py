from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.dating as dating_api
from bijux_phylogenetics.phylo.dating.constraints import (
    require_feasible_dating_calibration_constraints,
    solve_dating_calibration_constraints,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_dating_gateway_exports_constraint_solver_surface() -> None:
    assert (
        dating_api.solve_dating_calibration_constraints
        is solve_dating_calibration_constraints
    )


def test_solve_dating_calibration_constraints_maps_bounds_and_reports_feasible() -> (
    None
):
    report = solve_dating_calibration_constraints(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture("metadata", "dating_calibration_constraints_5_taxa.tsv"),
    )

    assert report.tip_count == 5
    assert report.internal_node_count == 4
    assert report.calibration_count == 4
    assert report.valid_calibration_count == 4
    assert report.invalid_calibration_count == 0
    assert report.resolved_calibration_count == 4
    assert report.contradictory_calibration_count == 0
    assert report.contradictory_node_count == 0
    assert report.feasible is True
    rows_by_calibration_id = {row.calibration_id: row for row in report.constraint_rows}
    assert rows_by_calibration_id["cal-ab"].fixed_date == pytest.approx(
        1996.0, abs=1e-12
    )
    assert rows_by_calibration_id["cal-root"].node_kind == "root"
    window_rows_by_taxa = {
        tuple(row.descendant_taxa): row for row in report.node_window_rows
    }
    assert window_rows_by_taxa[("A", "B", "C")].minimum_bound == pytest.approx(
        1988.0,
        abs=1e-12,
    )
    assert window_rows_by_taxa[("D", "E")].maximum_bound == pytest.approx(
        1997.0,
        abs=1e-12,
    )
    assert report.issue_rows == []


def test_solve_dating_calibration_constraints_reports_ancestor_descendant_conflicts() -> (
    None
):
    report = solve_dating_calibration_constraints(
        fixture(
            "trees",
            "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
        ),
        fixture(
            "metadata",
            "dating_calibration_constraints_contradictory_5_taxa.tsv",
        ),
    )

    assert report.feasible is False
    assert report.contradictory_calibration_count == 3
    assert report.contradictory_node_count == 3
    assert [row.code for row in report.issue_rows] == [
        "chronology-conflict",
        "chronology-conflict",
    ]
    window_rows_by_taxa = {
        tuple(row.descendant_taxa): row for row in report.node_window_rows
    }
    assert window_rows_by_taxa[("A", "B", "C")].contradictory is True
    assert window_rows_by_taxa[("A", "B")].contradictory is True


def test_require_feasible_dating_calibration_constraints_raises_for_contradictory_fixture() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="dating calibrations are infeasible",
    ):
        require_feasible_dating_calibration_constraints(
            fixture(
                "trees",
                "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
            ),
            fixture(
                "metadata",
                "dating_calibration_constraints_contradictory_5_taxa.tsv",
            ),
        )
