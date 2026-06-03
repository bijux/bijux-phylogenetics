from __future__ import annotations

import pytest

import bijux_phylogenetics.simulation as simulation_api
from bijux_phylogenetics.simulation import (
    CoalescentSkylineSummaryRow,
    simulate_coalescent_tree,
    simulate_coalescent_trees,
    write_coalescent_skyline_table,
)


def test_simulation_gateway_exports_coalescent_skyline_surface() -> None:
    assert simulation_api.CoalescentSkylineSummaryRow is CoalescentSkylineSummaryRow
    assert (
        simulation_api.write_coalescent_skyline_table is write_coalescent_skyline_table
    )


def test_simulate_coalescent_trees_reports_piecewise_effective_population_estimates() -> (
    None
):
    _trees, report = simulate_coalescent_trees(
        tree_count=256,
        tip_count=6,
        population_size=3.0,
        waiting_time_tolerance=0.15,
        seed=17,
    )

    assert [row.interval for row in report.coalescent_skyline_rows] == [
        "6->5",
        "5->4",
        "4->3",
        "3->2",
        "2->1",
    ]
    assert [row.lineage_count for row in report.coalescent_skyline_rows] == [
        6,
        5,
        4,
        3,
        2,
    ]
    assert [
        row.effective_population_size_estimate for row in report.coalescent_skyline_rows
    ] == pytest.approx(
        [
            2.87420076024687,
            2.85842047438072,
            3.40477055758602,
            2.970041631191979,
            3.071212929051539,
        ]
    )
    assert all(row.observation_count == 256 for row in report.coalescent_skyline_rows)
    assert {row.uncertainty_flag for row in report.coalescent_skyline_rows} == {"low"}


def test_simulate_coalescent_tree_preserves_skyline_rows_on_single_tree_surface() -> (
    None
):
    _tree, report = simulate_coalescent_tree(
        tip_count=4,
        population_size=2.5,
        seed=13,
    )
    _trees, batch_report = simulate_coalescent_trees(
        tree_count=1,
        tip_count=4,
        population_size=2.5,
        seed=13,
    )

    assert report.coalescent_skyline_rows == batch_report.coalescent_skyline_rows
