from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    CladeCompatibilityEdgeRow,
    CladeCompatibilityGraphReport,
    CladeCompatibilityNodeRow,
    compute_clade_compatibility_graph,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_tree_gateway_exports_clade_compatibility_surface() -> None:
    assert trees_api.CladeCompatibilityNodeRow is CladeCompatibilityNodeRow
    assert trees_api.CladeCompatibilityEdgeRow is CladeCompatibilityEdgeRow
    assert trees_api.CladeCompatibilityGraphReport is CladeCompatibilityGraphReport
    assert (
        trees_api.compute_clade_compatibility_graph is compute_clade_compatibility_graph
    )


def test_compute_clade_compatibility_graph_matches_exact_fixture() -> None:
    report = compute_clade_compatibility_graph(
        fixture("clade_compatibility_tree_set.nwk")
    )

    assert report.tree_count == 2
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.node_count == 4
    assert report.edge_count == 6
    assert report.compatible_edge_count == 2
    assert report.conflict_edge_count == 4
    assert report.nodes == [
        CladeCompatibilityNodeRow(
            clade="A|B",
            tree_count=1,
            frequency=0.5,
            compatible_neighbor_count=1,
            conflict_neighbor_count=2,
        ),
        CladeCompatibilityNodeRow(
            clade="A|C",
            tree_count=1,
            frequency=0.5,
            compatible_neighbor_count=1,
            conflict_neighbor_count=2,
        ),
        CladeCompatibilityNodeRow(
            clade="B|D",
            tree_count=1,
            frequency=0.5,
            compatible_neighbor_count=1,
            conflict_neighbor_count=2,
        ),
        CladeCompatibilityNodeRow(
            clade="C|D",
            tree_count=1,
            frequency=0.5,
            compatible_neighbor_count=1,
            conflict_neighbor_count=2,
        ),
    ]
    assert report.edges == [
        CladeCompatibilityEdgeRow(
            left_clade="A|B",
            right_clade="A|C",
            compatibility_relation="conflict",
            compatibility_reason="overlap-without-containment",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
        CladeCompatibilityEdgeRow(
            left_clade="A|B",
            right_clade="B|D",
            compatibility_relation="conflict",
            compatibility_reason="overlap-without-containment",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
        CladeCompatibilityEdgeRow(
            left_clade="A|B",
            right_clade="C|D",
            compatibility_relation="compatible",
            compatibility_reason="disjoint",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
        CladeCompatibilityEdgeRow(
            left_clade="A|C",
            right_clade="B|D",
            compatibility_relation="compatible",
            compatibility_reason="disjoint",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
        CladeCompatibilityEdgeRow(
            left_clade="A|C",
            right_clade="C|D",
            compatibility_relation="conflict",
            compatibility_reason="overlap-without-containment",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
        CladeCompatibilityEdgeRow(
            left_clade="B|D",
            right_clade="C|D",
            compatibility_relation="conflict",
            compatibility_reason="overlap-without-containment",
            left_tree_count=1,
            right_tree_count=1,
            left_frequency=0.5,
            right_frequency=0.5,
        ),
    ]


def test_compute_clade_compatibility_graph_uses_clades_not_taxa() -> None:
    report = compute_clade_compatibility_graph(
        fixture("clade_compatibility_tree_set.nwk")
    )

    assert all("|" in row.clade for row in report.nodes)
    assert {row.clade for row in report.nodes} == {"A|B", "A|C", "B|D", "C|D"}
