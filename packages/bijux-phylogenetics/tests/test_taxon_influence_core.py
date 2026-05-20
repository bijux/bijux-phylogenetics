from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.influence import analyze_taxon_influence


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_analyze_taxon_influence_ranks_taxa_by_leave_one_out_change() -> None:
    report = analyze_taxon_influence(
        fixture("example_tree_taxon_influence_left.nwk"),
        fixture("example_tree_taxon_influence_right.nwk"),
    )

    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert [row.influence_rank for row in report.rows] == [1, 2, 3, 4, 5]
    assert report.rows[0].taxon == "C"
    assert report.rows[0].retained_taxa == ["A", "B", "D", "E"]
    assert report.rows[0].baseline_topology_equal is False
    assert report.rows[0].leave_one_out_topology_equal is True
    assert report.rows[0].rooted_robinson_foulds_delta == -2
    assert report.rows[0].conflicting_clade_delta == -2
    assert report.rows[0].high_support_conflict_delta == -1
    assert report.rows[0].topology_changed is True
    assert report.rows[0].support_changed is True
    assert report.rows[0].influence_score > report.rows[1].influence_score


def test_analyze_taxon_influence_preserves_baseline_support_surface() -> None:
    report = analyze_taxon_influence(
        fixture("example_tree_taxon_influence_left.nwk"),
        fixture("example_tree_taxon_influence_right.nwk"),
    )

    assert report.baseline_topology.rooted_robinson_foulds_distance == 2
    assert report.baseline_topology.unrooted_robinson_foulds_distance == 0
    assert len(report.baseline_support.shared_clades) == 2
    assert len(report.baseline_support.conflicting_clades) == 2
