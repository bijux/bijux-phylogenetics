from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.trees import load_tree
import bijux_phylogenetics.phylo.dating as dating_api
from bijux_phylogenetics.phylo.dating import fit_least_squares_dating_from_metadata
from bijux_phylogenetics.runtime.errors import (
    PhylogeneticsError,
    UnrootedTreeError,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_package_dating_gateway_exports_least_squares_surface() -> None:
    assert (
        dating_api.fit_least_squares_dating_from_metadata
        is fit_least_squares_dating_from_metadata
    )


def test_least_squares_dating_recovers_simulated_node_ages_within_tolerance() -> None:
    report = fit_least_squares_dating_from_metadata(
        fixture("least_squares_dating_substitution_tree_4_taxa.nwk"),
        fixture("least_squares_dating_tip_dates_4_taxa.tsv"),
    )

    assert report.tree_newick == "(((A:0.25,B:0.25):0.5,C:1):1,D:2.25);"
    assert report.dated_tree_newick != report.tree_newick
    assert report.taxa == ["A", "B", "C", "D"]
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.branch_count == 6
    assert report.parameter_count == 4
    assert report.tree_path == str(
        fixture("least_squares_dating_substitution_tree_4_taxa.nwk")
    )
    assert report.metadata_path == str(
        fixture("least_squares_dating_tip_dates_4_taxa.tsv")
    )
    assert report.taxon_column == "taxon"
    assert report.date_column == "date"
    assert report.minimum_tip_date == 2007.0
    assert report.maximum_tip_date == 2009.0
    assert report.root_date == pytest.approx(2000.0, abs=1e-6)
    assert report.estimated_clock_rate == pytest.approx(0.25, abs=1e-9)
    assert report.residual_sum_squares == pytest.approx(0.0, abs=1e-12)
    assert report.condition_number > 1.0
    assert report.exact_fit is True
    assert report.optimizer_name == "closed-form-linear-least-squares"
    assert report.converged is True

    node_rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.node_rows
    }
    expected_node_dates = {
        ("A", "B", "C", "D"): 2000.0,
        ("A", "B", "C"): 2004.0,
        ("A", "B"): 2006.0,
        ("A",): 2007.0,
        ("B",): 2007.0,
        ("C",): 2008.0,
        ("D",): 2009.0,
    }
    assert set(node_rows_by_descendant_taxa) == set(expected_node_dates)
    for descendant_taxa, expected_date in expected_node_dates.items():
        row = node_rows_by_descendant_taxa[descendant_taxa]
        assert row.estimated_date == pytest.approx(expected_date, abs=1e-6)

    branch_rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.branch_rows
    }
    expected_time_tree = load_tree(fixture("least_squares_dating_time_tree_4_taxa.nwk"))
    expected_branch_lengths = {
        tuple(child.descendant_taxa): float(child.branch_length or 0.0)
        for _parent, child in expected_time_tree.iter_edges()
    }
    expected_observed_branch_lengths = {
        ("A", "B", "C"): 1.0,
        ("A", "B"): 0.5,
        ("A",): 0.25,
        ("B",): 0.25,
        ("C",): 1.0,
        ("D",): 2.25,
    }
    assert set(branch_rows_by_descendant_taxa) == set(expected_branch_lengths)
    for descendant_taxa, expected_duration in expected_branch_lengths.items():
        row = branch_rows_by_descendant_taxa[descendant_taxa]
        assert row.fitted_time_duration == pytest.approx(expected_duration, abs=1e-6)
        assert row.fitted_branch_length == pytest.approx(
            expected_observed_branch_lengths[descendant_taxa],
            abs=1e-8,
        )
        assert row.residual == pytest.approx(0.0, abs=1e-8)


def test_least_squares_dating_rejects_unrooted_tree(tmp_path: Path) -> None:
    metadata_path = tmp_path / "tip-dates.tsv"
    metadata_path.write_text(
        "taxon\tdate\nA\t0\nB\t0\nC\t1\nD\t2\n",
        encoding="utf-8",
    )

    with pytest.raises(UnrootedTreeError, match="tree is not rooted"):
        fit_least_squares_dating_from_metadata(
            fixture("example_tree_unrooted.nwk"),
            metadata_path,
        )


def test_least_squares_dating_requires_variation_in_tip_dates(tmp_path: Path) -> None:
    metadata_path = tmp_path / "constant-dates.tsv"
    metadata_path.write_text(
        "taxon\tdate\nA\t4\nB\t4\nC\t4\nD\t4\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PhylogeneticsError,
        match="least-squares dating requires variation in tip dates",
    ):
        fit_least_squares_dating_from_metadata(
            fixture("least_squares_dating_substitution_tree_4_taxa.nwk"),
            metadata_path,
        )
