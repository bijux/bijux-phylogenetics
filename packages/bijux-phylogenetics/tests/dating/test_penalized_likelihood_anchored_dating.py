from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.dating import (
    fit_penalized_likelihood_dating,
    load_fixed_dating_calibrations,
)
from bijux_phylogenetics.phylo.dating.inputs import load_tip_dates_for_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_penalized_likelihood_dating_honors_fixed_internal_node_dates() -> None:
    tree_path = fixture(
        "trees",
        "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
    )
    metadata_path = fixture(
        "metadata",
        "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
    )
    tree = load_tree(tree_path)
    tree.rooted = True
    tip_dates, _resolved_taxon_column = load_tip_dates_for_tree(
        metadata_path,
        tree_taxa=tree.tip_names,
        taxon_column=None,
        date_column="date",
    )
    calibration_rows = load_fixed_dating_calibrations(
        tree_path,
        fixture(
            "metadata",
            "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
        ),
    )
    fixed_node_dates = {
        row.node_id: row.fixed_date
        for row in calibration_rows
        if tuple(row.descendant_taxa) in {("A", "B"), ("D", "E")}
    }

    report = fit_penalized_likelihood_dating(
        tree,
        tip_dates,
        fixed_node_dates=fixed_node_dates,
        smoothing_parameter=1.0,
        tree_path=tree_path,
        metadata_path=metadata_path,
    )

    rows_by_taxa = {tuple(row.descendant_taxa): row for row in report.node_rows}
    assert rows_by_taxa[("A", "B")].estimated_date == pytest.approx(1996.0, abs=1e-9)
    assert rows_by_taxa[("D", "E")].estimated_date == pytest.approx(1994.0, abs=1e-9)
