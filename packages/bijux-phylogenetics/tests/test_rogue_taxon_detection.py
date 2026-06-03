from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    RogueTaxonDetectionReport,
    RogueTaxonScoreRow,
    detect_rogue_taxa,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_tree_gateway_exports_rogue_taxon_surface() -> None:
    assert trees_api.RogueTaxonDetectionReport is RogueTaxonDetectionReport
    assert trees_api.RogueTaxonScoreRow is RogueTaxonScoreRow
    assert trees_api.detect_rogue_taxa is detect_rogue_taxa


def test_detect_rogue_taxa_ranks_consensus_and_rf_improvement() -> None:
    report = detect_rogue_taxa(fixture("rogue_taxon_tree_set.nwk"))

    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.consensus_threshold == 0.5
    assert (
        report.ranking_objective
        == "maximize-consensus-resolution-then-support-then-rf-stability"
    )
    assert report.baseline_consensus_newick == (
        "((A:1,B:0.1)80:0.1,(C:0.1,D:0.1)80:0.1,E:0.1);"
    )
    assert report.baseline_consensus_resolution == 0.666666666666667
    assert report.baseline_mean_support_percent == 80.0
    assert report.baseline_mean_normalized_robinson_foulds == 0.533333333333333
    assert report.baseline_rooted_topology_count == 5
    assert report.baseline_dominant_topology_frequency == 0.2

    top_row = report.rows[0]
    assert top_row.taxon == "E"
    assert top_row.rank == 1
    assert top_row.consensus_resolution_delta == 0.333333333333333
    assert top_row.mean_support_percent_delta == 20.0
    assert top_row.normalized_robinson_foulds_stability_delta == 0.533333333333333
    assert top_row.pruned_consensus_newick == (
        "((A:1.02,B:0.1)100:0.12,(C:0.1,D:0.12)100:0.12);"
    )
    assert report.rows[1].taxon == "A"
    assert report.rows[1].mean_terminal_branch_length == 1.0
    assert top_row.mean_terminal_branch_length == 0.1
    assert report.rows[1].consensus_resolution_delta < 0.0
    assert report.rows[1].normalized_robinson_foulds_stability_delta < 0.0


def test_detect_rogue_taxa_requires_at_least_four_shared_taxa(
    tmp_path: Path,
) -> None:
    tree_set = tmp_path / "small-tree-set.nwk"
    tree_set.write_text(
        "\n".join(
            [
                "((A:0.1,B:0.1):0.1,C:0.1);",
                "(A:0.1,(B:0.1,C:0.1):0.1);",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(InvalidAlignmentError, match="requires at least four taxa"):
        detect_rogue_taxa(tree_set)


def test_detect_rogue_taxa_validates_consensus_threshold() -> None:
    with pytest.raises(ValueError, match="greater than 0 and at most 1"):
        detect_rogue_taxa(
            fixture("rogue_taxon_tree_set.nwk"),
            consensus_threshold=0.0,
        )
