from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare import (
    DuplicationLossTransferAssociationRow,
    DuplicationLossTransferEventRow,
    DuplicationLossTransferReport,
    reconcile_duplication_loss_transfer,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_compare_gateway_exports_dlt_reconciliation_surface() -> None:
    assert (
        compare_api.DuplicationLossTransferAssociationRow
        is DuplicationLossTransferAssociationRow
    )
    assert (
        compare_api.DuplicationLossTransferEventRow is DuplicationLossTransferEventRow
    )
    assert compare_api.DuplicationLossTransferReport is DuplicationLossTransferReport
    assert (
        compare_api.reconcile_duplication_loss_transfer
        is reconcile_duplication_loss_transfer
    )


def test_dlt_reconciliation_reports_duplication_transfer_and_loss_events() -> None:
    report = reconcile_duplication_loss_transfer(
        fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk"),
        fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "duplication_loss_transfer_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert report.observed_species_taxa == ["A", "C", "D"]
    assert report.species_only_taxa == ["B"]
    assert report.gene_tip_count == 4
    assert report.root_mapping_branch == "C|D"
    assert report.reconciliation_score == pytest.approx(6.0, abs=1e-12)
    assert report.duplication_event_count == 1
    assert report.loss_event_count == 1
    assert report.transfer_event_count == 1
    assert report.speciation_event_count == 1
    assert [(row.gene_taxon, row.species_taxon) for row in report.mapping_rows] == [
        ("A__1", "A"),
        ("A__2", "A"),
        ("C__1", "C"),
        ("D__1", "D"),
    ]
    assert [
        (
            row.event_type,
            row.gene_node_name,
            row.mapped_species_branch,
            row.left_child_species_branch,
            row.right_child_species_branch,
            row.transferred_child_side,
            row.transfer_recipient_branch,
            row.loss_branches,
            row.event_cost,
        )
        for row in report.event_rows[:4]
    ] == [
        (
            "origin",
            None,
            "A|B|C|D",
            "C|D",
            None,
            None,
            "C|D",
            ["A|B"],
            1.0,
        ),
        (
            "transfer",
            None,
            "C|D",
            "C|D",
            "A",
            "right",
            "A",
            [],
            3.0,
        ),
        (
            "speciation",
            None,
            "C|D",
            "C",
            "D",
            None,
            None,
            [],
            0.0,
        ),
        (
            "leaf",
            "C__1",
            "C",
            None,
            None,
            None,
            None,
            [],
            0.0,
        ),
    ]
    assert report.event_rows[5].event_type == "duplication"
    assert report.event_rows[5].mapped_species_branch == "A"
    assert report.event_rows[5].left_child_species_branch == "A"
    assert report.event_rows[5].right_child_species_branch == "A"


def test_dlt_reconciliation_requires_taxon_map_for_nonmatching_gene_tips() -> None:
    with pytest.raises(
        ValueError,
        match="requires --taxon-map when gene tips do not exactly match species-tree taxa",
    ):
        reconcile_duplication_loss_transfer(
            fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk"),
            fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk"),
        )


def test_dlt_reconciliation_rejects_negative_event_costs() -> None:
    with pytest.raises(ValueError, match="duplication_cost must be nonnegative"):
        reconcile_duplication_loss_transfer(
            fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk"),
            fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk"),
            taxon_map_path=fixture(
                "metadata",
                "duplication_loss_transfer_gene_taxon_map_4_tips.tsv",
            ),
            duplication_cost=-1.0,
        )
