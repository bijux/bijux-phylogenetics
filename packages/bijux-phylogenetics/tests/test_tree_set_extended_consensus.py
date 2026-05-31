from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    MajorityRuleExtendedAcceptedCladeRow,
    MajorityRuleExtendedConsensusReport,
    MajorityRuleExtendedRejectedCladeRow,
    compute_consensus_tree,
    compute_majority_rule_extended_consensus,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_tree_gateway_exports_extended_consensus_surface() -> None:
    assert (
        trees_api.MajorityRuleExtendedAcceptedCladeRow
        is MajorityRuleExtendedAcceptedCladeRow
    )
    assert trees_api.MajorityRuleExtendedRejectedCladeRow is (
        MajorityRuleExtendedRejectedCladeRow
    )
    assert trees_api.MajorityRuleExtendedConsensusReport is (
        MajorityRuleExtendedConsensusReport
    )
    assert (
        trees_api.compute_majority_rule_extended_consensus
        is compute_majority_rule_extended_consensus
    )


def test_compute_majority_rule_extended_consensus_reports_exact_fixture() -> None:
    tree, report = compute_majority_rule_extended_consensus(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )
    majority_tree, majority_report = compute_consensus_tree(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )

    assert dumps_newick(tree) == "(((A,B)60,C)40,(D,E)40);"
    assert dumps_newick(majority_tree) == "((A,B)60,C,D,E);"
    assert report.tree_count == 5
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.consensus_method == "majority-rule-extended"
    assert report.majority_threshold == 0.5
    assert report.included_clade_count == 3
    assert report.majority_included_clade_count == 1
    assert report.extension_included_clade_count == 2
    assert report.rejected_conflict_count == 6
    assert majority_report.included_clade_count == 1
    assert report.accepted_clades == [
        MajorityRuleExtendedAcceptedCladeRow(
            insertion_rank=1,
            clade="A|B",
            tree_count=3,
            frequency=0.6,
            inclusion_stage="majority",
        ),
        MajorityRuleExtendedAcceptedCladeRow(
            insertion_rank=2,
            clade="A|B|C",
            tree_count=2,
            frequency=0.4,
            inclusion_stage="compatible-extension",
        ),
        MajorityRuleExtendedAcceptedCladeRow(
            insertion_rank=3,
            clade="D|E",
            tree_count=2,
            frequency=0.4,
            inclusion_stage="compatible-extension",
        ),
    ]
    assert report.rejected_clades == [
        MajorityRuleExtendedRejectedCladeRow(
            clade="A|B|D",
            tree_count=2,
            frequency=0.4,
            blocking_clades=["A|B|C"],
        ),
        MajorityRuleExtendedRejectedCladeRow(
            clade="C|E",
            tree_count=2,
            frequency=0.4,
            blocking_clades=["A|B|C"],
        ),
        MajorityRuleExtendedRejectedCladeRow(
            clade="A|B|E",
            tree_count=1,
            frequency=0.2,
            blocking_clades=["A|B|C", "D|E"],
        ),
        MajorityRuleExtendedRejectedCladeRow(
            clade="A|C",
            tree_count=1,
            frequency=0.2,
            blocking_clades=["A|B"],
        ),
        MajorityRuleExtendedRejectedCladeRow(
            clade="A|D",
            tree_count=1,
            frequency=0.2,
            blocking_clades=["A|B", "A|B|C", "D|E"],
        ),
        MajorityRuleExtendedRejectedCladeRow(
            clade="C|D",
            tree_count=1,
            frequency=0.2,
            blocking_clades=["A|B|C", "D|E"],
        ),
    ]


def test_compute_majority_rule_extended_consensus_recovers_resolution_beyond_star(
    tmp_path: Path,
) -> None:
    tree_set = tmp_path / "extended-consensus-four-taxon.nwk"
    tree_set.write_text(
        "\n".join(
            [
                "((A,B),(C,D));",
                "((A,B),(C,D));",
                "((A,C),(B,D));",
                "((A,C),(B,D));",
                "((A,D),(B,C));",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    majority_tree, majority_report = compute_consensus_tree(tree_set)
    extended_tree, extended_report = compute_majority_rule_extended_consensus(tree_set)

    assert dumps_newick(majority_tree) == "(A,B,C,D);"
    assert majority_report.included_clade_count == 0
    assert dumps_newick(extended_tree) == "((A,B)40,(C,D)40);"
    assert extended_report.included_clade_count == 2
    assert [row.clade for row in extended_report.accepted_clades] == ["A|B", "C|D"]


def test_compute_majority_rule_extended_consensus_requires_exact_taxa() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="share the exact same taxon set",
    ):
        compute_majority_rule_extended_consensus(
            fixture("example_tree_set_mismatched.nwk")
        )
