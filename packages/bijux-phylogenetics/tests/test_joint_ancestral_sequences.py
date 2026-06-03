from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment,
    reconstruct_nucleotide_joint_ancestral_sequences_from_alignment,
    summarize_marginal_ancestral_sites,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_joint_ancestral_sequences_differ_from_marginal_fixture() -> None:
    tree_path = fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    alignment_path = fixture(
        "alignments",
        "jc69_joint_ancestral_difference_alignment_3_taxa.fasta",
    )

    joint_report = reconstruct_nucleotide_joint_ancestral_sequences_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
    )
    marginal_report = (
        evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
            tree_path,
            alignment_path,
            model_name="jc69",
        )
    )

    assert joint_report.model_name == "JC69"
    assert joint_report.site_count == 1
    assert joint_report.internal_node_count == 2
    assert [record.sequence for record in joint_report.sequence_records] == ["A", "A"]

    marginal_most_likely_by_node_id = {
        row.node_id: row.most_likely_state
        for row in summarize_marginal_ancestral_sites(marginal_report)
        if row.site_position == 1
    }
    joint_state_by_node_id = {
        row.node_id: row.state for row in joint_report.assignment_rows
    }

    root = load_tree(tree_path).root
    child = root.children[0]
    root_id = root.node_id or ""
    child_id = child.node_id or ""
    assert marginal_most_likely_by_node_id[root_id] == "C"
    assert marginal_most_likely_by_node_id[child_id] == "A"
    assert joint_state_by_node_id[root_id] == "A"
    assert joint_state_by_node_id[child_id] == "A"


def test_joint_ancestral_sequences_expand_internal_node_site_rows() -> None:
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    report = reconstruct_nucleotide_joint_ancestral_sequences_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )

    assert report.internal_node_count == 3
    assert report.site_count == 10
    assert report.pattern_count == 6
    assert len(report.sequence_records) == 3
    assert len(report.assignment_rows) == 30
    assert {row.node_id for row in report.assignment_rows} == {
        node.node_id for node in tree.iter_internal_nodes(order="preorder")
    }
    assert {row.site_position for row in report.assignment_rows} == set(range(1, 11))
    assert {len(record.sequence) for record in report.sequence_records} == {10}


def test_joint_ancestral_sequences_support_selected_nucleotide_models() -> None:
    cases = [
        ("jc69", {}),
        ("k80", {"kappa": 2.0}),
        ("f81", {}),
        ("hky85", {"kappa": 2.5}),
        (
            "gtr",
            {
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 2.0,
                    "AT": 1.5,
                    "CG": 0.8,
                    "CT": 1.7,
                    "GT": 1.2,
                }
            },
        ),
    ]
    tree_path = fixture("trees", "jc69_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")

    for model_name, parameters in cases:
        report = reconstruct_nucleotide_joint_ancestral_sequences_from_alignment(
            tree_path,
            alignment_path,
            model_name=model_name,
            **parameters,
        )
        assert report.site_count == 4
        assert report.internal_node_count == 1
        assert len(report.sequence_records) == 1
        assert len(report.sequence_records[0].sequence) == 4
        assert len(report.assignment_rows) == 4
