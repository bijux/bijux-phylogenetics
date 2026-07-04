from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    reconstruct_nucleotide_marginal_ancestral_sequences_from_alignment,
    write_marginal_ancestral_sequence_fasta,
    write_marginal_ancestral_sequence_uncertainty_table,
)
from bijux_phylogenetics.phylo.topology.tree import stable_node_label

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_marginal_ancestral_fasta_marks_low_confidence_sites_by_threshold() -> None:
    report = reconstruct_nucleotide_marginal_ancestral_sequences_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        model_name="jc69",
        posterior_probability_threshold=0.8,
    )

    assert report.model_name == "JC69"
    assert report.posterior_probability_threshold == 0.8
    assert report.low_confidence_state_symbol == "N"
    assert report.internal_node_count == 1
    assert len(report.sequence_records) == 1
    assert len(report.uncertainty_rows) == 4
    assert report.sequence_records[0].sequence == "ANCN"
    assert [row.exported_state for row in report.uncertainty_rows] == [
        "A",
        "N",
        "C",
        "N",
    ]
    assert [row.low_confidence for row in report.uncertainty_rows] == [
        False,
        True,
        False,
        True,
    ]


def test_marginal_ancestral_fasta_writes_fasta_and_uncertainty_table(
    tmp_path: Path,
) -> None:
    tree_path = fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk")
    report = reconstruct_nucleotide_marginal_ancestral_sequences_from_alignment(
        tree_path,
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
        posterior_probability_threshold=0.45,
    )
    fasta_path = tmp_path / "marginal-ancestral-sequences.fasta"
    uncertainty_path = tmp_path / "marginal-ancestral-uncertainty.tsv"

    write_marginal_ancestral_sequence_fasta(fasta_path, report)
    write_marginal_ancestral_sequence_uncertainty_table(uncertainty_path, report)

    records = load_fasta_alignment(fasta_path)
    assert len(records) == 3
    assert {len(record.sequence) for record in records} == {10}
    expected_ids = {
        node.node_id
        for node in load_tree(tree_path).iter_internal_nodes(order="preorder")
    }
    assert {record.identifier for record in records} == expected_ids

    lines = uncertainty_path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "model_name",
        "node_id",
        "node_name",
        "descendant_taxa",
        "pattern_id",
        "site_position",
        "exported_state",
        "most_likely_state",
        "max_posterior_probability",
        "low_confidence",
        "posterior_probability_threshold",
        "low_confidence_state_symbol",
        "posterior_probability_a",
        "posterior_probability_c",
        "posterior_probability_g",
        "posterior_probability_t",
    ]
    assert len(lines) == 31


def test_marginal_ancestral_fasta_uses_stable_internal_node_ids() -> None:
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    report = reconstruct_nucleotide_marginal_ancestral_sequences_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )

    assert {record.node_id for record in report.sequence_records} == {
        node.node_id for node in tree.iter_internal_nodes(order="preorder")
    }
    assert stable_node_label(tree.root) in report.sequence_records[0].node_id
