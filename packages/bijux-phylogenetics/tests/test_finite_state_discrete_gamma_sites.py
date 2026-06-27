from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_selected_nucleotide_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma,
    write_site_log_likelihood_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_finite_state_discrete_gamma_report_writes_site_likelihood_table(
    tmp_path: Path,
) -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="JC69",
        owner_name="finite-state shared gamma proof",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.8,
        category_count=4,
    )

    output_path = tmp_path / "jc69-gamma-site-likelihoods.tsv"
    write_site_log_likelihood_table(output_path, report)
    lines = output_path.read_text(encoding="utf-8").strip().splitlines()

    assert lines[0] == (
        "model_name\ttaxon_order\tpattern_id\tpattern_weight\tsite_position\t"
        "site_states\tlog_likelihood"
    )
    assert lines[1].startswith("JC69\tA|B\tpattern-1\t1\t1\tA|A\t")
