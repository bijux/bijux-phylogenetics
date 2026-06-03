from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    summarize_gene_tree_conflicts,
    write_gene_tree_conflict_quartet_table,
    write_gene_tree_conflict_summary_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_write_gene_tree_conflict_tables_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = summarize_gene_tree_conflicts(
        fixture("example_tree_set_left.nwk"),
        credibility_threshold=0.3,
    )
    summary_path = tmp_path / "gene-tree-conflicts.summary.tsv"
    quartet_path = tmp_path / "gene-tree-conflicts.quartet-concordance.tsv"

    write_gene_tree_conflict_summary_table(summary_path, report)
    write_gene_tree_conflict_quartet_table(quartet_path, report)

    assert summary_path.read_text(encoding="utf-8").splitlines() == [
        "tree_count\truntime_seconds\tpeak_memory_bytes\tskipped_malformed_tree_count\tshared_taxon_count\treference_rooted_topology_id\treference_tree_frequency\tclade_count\tquartet_branch_count\ttotal_quartet_count\tconflict_count\tconflicting_clade_count\trogue_taxon_count\ttop_ranked_rogue_taxon\tcredibility_threshold\trogue_consensus_threshold",
        "\t".join(
            [
                "3",
                format(report.processing.runtime_seconds, ".15g"),
                str(report.processing.peak_memory_bytes),
                "0",
                "4",
                report.reference_tree.rooted_topology_id,
                "0.666666666666667",
                "4",
                "1",
                "3",
                "4",
                "4",
                "4",
                "A",
                "0.3",
                "0.5",
            ]
        ),
    ]
    assert quartet_path.read_text(encoding="utf-8").splitlines() == [
        "branch_id\tleft_taxa\tright_taxa\tquartet_count_per_tree\tconcordant_quartet_count\tdiscordant_first_quartet_count\tdiscordant_second_quartet_count\tuninformative_quartet_count\tinformative_quartet_count\tconcordance_factor\tconcordant_frequency\tdiscordant_first_frequency\tdiscordant_second_frequency\tuninformative_frequency",
        "A|B::C|D\tA|B\tC|D\t1\t2\t1\t0\t0\t3\t0.666666666666667\t0.666666666666667\t0.333333333333333\t0\t0",
    ]
