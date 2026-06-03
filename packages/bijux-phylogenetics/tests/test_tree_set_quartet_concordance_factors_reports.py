from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    compute_gene_tree_quartet_concordance_factors,
    write_gene_tree_quartet_concordance_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_write_gene_tree_quartet_concordance_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "quartet-concordance-factors.tsv"
    report = compute_gene_tree_quartet_concordance_factors(
        fixture("quartet_concordance_species_tree_4_taxa.nwk"),
        fixture("quartet_concordance_gene_trees_4_taxa.nwk"),
    )

    write_gene_tree_quartet_concordance_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "branch_id\tleft_taxa\tright_taxa\tquartet_count_per_tree\tconcordant_quartet_count\tdiscordant_first_quartet_count\tdiscordant_second_quartet_count\tuninformative_quartet_count\tinformative_quartet_count\tconcordance_factor\tconcordant_frequency\tdiscordant_first_frequency\tdiscordant_second_frequency\tuninformative_frequency",
        "A|B::C|D\tA|B\tC|D\t1\t2\t1\t1\t1\t4\t0.5\t0.4\t0.2\t0.2\t0.2",
    ]
