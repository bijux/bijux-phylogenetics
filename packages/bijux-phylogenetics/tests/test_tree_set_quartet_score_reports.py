from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    compute_candidate_tree_quartet_score,
    write_candidate_tree_quartet_score_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_write_candidate_tree_quartet_score_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "quartet-score.tsv"
    report = compute_candidate_tree_quartet_score(
        fixture("quartet_score_candidate_high_4_taxa.nwk"),
        fixture("quartet_concordance_gene_trees_4_taxa.nwk"),
    )

    write_candidate_tree_quartet_score_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "branch_id\tleft_taxa\tright_taxa\tquartet_count_per_tree\tconcordant_quartet_count\tdiscordant_first_quartet_count\tdiscordant_second_quartet_count\tuninformative_quartet_count\tinformative_quartet_count\tconcordance_factor\tquartet_score\tnormalized_quartet_score",
        "A|B::C|D\tA|B\tC|D\t1\t2\t1\t1\t1\t4\t0.5\t2\t0.5",
    ]
