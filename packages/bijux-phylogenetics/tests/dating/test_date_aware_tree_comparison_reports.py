from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import write_date_aware_tree_comparison_table


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_date_aware_tree_comparison_table_writes_matched_clade_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "clade-ages.tsv"

    write_date_aware_tree_comparison_table(
        output_path,
        fixture("strict_clock_time_tree_4_taxa.nwk"),
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk"),
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "clade_id\tnode_kind\ttaxon_count\tdescendant_taxa\tleft_age\tright_age\tage_difference\tabsolute_age_difference\tage_rmse\tunstable_age\tcomparison_scope\ttaxon_overlap_policy\ttopology_equal\trobinson_foulds_distance",
        "A|B\tinternal\t2\tA|B\t1.0\t2.0\t1.0\t1.0\t3.16227766016838\tfalse\tfull-taxa\tprune-to-shared\ttrue\t0",
        "A|B|C\tinternal\t3\tA|B|C\t2.0\t4.0\t2.0\t2.0\t3.16227766016838\tfalse\tfull-taxa\tprune-to-shared\ttrue\t0",
        "A|B|C|D\troot\t4\tA|B|C|D\t3.0\t8.0\t5.0\t5.0\t3.16227766016838\ttrue\tfull-taxa\tprune-to-shared\ttrue\t0",
    ]
