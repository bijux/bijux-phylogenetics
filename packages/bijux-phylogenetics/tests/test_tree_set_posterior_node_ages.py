from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    NonUltrametricTreeError,
)
from bijux_phylogenetics.trees import (
    summarize_posterior_node_ages,
    write_posterior_node_age_summary_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_summarize_posterior_node_ages_matches_clades_across_topologies() -> None:
    report = summarize_posterior_node_ages(
        fixture("posterior_node_age_summary_tree_set.nwk")
    )
    rows = {row.clade: row for row in report.rows}

    assert report.tree_count == 4
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.hpd_mass == 0.95
    assert len(report.rows) == 7

    a_b = rows["A|B"]
    assert a_b.node_kind == "internal"
    assert a_b.matched_tree_count == 4
    assert a_b.posterior_tree_count == 4
    assert a_b.clade_frequency == 1.0
    assert a_b.mean_node_age == pytest.approx(1.75)
    assert a_b.median_node_age == pytest.approx(1.75)
    assert a_b.hpd_95_lower == pytest.approx(1.0)
    assert a_b.hpd_95_upper == pytest.approx(2.5)
    assert a_b.effective_sample_size == pytest.approx(2.4)

    d_e = rows["D|E"]
    assert d_e.matched_tree_count == 3
    assert d_e.clade_frequency == pytest.approx(0.75)
    assert d_e.mean_node_age == pytest.approx(1.166666666666667)
    assert d_e.median_node_age == pytest.approx(1.0)
    assert d_e.hpd_95_lower == pytest.approx(1.0)
    assert d_e.hpd_95_upper == pytest.approx(1.5)
    assert d_e.effective_sample_size == pytest.approx(3.0)

    root = rows["A|B|C|D|E"]
    assert root.node_kind == "root"
    assert root.matched_tree_count == 4
    assert root.mean_node_age == pytest.approx(4.0)
    assert root.median_node_age == pytest.approx(4.0)
    assert root.hpd_95_lower == pytest.approx(4.0)
    assert root.hpd_95_upper == pytest.approx(4.0)
    assert root.effective_sample_size == pytest.approx(4.0)


def test_write_posterior_node_age_summary_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = summarize_posterior_node_ages(
        fixture("posterior_node_age_summary_tree_set.nwk")
    )

    output_path = tmp_path / "posterior-node-ages.tsv"
    write_posterior_node_age_summary_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines()[:3] == [
        "clade\tnode_kind\tmatched_tree_count\tposterior_tree_count\tclade_frequency\tmean_node_age\tmedian_node_age\thpd_95_lower\thpd_95_upper\teffective_sample_size",
        "A|B\tinternal\t4\t4\t1\t1.75\t1.75\t1\t2.5\t2.4",
        "A|B|C|D|E\troot\t4\t4\t1\t4\t4\t4\t4\t4",
    ]


def test_summarize_posterior_node_ages_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        summarize_posterior_node_ages(fixture("example_tree_set_mismatched.nwk"))


def test_summarize_posterior_node_ages_rejects_non_ultrametric_tree_sets(
    tmp_path: Path,
) -> None:
    tree_set_path = tmp_path / "non-ultrametric-tree-set.nwk"
    nonclock = (
        fixture("strict_clock_nonclock_tree_4_taxa.nwk")
        .read_text(encoding="utf-8")
        .strip()
    )
    tree_set_path.write_text(f"{nonclock}\n{nonclock}\n", encoding="utf-8")

    with pytest.raises(NonUltrametricTreeError) as error:
        summarize_posterior_node_ages(tree_set_path)

    assert error.value.code == "date_aware_tree_comparison_requires_ultrametric_tree"
    assert error.value.details["source_tree_index"] == 1
