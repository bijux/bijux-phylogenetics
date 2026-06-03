from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.trees import (
    compute_credible_clade_set,
    write_credible_clade_set_artifacts,
    write_credible_clade_set_excluded_table,
    write_credible_clade_set_included_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compute_credible_clade_set_reports_exact_95_percent_set() -> None:
    report = compute_credible_clade_set(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )

    assert report.tree_count == 5
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.credible_threshold == 0.95
    assert report.included_clade_count == 2
    assert report.excluded_clade_count == 7
    assert report.included_cumulative_frequency == 1.0
    assert [
        (
            row.inclusion_rank,
            row.clade,
            row.tree_count,
            row.frequency,
            row.cumulative_frequency,
        )
        for row in report.included_clades
    ] == [
        (1, "A|B", 3, 0.6, 0.6),
        (2, "A|B|C", 2, 0.4, 1.0),
    ]
    assert [row.clade for row in report.excluded_clades] == [
        "A|B|D",
        "C|E",
        "D|E",
        "A|B|E",
        "A|C",
        "A|D",
        "C|D",
    ]


def test_write_credible_clade_set_tables_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = compute_credible_clade_set(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )

    included_path = tmp_path / "credible-clades.tsv"
    excluded_path = tmp_path / "excluded-clades.tsv"
    write_credible_clade_set_included_table(included_path, report)
    write_credible_clade_set_excluded_table(excluded_path, report)

    assert included_path.read_text(encoding="utf-8").splitlines() == [
        "inclusion_rank\tclade\ttree_count\tfrequency\tcumulative_frequency",
        "1\tA|B\t3\t0.6\t0.6",
        "2\tA|B|C\t2\t0.4\t1",
    ]
    assert excluded_path.read_text(encoding="utf-8").splitlines()[0] == (
        "inclusion_rank\tclade\ttree_count\tfrequency\tcumulative_frequency"
    )
    assert (
        excluded_path.read_text(encoding="utf-8")
        .splitlines()[1]
        .startswith("3\tA|B|D\t2\t0.4\t1.4")
    )


def test_write_credible_clade_set_artifacts_writes_included_and_excluded_ledgers(
    tmp_path: Path,
) -> None:
    report = compute_credible_clade_set(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )

    paths = write_credible_clade_set_artifacts(tmp_path, report)

    assert sorted(paths) == ["credible_clades_path", "excluded_clades_path"]
    assert (
        paths["credible_clades_path"]
        .read_text(encoding="utf-8")
        .startswith("inclusion_rank\tclade\ttree_count")
    )
    assert (
        paths["excluded_clades_path"]
        .read_text(encoding="utf-8")
        .startswith("inclusion_rank\tclade\ttree_count")
    )


def test_compute_credible_clade_set_rejects_invalid_threshold() -> None:
    with pytest.raises(ValueError, match="credible_threshold must be greater than 0"):
        compute_credible_clade_set(
            fixture("majority_rule_extended_consensus_tree_set.nwk"),
            credible_threshold=0.0,
        )
