from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    compute_majority_rule_extended_consensus,
    write_majority_rule_extended_consensus_artifacts,
    write_majority_rule_extended_consensus_inclusion_table,
    write_majority_rule_extended_consensus_rejected_conflict_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_majority_rule_extended_consensus_inclusion_table_writes_expected_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "inclusion-order.tsv"
    report = compute_majority_rule_extended_consensus(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )[1]

    write_majority_rule_extended_consensus_inclusion_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "insertion_rank\tclade\ttree_count\tfrequency\tinclusion_stage",
        "1\tA|B\t3\t0.6\tmajority",
        "2\tA|B|C\t2\t0.4\tcompatible-extension",
        "3\tD|E\t2\t0.4\tcompatible-extension",
    ]


def test_write_majority_rule_extended_consensus_rejected_conflict_table_writes_expected_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "rejected-conflicts.tsv"
    report = compute_majority_rule_extended_consensus(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )[1]

    write_majority_rule_extended_consensus_rejected_conflict_table(
        output_path,
        report,
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "clade\ttree_count\tfrequency\tblocking_clades",
        "A|B|D\t2\t0.4\tA|B|C",
        "C|E\t2\t0.4\tA|B|C",
        "A|B|E\t1\t0.2\tA|B|C||D|E",
        "A|C\t1\t0.2\tA|B",
        "A|D\t1\t0.2\tA|B||A|B|C||D|E",
        "C|D\t1\t0.2\tA|B|C||D|E",
    ]


def test_write_majority_rule_extended_consensus_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    _tree, report = compute_majority_rule_extended_consensus(
        fixture("majority_rule_extended_consensus_tree_set.nwk")
    )

    outputs = write_majority_rule_extended_consensus_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "consensus_tree_path",
        "inclusion_table_path",
        "rejected_conflicts_path",
    ]
    assert outputs["consensus_tree_path"].read_text(encoding="utf-8").strip() == (
        "(((A,B)60,C)40,(D,E)40);"
    )
    assert (
        outputs["inclusion_table_path"]
        .read_text(encoding="utf-8")
        .startswith("insertion_rank\tclade\ttree_count\tfrequency\tinclusion_stage\n")
    )
    assert (
        outputs["rejected_conflicts_path"]
        .read_text(encoding="utf-8")
        .startswith("clade\ttree_count\tfrequency\tblocking_clades\n")
    )
