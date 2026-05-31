from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    write_nucleotide_likelihood_starting_tree_score_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_nucleotide_likelihood_starting_tree_score_table_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_nucleotide_likelihood_starting_tree_pool_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        random_start_tree_count=2,
        random_start_tree_seed=17,
    )

    output_path = write_nucleotide_likelihood_starting_tree_score_table(
        tmp_path / "starting_trees.tsv",
        report,
    )

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert [row["tree_id"] for row in rows] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row["source_strategy"] for row in rows] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree",
        "random-tree",
    ]
    assert [row["generation_seed"] for row in rows] == ["", "", "17", "18"]
    assert all(row["model_name"] == "JC69" for row in rows)
    assert all(len(row["topology_hash"]) == 64 for row in rows)
    assert all(
        math.isfinite(float(row["starting_log_likelihood"]))
        for row in rows
    )
    assert all(
        row["substitution_parameter_policy"] == "fixed-from-model"
        for row in rows
    )
    assert all(row["substitution_parameter_values"] == "{}" for row in rows)
    assert all(row["substitution_parameter_warnings"] == "" for row in rows)
    assert rows[0]["tree_newick"].endswith(";")
