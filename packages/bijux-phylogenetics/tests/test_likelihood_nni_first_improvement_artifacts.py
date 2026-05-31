from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
    write_nucleotide_likelihood_nni_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_nucleotide_likelihood_nni_artifacts_records_first_improvement_policy(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        improvement_policy="first-improvement",
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_nni_artifacts(
        tmp_path / "likelihood-nni-first-improvement-run",
        report,
    )

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    selected_rows = [
        row for row in payload["candidate_rows"] if row["selected_best_move"]
    ]

    assert payload["improvement_policy"] == "first-improvement"
    assert len(payload["candidate_rows"]) == 10
    assert len(selected_rows) == 2
    assert [(row["iteration"], row["candidate_order"]) for row in selected_rows] == [
        (1, 4),
        (2, 2),
    ]
    assert [
        row["candidate_order"]
        for row in payload["candidate_rows"]
        if row["iteration"] == 2
    ] == [
        1,
        2,
    ]
    assert (
        outputs["candidate_table_path"]
        .read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("iteration\tcandidate_order\tpivot_branch_id")
    )
