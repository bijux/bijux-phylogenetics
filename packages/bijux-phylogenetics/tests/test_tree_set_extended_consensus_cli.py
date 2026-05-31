from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_tree_set_majority_rule_extended_consensus_writes_expected_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "extended-consensus"

    exit_code = main(
        [
            "tree-set",
            "majority-rule-extended-consensus",
            str(fixture("majority_rule_extended_consensus_tree_set.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 5
    assert payload["metrics"]["shared_taxon_count"] == 5
    assert payload["metrics"]["included_clade_count"] == 3
    assert payload["metrics"]["majority_included_clade_count"] == 1
    assert payload["metrics"]["extension_included_clade_count"] == 2
    assert payload["metrics"]["rejected_conflict_count"] == 6
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "consensus_tree.nwk",
        "inclusion_order.tsv",
        "rejected_conflicts.tsv",
    ]

    assert (out_dir / "consensus_tree.nwk").read_text(encoding="utf-8").strip() == (
        "(((A,B)60,C)40,(D,E)40);"
    )
    assert (
        (out_dir / "inclusion_order.tsv")
        .read_text(encoding="utf-8")
        .startswith("insertion_rank\tclade\ttree_count\tfrequency\tinclusion_stage\n")
    )
    assert (
        (out_dir / "rejected_conflicts.tsv")
        .read_text(encoding="utf-8")
        .startswith("clade\ttree_count\tfrequency\tblocking_clades\n")
    )
