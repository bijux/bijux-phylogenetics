from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def test_cli_tree_set_quartet_puzzling_writes_expected_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    tree_set_path = tmp_path / "quartet-puzzling-tree-set.nwk"
    out_dir = tmp_path / "quartet-puzzling"
    tree_set_path.write_text(
        "\n".join(
            [
                "(((A,B),C),(D,E));",
                "((E,D),(C,(B,A)));",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "tree-set",
            "quartet-puzzling",
            str(tree_set_path),
            "--out-dir",
            str(out_dir),
            "--max-order-count",
            "4",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 2
    assert payload["metrics"]["shared_taxon_count"] == 5
    assert payload["metrics"]["quartet_count"] == 5
    assert payload["metrics"]["assembly_count"] == 4
    assert payload["metrics"]["unique_assembled_topology_count"] == 1
    assert payload["metrics"]["canonical_root_taxon"] == "A"
    assert (out_dir / "consensus_tree.nwk").is_file()
    assert (out_dir / "assembled_trees.nwk").is_file()
    assert (out_dir / "quartet_scores.tsv").is_file()
    assert (out_dir / "assembly_scores.tsv").is_file()
    assert (out_dir / "run.json").is_file()
