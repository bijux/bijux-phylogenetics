from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def test_phylo_parsimony_bremer_support_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    matrix_path = tmp_path / "bremer_matrix.tsv"
    tree_path = tmp_path / "reference_tree.nwk"
    out_dir = tmp_path / "bremer-support"
    matrix_path.write_text(
        (
            "taxon\tchar01_terminal_a\tchar02_clade_bd\n"
            "A\t1\t0\n"
            "B\t0\t1\n"
            "C\t0\t0\n"
            "D\t0\t1\n"
        ),
        encoding="utf-8",
    )
    tree_path.write_text("((A,(B,D)),C);\n", encoding="utf-8")

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "bremer-support",
            str(tree_path),
            str(matrix_path),
            "--method",
            "dollo",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"] == {
        "algorithm": "parsimony-bremer-support",
        "method": "dollo",
        "taxon_count": 4,
        "character_count": 2,
        "candidate_tree_count": 15,
        "reference_tree_score": 2.0,
        "optimal_score": 2.0,
        "reference_tree_is_optimal": True,
        "bremer_row_count": 2,
    }
    assert (out_dir / "reference_tree.nwk").is_file()
    assert (out_dir / "optimal_tree.nwk").is_file()
    assert (out_dir / "bremer_support.tsv").is_file()
    assert (out_dir / "run.json").is_file()
