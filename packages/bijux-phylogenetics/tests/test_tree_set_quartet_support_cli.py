from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def test_cli_tree_set_quartet_support_writes_expected_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    reference_tree = tmp_path / "quartet-reference-tree.nwk"
    comparison_tree_set = tmp_path / "quartet-support-tree-set.nwk"
    output_path = tmp_path / "quartet-support.tsv"
    reference_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    comparison_tree_set.write_text(
        "\n".join(
            [
                "((A,B),(C,D));",
                "((A,C),(B,D));",
                "((A,D),(B,C));",
                "(A,B,C,D);",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "tree-set",
            "quartet-support",
            str(reference_tree),
            str(comparison_tree_set),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 4
    assert payload["metrics"]["shared_taxon_count"] == 4
    assert payload["metrics"]["branch_count"] == 1
    assert payload["metrics"]["total_quartet_count"] == 4
    assert payload["metrics"]["concordant_quartet_count"] == 1
    assert payload["metrics"]["discordant_quartet_count"] == 2
    assert payload["metrics"]["uninformative_quartet_count"] == 1
    assert output_path.is_file()
