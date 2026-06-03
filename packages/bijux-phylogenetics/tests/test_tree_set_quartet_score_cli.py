from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_tree_set_quartet_score_writes_expected_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    output_path = tmp_path / "quartet-score.tsv"

    exit_code = main(
        [
            "tree-set",
            "quartet-score",
            str(fixture("quartet_score_candidate_high_4_taxa.nwk")),
            str(fixture("quartet_concordance_gene_trees_4_taxa.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 5
    assert payload["metrics"]["shared_taxon_count"] == 4
    assert payload["metrics"]["branch_count"] == 1
    assert payload["metrics"]["quartet_score"] == 2
    assert payload["metrics"]["normalized_quartet_score"] == 0.5
    assert payload["metrics"]["concordant_quartet_count"] == 2
    assert output_path.is_file()
