from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_nni_search_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "nni-search"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "nni-search",
            str(fixture("nni_search_start_tree_4_taxa.nwk")),
            str(fixture("nni_search_matrix.tsv")),
            "--method",
            "fitch",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"] == {
        "algorithm": "parsimony-nni-search",
        "method": "fitch",
        "taxon_count": 4,
        "character_count": 2,
        "start_score": 4.0,
        "final_score": 2.0,
        "accepted_move_count": 1,
        "evaluated_neighbor_count": 8,
        "stopping_reason": "no-improving-neighbor",
    }
    assert (out_dir / "start_tree.nwk").is_file()
    assert (out_dir / "final_tree.nwk").is_file()
    assert (out_dir / "search_trace.tsv").is_file()
    assert (out_dir / "run.json").is_file()
