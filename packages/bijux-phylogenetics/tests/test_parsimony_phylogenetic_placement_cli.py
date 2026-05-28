from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.io.newick import load_newick_tree_set

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_placement_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "parsimony-placement"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "placement",
            str(fixture("placement_reference_tree_4_taxa.nwk")),
            str(fixture("placement_reference_matrix.tsv")),
            str(fixture("placement_query_matrix.tsv")),
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
        "algorithm": "parsimony-placement",
        "method": "unordered-fitch",
        "reference_taxon_count": 4,
        "character_count": 2,
        "edge_count": 6,
        "query_count": 2,
        "reference_total_steps": 2,
        "placement_count": 12,
        "equally_best_placement_count": 5,
    }
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "alternative_placements.tsv").is_file()
    assert (out_dir / "equally_best_placements.nwk").is_file()
    assert (out_dir / "run.json").is_file()
    assert len(load_newick_tree_set(out_dir / "equally_best_placements.nwk")) == 5
