from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_topology_rooted_tbr_neighbors_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "rooted-tbr-neighbors"

    exit_code = main(
        [
            "topology",
            "rooted-tbr-neighbors",
            str(fixture("spr_search_start_tree_5_taxa.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"] == {
        "algorithm": "rooted-tbr-neighbor-enumeration",
        "tip_count": 5,
        "internal_node_count": 4,
        "generated_cut_edge_count": 3,
        "generated_reconnection_count": 52,
        "identity_reconnection_count": 6,
        "generated_neighbor_count": 10,
        "unique_neighbor_topology_count": 10,
        "duplicate_reconnection_neighbor_topology_count": 10,
    }
    assert (out_dir / "input_tree.nwk").is_file()
    assert (out_dir / "neighbors.tsv").is_file()
    assert (out_dir / "run.json").is_file()
