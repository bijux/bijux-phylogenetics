from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_topology_rooted_spr_neighbors_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "rooted-spr-neighbors"

    exit_code = main(
        [
            "topology",
            "rooted-spr-neighbors",
            str(fixture("example_tree.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"] == {
        "algorithm": "rooted-spr-neighbor-enumeration",
        "tip_count": 4,
        "internal_node_count": 3,
        "max_pruned_clade_count": None,
        "max_regraft_target_count_per_pruned_clade": None,
        "skipped_pruned_clade_count": 0,
        "skipped_regraft_target_count": 0,
        "generated_move_candidate_count": 32,
        "identity_move_candidate_count": 8,
        "self_regraft_candidate_count": 0,
        "generated_neighbor_count": 12,
        "unique_neighbor_topology_count": 12,
        "duplicate_move_neighbor_topology_count": 4,
    }
    assert (out_dir / "input_tree.nwk").is_file()
    assert (out_dir / "neighbors.tsv").is_file()
    assert (out_dir / "run.json").is_file()
