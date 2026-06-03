from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_topology_rooted_nni_neighbors_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "rooted-nni-neighbors"

    exit_code = main(
        [
            "topology",
            "rooted-nni-neighbors",
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
        "algorithm": "rooted-nni-neighbor-enumeration",
        "tip_count": 4,
        "internal_node_count": 3,
        "expected_neighbor_count": 4,
        "generated_neighbor_count": 4,
        "unique_neighbor_topology_count": 4,
        "duplicate_neighbor_topology_count": 0,
        "missing_tip_taxa": 0,
        "unexpected_tip_taxa": 0,
    }
    assert (out_dir / "input_tree.nwk").is_file()
    assert (out_dir / "neighbors.tsv").is_file()
    assert (out_dir / "run.json").is_file()
