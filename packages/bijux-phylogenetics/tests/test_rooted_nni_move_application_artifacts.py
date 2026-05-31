from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    summarize_rooted_nni_move_application,
    write_rooted_nni_move_artifacts,
    write_rooted_nni_move_run_json,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_rooted_nni_move_artifact_surface() -> None:
    assert topology_api.write_rooted_nni_move_artifacts is write_rooted_nni_move_artifacts
    assert topology_api.write_rooted_nni_move_run_json is write_rooted_nni_move_run_json


def test_write_rooted_nni_move_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = summarize_rooted_nni_move_application(fixture("example_tree.nwk"), 1)

    outputs = write_rooted_nni_move_artifacts(tmp_path / "rooted-nni-move", report)

    assert set(outputs) == {
        "input_tree_path",
        "moved_tree_path",
        "reversed_tree_path",
        "run_json_path",
    }
    assert outputs["input_tree_path"].read_text(encoding="utf-8").strip().endswith(";")
    assert outputs["moved_tree_path"].read_text(encoding="utf-8").strip().endswith(";")
    assert outputs["reversed_tree_path"].read_text(encoding="utf-8").strip().endswith(";")

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "rooted-nni-move-application"
    assert payload["input_tree_path"].endswith("example_tree.nwk")
    assert payload["selected_move_index"] == 1
    assert payload["available_move_count"] == 4
    assert payload["moved_topology_changed"] is True
    assert payload["reverse_restores_original_topology"] is True
    assert payload["node_metadata_preserved"] is True
    assert payload["edge_metadata_preserved"] is True
    assert payload["branch_lengths_preserved"] is True

