from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprEnumerationBudget,
    summarize_rooted_spr_move_application,
    write_rooted_spr_move_artifacts,
    write_rooted_spr_move_run_json,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_rooted_spr_move_artifact_surface() -> None:
    assert (
        topology_api.write_rooted_spr_move_artifacts is write_rooted_spr_move_artifacts
    )
    assert topology_api.write_rooted_spr_move_run_json is write_rooted_spr_move_run_json


def test_write_rooted_spr_move_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = summarize_rooted_spr_move_application(
        fixture("example_tree.nwk"),
        1,
        budget=RootedSprEnumerationBudget(
            max_pruned_clade_count=1,
            max_regraft_target_count_per_pruned_clade=3,
        ),
    )

    outputs = write_rooted_spr_move_artifacts(tmp_path / "rooted-spr-move", report)

    assert set(outputs) == {
        "input_tree_path",
        "moved_tree_path",
        "run_json_path",
    }
    assert outputs["input_tree_path"].read_text(encoding="utf-8").strip().endswith(";")
    assert outputs["moved_tree_path"].read_text(encoding="utf-8").strip().endswith(";")

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "rooted-spr-move-application"
    assert payload["input_tree_path"].endswith("example_tree.nwk")
    assert payload["selected_move_index"] == 1
    assert payload["available_move_count"] == 3
    assert payload["max_pruned_clade_count"] == 1
    assert payload["max_regraft_target_count_per_pruned_clade"] == 3
    assert payload["moved_topology_changed"] is True
    assert payload["affected_subtrees"]["affected_branch_clade_ids"] == (
        report.affected_subtree_report.affected_branch_clade_ids
    )
    assert payload["affected_subtrees"]["unaffected_branch_clade_ids"] == (
        report.affected_subtree_report.unaffected_branch_clade_ids
    )
    assert payload["selected_pruned_clade_id"] == report.selected_pruned_clade_id
    assert (
        payload["selected_regraft_target_branch_id"]
        == report.selected_regraft_target_branch_id
    )
    assert payload["affected_clade_ids"] == report.affected_clade_ids
