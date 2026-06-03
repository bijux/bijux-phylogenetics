from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.phylo.topology.rooted_tbr import (
    summarize_rooted_tbr_move_application,
    write_rooted_tbr_move_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_rooted_tbr_move_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = summarize_rooted_tbr_move_application(
        fixture("parsimony", "spr_search_start_tree_5_taxa.nwk"),
        39,
    )

    outputs = write_rooted_tbr_move_artifacts(tmp_path / "rooted-tbr-move", report)

    assert set(outputs) == {
        "input_tree_path",
        "moved_tree_path",
        "run_json_path",
    }
    assert outputs["input_tree_path"].read_text(encoding="utf-8").strip().endswith(";")
    assert outputs["moved_tree_path"].read_text(encoding="utf-8").strip().endswith(";")

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "rooted-tbr-move-application"
    assert payload["input_tree_path"].endswith("spr_search_start_tree_5_taxa.nwk")
    assert payload["selected_move_index"] == 39
    assert payload["available_move_count"] == 46
    assert payload["selected_cut_edge_id"] == "A|B|C|D"
    assert payload["left_component_tip_count"] == 4
    assert payload["right_component_tip_count"] == 1
    assert payload["selected_left_attachment_branch_id"] == "A"
    assert payload["selected_right_attachment_branch_id"] == "interface"
    assert payload["moved_topology_changed"] is True
    assert payload["reverse_move_available"] is True
    assert payload["reverse_available_move_count"] == 2
    assert payload["affected_subtrees"]["affected_branch_clade_ids"] == (
        report.affected_subtree_report.affected_branch_clade_ids
    )
