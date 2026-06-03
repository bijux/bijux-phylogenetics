from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    resolve_rooted_spr_move_candidate,
    summarize_rooted_spr_move_validity,
    write_topology_move_validity_artifacts,
    write_topology_move_validity_run_json,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_topology_gateway_exports_move_validity_artifact_surface() -> None:
    assert (
        topology_api.write_topology_move_validity_artifacts
        is write_topology_move_validity_artifacts
    )
    assert (
        topology_api.write_topology_move_validity_run_json
        is write_topology_move_validity_run_json
    )


def test_write_topology_move_validity_artifacts_materializes_rejected_run_json(
    tmp_path: Path,
) -> None:
    tree = loads_newick("(((A,C),B),D);")
    resolved_candidate, _available_move_count = resolve_rooted_spr_move_candidate(
        tree, 1
    )
    candidate = RootedSprMoveCandidate(
        pruned_node_id=resolved_candidate.pruned_node_id,
        pruned_clade_id=resolved_candidate.pruned_clade_id,
        pruned_descendant_taxa=resolved_candidate.pruned_descendant_taxa,
        regraft_target_branch_id=resolved_candidate.pruned_clade_id,
        regraft_target_descendant_taxa=resolved_candidate.pruned_descendant_taxa,
    )
    report = summarize_rooted_spr_move_validity(fixture("example_tree.nwk"), candidate)

    outputs = write_topology_move_validity_artifacts(
        tmp_path / "topology-move-validity",
        report,
    )

    assert set(outputs) == {"input_tree_path", "run_json_path"}
    assert outputs["input_tree_path"].read_text(encoding="utf-8").strip().endswith(";")
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "topology-move-validity"
    assert payload["move_family"] == "rooted-spr"
    assert payload["input_tree_path"].endswith("example_tree.nwk")
    assert payload["validity_decision"] == "rejected"
    assert payload["rejection_code"] == "topology_move_self_regraft"
    assert payload["evidence"] == report.evidence
    assert (
        payload["candidate_payload"]["regraft_target_branch_id"]
        == (report.candidate_payload["regraft_target_branch_id"])
    )
