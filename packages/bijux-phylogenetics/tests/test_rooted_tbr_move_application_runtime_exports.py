from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedTbrMoveApplicationReport,
    RootedTbrMoveCandidate,
    apply_rooted_tbr_move,
    resolve_rooted_tbr_move_candidate,
    summarize_rooted_tbr_move_application,
    write_rooted_tbr_move_artifacts,
    write_rooted_tbr_move_run_json,
)


def test_public_runtime_exports_rooted_tbr_move_application_surface() -> None:
    assert topology_api.RootedTbrMoveApplicationReport is RootedTbrMoveApplicationReport
    assert topology_api.RootedTbrMoveCandidate is RootedTbrMoveCandidate
    assert topology_api.apply_rooted_tbr_move is apply_rooted_tbr_move
    assert (
        topology_api.resolve_rooted_tbr_move_candidate
        is resolve_rooted_tbr_move_candidate
    )
    assert (
        topology_api.summarize_rooted_tbr_move_application
        is summarize_rooted_tbr_move_application
    )
    assert (
        topology_api.write_rooted_tbr_move_artifacts is write_rooted_tbr_move_artifacts
    )
    assert topology_api.write_rooted_tbr_move_run_json is write_rooted_tbr_move_run_json
