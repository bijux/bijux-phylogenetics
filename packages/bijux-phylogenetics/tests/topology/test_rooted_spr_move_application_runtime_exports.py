from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveApplicationReport,
    resolve_rooted_spr_move_candidate,
    summarize_rooted_spr_move_application,
    write_rooted_spr_move_artifacts,
    write_rooted_spr_move_run_json,
)


def test_public_runtime_exports_rooted_spr_move_application_surface() -> None:
    assert topology_api.RootedSprMoveApplicationReport is RootedSprMoveApplicationReport
    assert (
        topology_api.resolve_rooted_spr_move_candidate
        is resolve_rooted_spr_move_candidate
    )
    assert (
        topology_api.summarize_rooted_spr_move_application
        is summarize_rooted_spr_move_application
    )
    assert (
        topology_api.write_rooted_spr_move_artifacts is write_rooted_spr_move_artifacts
    )
    assert topology_api.write_rooted_spr_move_run_json is write_rooted_spr_move_run_json
