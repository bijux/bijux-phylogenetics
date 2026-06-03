from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    TopologyMoveValidityReport,
    summarize_rooted_nni_move_validity,
    summarize_rooted_spr_move_validity,
    summarize_rooted_tbr_move_validity,
    write_topology_move_validity_artifacts,
    write_topology_move_validity_run_json,
)


def test_public_runtime_exports_topology_move_validity_surface() -> None:
    assert topology_api.TopologyMoveValidityReport is TopologyMoveValidityReport
    assert (
        topology_api.summarize_rooted_nni_move_validity
        is summarize_rooted_nni_move_validity
    )
    assert (
        topology_api.summarize_rooted_spr_move_validity
        is summarize_rooted_spr_move_validity
    )
    assert (
        topology_api.summarize_rooted_tbr_move_validity
        is summarize_rooted_tbr_move_validity
    )
    assert (
        topology_api.write_topology_move_validity_artifacts
        is write_topology_move_validity_artifacts
    )
    assert (
        topology_api.write_topology_move_validity_run_json
        is write_topology_move_validity_run_json
    )
