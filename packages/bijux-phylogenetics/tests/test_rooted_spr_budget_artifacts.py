from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology import (
    RootedSprEnumerationBudget,
    enumerate_rooted_spr_neighbors,
    write_rooted_spr_artifacts,
)


def test_write_rooted_spr_artifacts_records_budget_fields(tmp_path: Path) -> None:
    report = enumerate_rooted_spr_neighbors(
        loads_newick("(((A,C),B),D);"),
        budget=RootedSprEnumerationBudget(
            max_pruned_clade_count=1,
            max_regraft_target_count_per_pruned_clade=3,
        ),
    )

    outputs = write_rooted_spr_artifacts(tmp_path / "rooted-spr-neighbors", report)
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))

    assert set(outputs) == {
        "input_tree_path",
        "neighbors_path",
        "summary_path",
        "run_json_path",
    }
    assert payload["max_pruned_clade_count"] == 1
    assert payload["max_regraft_target_count_per_pruned_clade"] == 3
    assert payload["skipped_pruned_clade_count"] == 5
    assert payload["skipped_regraft_target_count"] == 3
    assert payload["skipped_budget_move_candidate_count"] == 27
    assert payload["generated_move_candidate_count"] == 3
    assert payload["identity_move_candidate_count"] == 1
    assert payload["self_regraft_candidate_count"] == 0
    assert payload["generated_neighbor_count"] == 2
    assert payload["unique_neighbor_topology_count"] == 2
    assert outputs["summary_path"].read_text(encoding="utf-8").splitlines()[1].startswith(
        "rooted-spr\trooted-spr-neighbor-enumeration\t30\t2\t0\t28\t"
    )
