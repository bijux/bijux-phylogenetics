from __future__ import annotations

import json

import pytest

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    STEPWISE_ADDITION_ROOT_BRANCH_ID,
    apply_stepwise_addition_candidate,
    build_greedy_stepwise_addition_tree,
    iter_stepwise_addition_edge_candidates,
    rooted_topology_signature_ids,
    validate_stepwise_addition_taxa,
    validate_stepwise_objective_direction,
    write_stepwise_addition_artifacts,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_topology_gateway_exports_stepwise_addition_validation_surface() -> None:
    assert (
        topology_api.validate_stepwise_addition_taxa is validate_stepwise_addition_taxa
    )
    assert (
        topology_api.validate_stepwise_objective_direction
        is validate_stepwise_objective_direction
    )
    assert (
        topology_api.STEPWISE_ADDITION_ROOT_BRANCH_ID
        == STEPWISE_ADDITION_ROOT_BRANCH_ID
    )
    assert (
        topology_api.iter_stepwise_addition_edge_candidates
        is iter_stepwise_addition_edge_candidates
    )
    assert (
        topology_api.apply_stepwise_addition_candidate
        is apply_stepwise_addition_candidate
    )
    assert (
        topology_api.build_greedy_stepwise_addition_tree
        is build_greedy_stepwise_addition_tree
    )
    assert (
        topology_api.write_stepwise_addition_artifacts
        is write_stepwise_addition_artifacts
    )


def test_validate_stepwise_addition_taxa_preserves_insertion_order() -> None:
    assert validate_stepwise_addition_taxa(["Beta", "Alpha", "Gamma"]) == [
        "Beta",
        "Alpha",
        "Gamma",
    ]


@pytest.mark.parametrize(
    ("taxa", "message"),
    [
        (["Alpha"], "stepwise addition requires at least two taxa"),
        (
            ["Alpha", "", "Gamma"],
            "stepwise addition does not allow blank taxon labels",
        ),
        (
            ["Alpha", "Beta", "Alpha"],
            "stepwise addition requires distinct taxa; duplicates: Alpha",
        ),
    ],
)
def test_validate_stepwise_addition_taxa_rejects_invalid_taxon_sets(
    taxa: list[str],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        validate_stepwise_addition_taxa(taxa)


def test_validate_stepwise_objective_direction_accepts_supported_values() -> None:
    assert validate_stepwise_objective_direction(" minimize ") == "minimize"
    assert validate_stepwise_objective_direction("MAXIMIZE") == "maximize"


def test_validate_stepwise_objective_direction_rejects_unknown_values() -> None:
    with pytest.raises(
        ValueError,
        match="objective_direction must be one of maximize, minimize",
    ):
        validate_stepwise_objective_direction("median")


def test_iter_stepwise_addition_edge_candidates_lists_root_and_all_branches() -> None:
    tree = PhyloTree(
        root=TreeNode(children=[TreeNode(name="Alpha"), TreeNode(name="Beta")]),
        rooted=True,
    )

    candidates = list(iter_stepwise_addition_edge_candidates(tree))

    assert [
        (candidate.branch_id, candidate.descendant_taxa) for candidate in candidates
    ] == [
        ("root", ("Alpha", "Beta")),
        ("Alpha", ("Alpha",)),
        ("Beta", ("Beta",)),
    ]


def test_apply_stepwise_addition_candidate_inserts_taxon_on_selected_edge() -> None:
    tree = PhyloTree(
        root=TreeNode(children=[TreeNode(name="Alpha"), TreeNode(name="Beta")]),
        rooted=True,
    )
    candidate = list(iter_stepwise_addition_edge_candidates(tree))[1]

    inserted_tree = apply_stepwise_addition_candidate(tree, candidate, "Gamma")

    assert inserted_tree.tip_names == ["Alpha", "Gamma", "Beta"]
    assert rooted_topology_signature_ids(inserted_tree) == ("Alpha|Gamma",)


def test_build_greedy_stepwise_addition_tree_records_best_scoring_edge_per_step() -> (
    None
):
    def score_tree(tree: PhyloTree) -> float:
        clade_ids = set(rooted_topology_signature_ids(tree))
        score = 0.0
        if tree.tip_count >= 3 and "Alpha|Gamma" not in clade_ids:
            score += 10.0
        if tree.tip_count >= 4 and "Beta|Delta" not in clade_ids:
            score += 5.0
        return score

    tree, report = build_greedy_stepwise_addition_tree(
        ["Alpha", "Beta", "Gamma", "Delta"],
        score_tree=score_tree,
        objective_name="clade-presence",
    )

    assert rooted_topology_signature_ids(tree) == ("Alpha|Gamma", "Beta|Delta")
    assert report.algorithm == "greedy-stepwise-addition-tree"
    assert report.objective_name == "clade-presence"
    assert report.objective_direction == "minimize"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.strictly_bifurcating is True
    assert report.all_requested_taxa_present_once is True
    assert report.final_score == 0.0
    assert len(report.trace_rows) == 2
    assert report.trace_rows[0].taxon == "Gamma"
    assert report.trace_rows[0].best_edge_id == "Alpha"
    assert len(report.trace_rows[0].tested_edge_rows) == 3
    assert report.trace_rows[1].taxon == "Delta"
    assert report.trace_rows[1].best_edge_id == "Beta"
    assert len(report.trace_rows[1].tested_edge_rows) == 5


def test_write_stepwise_addition_artifacts_materializes_governed_outputs(
    tmp_path,
) -> None:
    def score_tree(tree: PhyloTree) -> float:
        clade_ids = set(rooted_topology_signature_ids(tree))
        score = 0.0
        if tree.tip_count >= 3 and "Alpha|Gamma" not in clade_ids:
            score += 10.0
        if tree.tip_count >= 4 and "Beta|Delta" not in clade_ids:
            score += 5.0
        return score

    _tree, report = build_greedy_stepwise_addition_tree(
        ["Alpha", "Beta", "Gamma", "Delta"],
        score_tree=score_tree,
        objective_name="clade-presence",
    )

    outputs = write_stepwise_addition_artifacts(
        tmp_path / "stepwise-addition-run", report
    )

    assert set(outputs) == {"tree_path", "trace_path", "run_json_path"}
    assert (
        outputs["trace_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "step_index\ttaxon\tinserted_taxa\ttested_edge_id\ttested_edge_descendant_taxa\ttested_edge_score\tbest_edge_id\tbest_edge_descendant_taxa\tbest_score\tselected\tcandidate_tree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "greedy-stepwise-addition-tree"
    assert payload["objective_name"] == "clade-presence"
    assert payload["final_score"] == 0.0
