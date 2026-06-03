from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    generate_random_bifurcating_tree,
    validate_random_bifurcating_taxa,
    write_random_bifurcating_tree_artifacts,
)


def test_topology_gateway_exports_random_bifurcating_tree_surface() -> None:
    assert (
        topology_api.generate_random_bifurcating_tree
        is generate_random_bifurcating_tree
    )
    assert (
        topology_api.validate_random_bifurcating_taxa
        is validate_random_bifurcating_taxa
    )
    assert (
        topology_api.write_random_bifurcating_tree_artifacts
        is write_random_bifurcating_tree_artifacts
    )


def test_generate_random_bifurcating_tree_uses_requested_taxa_once_with_seeded_topology() -> (
    None
):
    left_tree, left_report = generate_random_bifurcating_tree(
        ["Gamma", "Alpha", "Delta", "Beta"],
        seed=19,
    )
    right_tree, right_report = generate_random_bifurcating_tree(
        ["Beta", "Delta", "Alpha", "Gamma"],
        seed=19,
    )

    assert left_tree.to_newick() == right_tree.to_newick()
    assert sorted(left_tree.tip_names) == ["Alpha", "Beta", "Delta", "Gamma"]
    assert left_report.algorithm == "random-bifurcating-tree-generation"
    assert left_report.seed == 19
    assert left_report.branch_length_policy == "none"
    assert left_report.requested_taxa == ["Alpha", "Beta", "Delta", "Gamma"]
    assert left_report.tip_count == 4
    assert left_report.internal_node_count == 3
    assert left_report.rooted is True
    assert left_report.strictly_bifurcating is True
    assert left_report.all_requested_taxa_present_once is True
    assert left_report.missing_requested_taxa == []
    assert left_report.duplicate_generated_taxa == []
    assert left_report.unexpected_generated_taxa == []
    assert left_report.validation_errors == []
    assert ":" not in left_report.tree_newick


def test_generate_random_bifurcating_tree_supports_uniform_branch_lengths() -> None:
    tree, report = generate_random_bifurcating_tree(
        ["Alpha", "Beta", "Gamma", "Delta"],
        seed=19,
        branch_length_policy="uniform",
    )

    assert report.branch_length_policy == "uniform"
    assert all(
        node.branch_length is None or 0.0 <= node.branch_length <= 1.0
        for node in tree.iter_nodes()
    )
    assert ":" in report.tree_newick


def test_write_random_bifurcating_tree_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    _tree, report = generate_random_bifurcating_tree(
        ["Alpha", "Beta", "Gamma", "Delta"],
        seed=19,
    )

    outputs = write_random_bifurcating_tree_artifacts(
        tmp_path / "random-bifurcating-tree-run",
        report,
    )

    assert set(outputs) == {"tree_path", "summary_path", "run_json_path"}
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "algorithm\tseed\tbranch_length_policy\trequested_taxa\ttip_order\ttip_count\tinternal_node_count\trooted\tstrictly_bifurcating\tall_requested_taxa_present_once\tmissing_requested_taxa\tduplicate_generated_taxa\tunexpected_generated_taxa\tvalidation_errors\ttree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "random-bifurcating-tree-generation"
    assert payload["seed"] == 19
    assert payload["strictly_bifurcating"] is True
    assert payload["all_requested_taxa_present_once"] is True


def test_validate_random_bifurcating_taxa_rejects_duplicates() -> None:
    try:
        validate_random_bifurcating_taxa(["Alpha", "Beta", "Alpha"])
    except ValueError as error:
        assert (
            str(error)
            == "random bifurcating tree generation requires distinct taxa; duplicates: Alpha"
        )
    else:
        raise AssertionError("duplicate taxa must be rejected")
