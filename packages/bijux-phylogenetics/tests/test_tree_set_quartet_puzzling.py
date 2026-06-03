from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology.clades import informative_unrooted_splits
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    QuartetPuzzlingAssemblyRow,
    QuartetPuzzlingReport,
    QuartetTopologyScoreRow,
    build_quartet_puzzling_consensus,
    write_quartet_puzzling_artifacts,
)


def test_package_tree_gateway_exports_quartet_puzzling_surface() -> None:
    assert trees_api.QuartetTopologyScoreRow is QuartetTopologyScoreRow
    assert trees_api.QuartetPuzzlingAssemblyRow is QuartetPuzzlingAssemblyRow
    assert trees_api.QuartetPuzzlingReport is QuartetPuzzlingReport
    assert (
        trees_api.build_quartet_puzzling_consensus is build_quartet_puzzling_consensus
    )
    assert (
        trees_api.write_quartet_puzzling_artifacts is write_quartet_puzzling_artifacts
    )


def test_quartet_puzzling_recovers_expected_tree_from_compatible_quartets(
    tmp_path: Path,
) -> None:
    tree_set_path = tmp_path / "compatible-quartet-tree-set.nwk"
    tree_set_path.write_text(
        "\n".join(
            [
                "(((A,B),C),(D,E));",
                "((E,D),(C,(B,A)));",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    consensus_tree, report = build_quartet_puzzling_consensus(
        tree_set_path,
        max_order_count=4,
    )

    assert report.tree_count == 2
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.quartet_count == 5
    assert report.assembly_count == 4
    assert report.unique_assembled_topology_count == 1
    assert report.canonical_root_taxon == "A"
    assert report.canonical_rooting_strategy == "lexicographic-tip-outgroup"
    assert len(report.quartet_rows) == 5
    assert all(row.best_split_support_frequency == 1.0 for row in report.quartet_rows)
    assert all(
        row.tied_best_split_taxa == [row.best_split_taxa] for row in report.quartet_rows
    )
    assert len({row.assembled_tree_newick for row in report.assembly_rows}) == 1
    assert {
        tuple(sorted(split))
        for split in informative_unrooted_splits(
            consensus_tree, set(report.shared_taxa)
        )
    } == {
        ("A", "B"),
        ("D", "E"),
    }
    assert {
        tuple(sorted(split))
        for split in informative_unrooted_splits(
            loads_newick(report.consensus_newick),
            set(report.shared_taxa),
        )
    } == {
        ("A", "B"),
        ("D", "E"),
    }


def test_write_quartet_puzzling_artifacts_materializes_bundle(tmp_path: Path) -> None:
    tree_set_path = tmp_path / "compatible-quartet-tree-set.nwk"
    tree_set_path.write_text(
        "\n".join(
            [
                "(((A,B),C),(D,E));",
                "((E,D),(C,(B,A)));",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _consensus_tree, report = build_quartet_puzzling_consensus(
        tree_set_path,
        max_order_count=4,
    )

    outputs = write_quartet_puzzling_artifacts(tmp_path / "quartet-puzzling", report)

    assert set(outputs) == {
        "consensus_tree_path",
        "assembled_trees_path",
        "quartet_scores_path",
        "assembly_scores_path",
        "run_json_path",
    }
    assert (
        outputs["quartet_scores_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "quartet_taxa\tfirst_split_taxa\tfirst_split_tree_count\tsecond_split_taxa\tsecond_split_tree_count\tthird_split_taxa\tthird_split_tree_count\tuninformative_tree_count\tbest_split_taxa\tbest_split_support_frequency\ttied_best_split_taxa\n"
        )
    )
    assert (
        outputs["assembly_scores_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "order_index\ttaxon_order\tquartet_score\tassembled_topology_id\tassembled_tree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["quartet_count"] == 5
    assert payload["assembly_count"] == 4
    assert payload["unique_assembled_topology_count"] == 1
    assert payload["canonical_root_taxon"] == "A"
