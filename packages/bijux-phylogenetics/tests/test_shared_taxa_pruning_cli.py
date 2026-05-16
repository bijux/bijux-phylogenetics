from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_compare_prune_writes_review_bundle(tmp_path: Path, capsys) -> None:
    output_dir = tmp_path / "shared"

    exit_code = main(
        [
            "compare",
            "prune",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_overlap.nwk")),
            "--out",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["shared_taxa"] == 3
    assert payload["metrics"]["left_removed_taxa"] == 1
    assert payload["metrics"]["right_removed_taxa"] == 1
    assert payload["metrics"]["topology_equal_after_pruning"] is True
    assert payload["metrics"]["post_pruning_robinson_foulds_distance"] == 0
    assert payload["data"]["post_pruning_comparison"]["left_only_taxa"] == []
    assert payload["data"]["left_pruning"]["removed_taxa"] == ["D"]
    assert payload["data"]["right_pruning"]["removed_taxa"] == ["E"]
    assert payload["outputs"] == [
        str(output_dir / "left-shared.nwk"),
        str(output_dir / "right-shared.nwk"),
        str(output_dir / "shared-taxa-pruning.tsv"),
        str(output_dir / "shared-taxa-removed.tsv"),
        str(output_dir / "shared-taxa-comparison.tsv"),
    ]
    assert (output_dir / "left-shared.nwk").read_text(
        encoding="utf-8"
    ) == "((A:0.1,B:0.1):0.2,C:0.3);\n"
    assert (output_dir / "right-shared.nwk").read_text(
        encoding="utf-8"
    ) == "((A:0.1,B:0.1):0.2,C:0.3);\n"
    assert (
        (output_dir / "shared-taxa-pruning.tsv")
        .read_text(encoding="utf-8")
        .startswith("tree_side\ttree_path\toriginal_tip_count\tretained_tip_count\t")
    )
    assert (output_dir / "shared-taxa-removed.tsv").read_text(encoding="utf-8") == (
        "tree_side\ttree_path\ttaxon\treason\n"
        f"left\t{fixture('example_tree.nwk')}\tD\tnot_requested\n"
        f"right\t{fixture('example_tree_overlap.nwk')}\tE\tnot_requested\n"
    )
    assert (
        (output_dir / "shared-taxa-comparison.tsv")
        .read_text(encoding="utf-8")
        .startswith("split_id\tcomparison_status\tshared_clade\tleft_support\t")
    )
