from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_compare_agreement_subtree_writes_review_bundle(
    tmp_path: Path,
    capsys,
) -> None:
    output_dir = tmp_path / "agreement"

    exit_code = main(
        [
            "compare",
            "agreement-subtree",
            str(fixture("agreement_subtree_left.nwk")),
            str(fixture("agreement_subtree_right.nwk")),
            "--out",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["shared_taxa"] == 5
    assert payload["metrics"]["retained_taxa"] == 4
    assert payload["metrics"]["agreement_removed_taxa"] == 1
    assert payload["metrics"]["evaluated_candidate_count"] == 4
    assert payload["metrics"]["topology_equal_after_pruning"] is True
    assert payload["metrics"]["post_pruning_robinson_foulds_distance"] == 0
    assert payload["data"]["agreement_removed_taxa"] == ["C"]
    assert payload["data"]["retained_taxa"] == ["A", "B", "D", "E"]
    assert payload["outputs"] == [
        str(output_dir / "left-agreement.nwk"),
        str(output_dir / "right-agreement.nwk"),
        str(output_dir / "agreement-subtree-pruning.tsv"),
        str(output_dir / "agreement-subtree-removed.tsv"),
        str(output_dir / "agreement-subtree-search.tsv"),
        str(output_dir / "agreement-subtree-comparison.tsv"),
    ]
    assert (output_dir / "left-agreement.nwk").read_text(encoding="utf-8") == (
        "((A:0.1,B:0.1):0.2,(D:0.1,E:0.1):0.2);\n"
    )
    assert (output_dir / "right-agreement.nwk").read_text(encoding="utf-8") == (
        "((A:0.1,B:0.1):0.2,(D:0.1,E:0.1):0.2);\n"
    )
    assert (
        (output_dir / "agreement-subtree-pruning.tsv")
        .read_text(encoding="utf-8")
        .startswith("tree_side\ttree_path\trf_mode\tsearch_strategy\t")
    )
    assert (output_dir / "agreement-subtree-removed.tsv").read_text(
        encoding="utf-8"
    ) == (
        "tree_side\ttree_path\ttaxon\treason\tshared_taxon\tremoved_for_agreement_subtree\n"
        f"left\t{fixture('agreement_subtree_left.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
        f"right\t{fixture('agreement_subtree_right.nwk')}\tC\tnot_requested\ttrue\ttrue\n"
    )
    assert (
        (output_dir / "agreement-subtree-search.tsv")
        .read_text(encoding="utf-8")
        .startswith(
            "candidate_index\tretained_taxon_count\tretained_taxa\tremoved_taxa\t"
        )
    )
    assert (
        (output_dir / "agreement-subtree-comparison.tsv")
        .read_text(encoding="utf-8")
        .startswith(
            "split_id\tcomparison_status\tshared_clade\tleft_support\tright_support\t"
        )
    )
