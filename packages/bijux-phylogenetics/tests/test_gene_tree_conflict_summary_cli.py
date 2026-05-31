from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_tree_set_gene_tree_conflicts_writes_expected_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    output_dir = tmp_path / "gene-tree-conflicts"

    exit_code = main(
        [
            "tree-set",
            "gene-tree-conflicts",
            str(fixture("example_tree_set_left.nwk")),
            "--out-dir",
            str(output_dir),
            "--prefix",
            "gene-tree-review",
            "--credibility-threshold",
            "0.3",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["shared_taxon_count"] == 4
    assert payload["metrics"]["reference_tree_frequency"] == 0.666666666666667
    assert payload["metrics"]["clade_count"] == 4
    assert payload["metrics"]["quartet_branch_count"] == 1
    assert payload["metrics"]["conflict_count"] == 4
    assert payload["metrics"]["rogue_taxon_count"] == 4
    assert payload["metrics"]["top_ranked_rogue_taxon"] == "A"
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "gene-tree-review.clade-conflicts.tsv",
        "gene-tree-review.clade-frequencies.tsv",
        "gene-tree-review.quartet-concordance.tsv",
        "gene-tree-review.reference-tree.nwk",
        "gene-tree-review.rogue-taxa.tsv",
        "gene-tree-review.summary.tsv",
    ]
    assert (
        (output_dir / "gene-tree-review.clade-conflicts.tsv")
        .read_text(encoding="utf-8")
        .splitlines()[1]
        .startswith("A|B\t0.666666666666667")
    )
