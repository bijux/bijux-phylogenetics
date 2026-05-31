from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_tree_set_rogue_taxa_writes_review_bundle(
    tmp_path: Path,
    capsys,
) -> None:
    output_dir = tmp_path / "rogue-taxa"

    exit_code = main(
        [
            "tree-set",
            "rogue-taxa",
            str(fixture("rogue_taxon_tree_set.nwk")),
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["candidate_taxon_count"] == 5
    assert payload["metrics"]["top_ranked_taxon"] == "E"
    assert payload["metrics"]["consensus_threshold"] == 0.5
    assert payload["metrics"]["baseline_consensus_resolution"] == 0.666666666666667
    assert payload["metrics"]["top_ranked_consensus_resolution"] == 1.0
    assert payload["metrics"]["baseline_mean_support_percent"] == 80.0
    assert payload["metrics"]["top_ranked_mean_support_percent"] == 100.0
    assert payload["metrics"]["baseline_mean_normalized_robinson_foulds"] == (
        0.533333333333333
    )
    assert payload["metrics"]["top_ranked_mean_normalized_robinson_foulds"] == 0.0
    assert payload["data"]["rows"][0]["taxon"] == "E"
    assert payload["outputs"] == [
        str(output_dir / "rogue-taxon-ranking.tsv"),
        str(output_dir / "baseline-consensus.nwk"),
        str(output_dir / "best-rogue-taxon-consensus.nwk"),
    ]
    assert (output_dir / "baseline-consensus.nwk").read_text(
        encoding="utf-8"
    ) == "((A:1,B:0.1)80:0.1,(C:0.1,D:0.1)80:0.1,E:0.1);\n"
    assert (output_dir / "best-rogue-taxon-consensus.nwk").read_text(
        encoding="utf-8"
    ) == "((A:1.02,B:0.1)100:0.12,(C:0.1,D:0.12)100:0.12);\n"
    assert (
        (output_dir / "rogue-taxon-ranking.tsv")
        .read_text(encoding="utf-8")
        .startswith(
            "rank\ttaxon\tmean_terminal_branch_length\tbaseline_consensus_resolution\t"
        )
    )
