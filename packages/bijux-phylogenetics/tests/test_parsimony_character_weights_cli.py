from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


@pytest.mark.parametrize(
    ("command", "expected_metric", "expected_total"),
    [
        (
            [
                "phylo",
                "parsimony",
                "fitch",
                str(fixture("fitch_tree.nwk")),
                str(fixture("fitch_binary_matrix.tsv")),
                "--character-weights",
                str(fixture("fitch_character_weights.tsv")),
            ],
            "total_weighted_score",
            3.5,
        ),
        (
            [
                "phylo",
                "parsimony",
                "wagner",
                str(fixture("fitch_tree.nwk")),
                str(fixture("wagner_ordinal_matrix.tsv")),
                "--character-weights",
                str(fixture("wagner_character_weights.tsv")),
            ],
            "total_weighted_score",
            7.5,
        ),
        (
            [
                "phylo",
                "parsimony",
                "sankoff",
                str(fixture("sankoff_tree_5_taxa.nwk")),
                str(fixture("sankoff_character_matrix.tsv")),
                str(fixture("sankoff_cost_matrix.tsv")),
                "--character-weights",
                str(fixture("sankoff_character_weights.tsv")),
            ],
            "total_weighted_score",
            7.5,
        ),
        (
            [
                "phylo",
                "parsimony",
                "dollo",
                str(fixture("dollo_tree_5_taxa.nwk")),
                str(fixture("dollo_binary_matrix.tsv")),
                "--character-weights",
                str(fixture("dollo_character_weights.tsv")),
            ],
            "total_weighted_score",
            9.5,
        ),
        (
            [
                "phylo",
                "parsimony",
                "camin-sokal",
                str(fixture("camin_sokal_tree_5_taxa.nwk")),
                str(fixture("camin_sokal_binary_matrix.tsv")),
                "--character-weights",
                str(fixture("camin_sokal_character_weights.tsv")),
            ],
            "total_weighted_score",
            5.5,
        ),
        (
            [
                "phylo",
                "parsimony",
                "acctran",
                str(fixture("acctran_tree_5_taxa.nwk")),
                str(fixture("acctran_ambiguous_matrix.tsv")),
                "--character-weights",
                str(fixture("acctran_character_weights.tsv")),
            ],
            "total_weighted_score",
            5.0,
        ),
        (
            [
                "phylo",
                "parsimony",
                "deltran",
                str(fixture("acctran_tree_5_taxa.nwk")),
                str(fixture("acctran_ambiguous_matrix.tsv")),
                "--character-weights",
                str(fixture("acctran_character_weights.tsv")),
            ],
            "total_weighted_score",
            5.0,
        ),
    ],
)
def test_phylo_parsimony_cli_surfaces_forward_character_weights(
    command: list[str],
    expected_metric: str,
    expected_total: float,
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main([*command, "--out-dir", str(tmp_path), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"][expected_metric] == expected_total


def test_phylo_parsimony_cli_reports_negative_weight_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "fitch",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_binary_matrix.tsv")),
            "--character-weights",
            str(fixture("parsimony_negative_character_weight.tsv")),
            "--out-dir",
            str(tmp_path / "fitch-negative"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_character_weight_negative_value"
