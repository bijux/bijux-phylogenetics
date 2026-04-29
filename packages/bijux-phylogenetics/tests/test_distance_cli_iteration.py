from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_cli_alignment_distance_quality_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-quality",
            str(fixture("example_alignment_distance_saturated.fasta")),
            "--model",
            "jukes-cantor",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["decision"] == "risky"
    assert payload["metrics"]["saturated_pair_count"] > 0


def test_cli_alignment_bootstrap_tree_writes_outputs(tmp_path: Path, capsys) -> None:
    support_path = tmp_path / "support.tsv"
    tree_set_path = tmp_path / "bootstrap.trees"
    exit_code = main(
        [
            "alignment",
            "bootstrap-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "5",
            "--support-out",
            str(support_path),
            "--tree-set-out",
            str(tree_set_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert support_path.exists()
    assert tree_set_path.exists()
    assert payload["metrics"]["replicate_count"] == 5


def test_cli_distance_reference_json_output(capsys) -> None:
    exit_code = main(["distance", "reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    assert payload["metrics"]["case_count"] == 3
