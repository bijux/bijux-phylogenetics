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


def test_cli_alignment_distance_assumptions_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-assumptions",
            str(fixture("example_alignment_distance.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["ultrametric_compatible"] is False
    assert payload["metrics"]["upgma_violation_count"] > 0


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


def test_cli_alignment_distance_support_summary_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-support-summary",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "3",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 5
    assert payload["metrics"]["clade_count"] > 0


def test_cli_alignment_distance_models_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-models",
            str(fixture("example_alignment_distance.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["model_count"] == 3


def test_cli_alignment_distance_gap_sensitivity_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-gap-sensitivity",
            str(fixture("example_alignment_distance_gaps.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["changed_pair_count"] > 0


def test_cli_alignment_distance_maturity_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-maturity",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "3",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["check_count"] > 0
    assert payload["metrics"]["decision"] in {
        "production_candidate",
        "validated_with_limits",
    }


def test_cli_distance_reference_json_output(capsys) -> None:
    exit_code = main(["distance", "reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    assert payload["metrics"]["case_count"] == 9


def test_cli_distance_assumptions_json_output(capsys) -> None:
    exit_code = main(
        [
            "distance",
            "assumptions",
            str(fixture("example_distance_matrix_nonultrametric.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["ultrametric_compatible"] is False
    assert payload["metrics"]["upgma_violation_count"] > 0


def test_cli_distance_quality_json_output(capsys) -> None:
    exit_code = main(
        ["distance", "quality", str(fixture("example_distance_matrix.tsv")), "--json"]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["saturation_audit_scale"] == "unit-interval-like"
    assert payload["metrics"]["low_information_pair_count"] == 3
