from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

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


def test_cli_distance_build_tree_supports_wpgma(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")),
            "--method",
            "wpgma",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.75,C:3.5)Inner3:0.5625,B:4.0625)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "wpgma"


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
    model_rows = payload["data"]["rows"]
    assert payload["metrics"]["model_count"] == len(model_rows)
    assert {row["model"] for row in model_rows} >= {
        "p-distance",
        "jukes-cantor",
        "kimura-2-parameter",
        "felsenstein-81",
        "tamura-nei-93",
    }


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


def test_cli_alignment_build_tree_reports_explicit_bionj_exclusion_json(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "build-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "bionj",
            "--out",
            "artifacts/distance-tree.nwk",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "unsupported_distance_tree_method_error"
    assert payload["errors"][0]["message"].startswith(
        "BIONJ is explicitly out of scope"
    )
    assert payload["errors"][0]["details"]["reference_surface"] == "ape::bionj"


def test_cli_distance_build_tree_reports_explicit_bionj_exclusion(capsys) -> None:
    with pytest.raises(SystemExit) as error:
        main(
            [
                "distance",
                "build-tree",
                str(fixture("example_distance_matrix.tsv")),
                "--method",
                "bionj",
                "--out",
                "artifacts/imported-distance-tree.nwk",
            ]
        )
    error_text = capsys.readouterr().err
    assert error.value.code == 2
    assert "unsupported_distance_tree_method_error" in error_text
    assert "ape::bionj" in error_text


def test_cli_distance_reference_json_output(capsys) -> None:
    exit_code = main(["distance", "reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    observations = payload["data"]["observations"]
    assert payload["metrics"]["case_count"] == len(observations)
    assert {row["case"] for row in observations} >= {
        "dna-p-distance",
        "dna-jukes-cantor",
        "dna-kimura-2-parameter",
        "dna-felsenstein-81",
        "dna-tamura-nei-93",
        "protein-p-distance",
        "ambiguity-partial-match",
        "ambiguity-strict-mismatch",
        "ambiguity-report-only",
    }


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
