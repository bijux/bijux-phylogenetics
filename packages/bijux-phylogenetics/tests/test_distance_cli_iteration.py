from __future__ import annotations

import json
from pathlib import Path

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


def test_cli_distance_build_tree_supports_single_linkage(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_single_linkage_chain.tsv")),
            "--method",
            "single-linkage",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((((A:0.5,B:0.5)Inner1:0.5,C:1)Inner2:0.5,D:1.5)Inner3:0.5,E:2)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "single-linkage"


def test_cli_distance_build_tree_supports_complete_linkage(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")),
            "--method",
            "complete-linkage",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "complete-linkage"


def test_cli_distance_minimum_evolution_writes_fitted_tree(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "minimum-evolution-tree.nwk"
    exit_code = main(
        [
            "distance",
            "minimum-evolution",
            str(fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv")),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:1,B:2):3,C:4,(D:5,E:6):7);\n"
    )
    assert payload["metrics"]["criterion"] == "minimum-evolution"
    assert payload["metrics"]["minimum_evolution_score"] == 28.0


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


def test_cli_alignment_build_tree_supports_bionj_json(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "alignment",
            "build-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "bionj",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;\n"
    )
    assert payload["metrics"]["method"] == "bionj"


def test_cli_distance_build_tree_supports_bionj(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_bionj_noisy.tsv")),
            "--method",
            "bionj",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;\n"
    )
    assert payload["metrics"]["method"] == "bionj"


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
