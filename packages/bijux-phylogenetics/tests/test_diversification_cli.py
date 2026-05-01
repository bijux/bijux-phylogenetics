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


def test_diversification_ltt_sampling_and_estimate_cli_outputs_json_and_tables(
    tmp_path: Path, capsys
) -> None:
    ltt_path = tmp_path / "ltt.tsv"

    ltt_exit = main(
        [
            "diversification",
            "ltt",
            str(fixture("example_tree.nwk")),
            "--out",
            str(ltt_path),
            "--json",
        ]
    )
    ltt_payload = json.loads(capsys.readouterr().out)
    assert ltt_exit == 0
    assert ltt_payload["metrics"]["point_count"] == 4
    assert "lineage_count" in ltt_path.read_text(encoding="utf-8")

    sampling_exit = main(
        [
            "diversification",
            "sampling",
            str(fixture("example_tree.nwk")),
            str(fixture("example_sampling_fractions_incomplete.tsv")),
            "--json",
        ]
    )
    sampling_payload = json.loads(capsys.readouterr().out)
    assert sampling_exit == 0
    assert sampling_payload["metrics"]["missing_taxon_count"] == 1
    assert sampling_payload["metrics"]["invalid_row_count"] == 2

    estimate_exit = main(
        [
            "diversification",
            "estimate",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--model",
            "birth-death",
            "--json",
        ]
    )
    estimate_payload = json.loads(capsys.readouterr().out)
    assert estimate_exit == 0
    assert estimate_payload["metrics"]["model"] == "birth-death"
    assert estimate_payload["metrics"]["sampling_fraction"] == 0.75


def test_diversification_compare_clades_trait_and_report_cli_write_outputs(
    tmp_path: Path, capsys
) -> None:
    clades_path = tmp_path / "clades.tsv"
    trait_path = tmp_path / "trait.tsv"
    report_path = tmp_path / "diversification-report.html"

    compare_exit = main(
        [
            "diversification",
            "compare-models",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--json",
        ]
    )
    compare_payload = json.loads(capsys.readouterr().out)
    assert compare_exit == 0
    assert compare_payload["metrics"]["model_count"] == 2

    clades_exit = main(
        [
            "diversification",
            "clades",
            str(fixture("example_tree.nwk")),
            "--out",
            str(clades_path),
            "--json",
        ]
    )
    clades_payload = json.loads(capsys.readouterr().out)
    assert clades_exit == 0
    assert clades_payload["metrics"]["high_clade_count"] == 1
    assert "classification" in clades_path.read_text(encoding="utf-8")

    trait_exit = main(
        [
            "diversification",
            "trait-dependent",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_diversification.tsv")),
            "--trait",
            "habitat",
            "--out",
            str(trait_path),
            "--json",
        ]
    )
    trait_payload = json.loads(capsys.readouterr().out)
    assert trait_exit == 0
    assert trait_payload["metrics"]["state_count"] == 2
    assert "monophyletic" in trait_path.read_text(encoding="utf-8")

    report_exit = main(
        [
            "diversification",
            "report",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--traits",
            str(fixture("example_traits_diversification.tsv")),
            "--trait",
            "habitat",
            "--out",
            str(report_path),
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["report_kind"] == "diversification"
    assert "diversification-model-comparison" in report_path.read_text(encoding="utf-8")
