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


def test_diversification_ltt_sampling_and_estimate_cli_outputs_json_and_tables(
    tmp_path: Path, capsys
) -> None:
    ltt_path = tmp_path / "ltt.tsv"
    gamma_path = tmp_path / "gamma-statistic.tsv"

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

    gamma_exit = main(
        [
            "diversification",
            "gamma-stat",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--out",
            str(gamma_path),
            "--json",
        ]
    )
    gamma_payload = json.loads(capsys.readouterr().out)
    assert gamma_exit == 0
    assert gamma_payload["metrics"]["tip_count"] == 4
    assert gamma_payload["metrics"]["sampling_fraction"] == 0.75
    assert gamma_payload["metrics"]["gamma_statistic"] < 0.0
    assert any(
        "assumes complete taxon sampling" in warning
        for warning in gamma_payload["warnings"]
    )
    assert "gamma_statistic" in gamma_path.read_text(encoding="utf-8")


def test_diversification_compare_clades_trait_and_report_cli_write_outputs(
    tmp_path: Path, capsys
) -> None:
    clades_path = tmp_path / "clades.tsv"
    trait_path = tmp_path / "trait.tsv"
    report_path = tmp_path / "diversification-report.html"
    methods_summary_path = tmp_path / "diversification-methods-summary.md"

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
            "--methods-summary-out",
            str(methods_summary_path),
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["report_kind"] == "diversification"
    assert report_payload["metrics"]["methods_summary_warning_count"] >= 0
    assert report_payload["metrics"]["better_model"] == "yule"
    assert "diversification-model-comparison" in report_path.read_text(encoding="utf-8")
    assert methods_summary_path.exists()


def test_diversification_methods_summary_cli_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "diversification-methods-summary.md"

    exit_code = main(
        [
            "diversification",
            "methods-summary",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions_incomplete.tsv")),
            "--traits",
            str(fixture("example_traits_diversification_polyphyletic.tsv")),
            "--trait",
            "habitat",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["warning_count"] >= 1
    assert payload["metrics"]["better_model"] == "yule"
    assert payload["metrics"]["sampling_metadata_complete"] is False
    assert payload["metrics"]["clade_observation_count"] == 3
    assert payload["metrics"]["trait_state_count"] == 2
    assert output_path.exists()


def test_diversification_package_cli_writes_publication_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "diversification-figure-package"

    exit_code = main(
        [
            "diversification",
            "package",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["sampling_metadata_complete"] is True
    assert payload["metrics"]["plotted_ltt_point_count"] == 4
    assert payload["metrics"]["plotted_clade_count"] == 3
    assert payload["metrics"]["highlighted_outlier_count"] == 2
    assert payload["metrics"]["plotted_model_count"] == 2
    assert payload["metrics"]["better_model"] in {"yule", "birth-death"}
    assert payload["metrics"]["methods_summary_warning_count"] >= 0
    assert len(payload["outputs"]) == 12
    assert (output_dir / "lineage-through-time.svg").exists()
    assert (output_dir / "clade-diversification-outliers.svg").exists()
    assert (output_dir / "diversification-model-comparison.svg").exists()
    assert (output_dir / "diversification-methods-summary.md").exists()
    assert (output_dir / "diversification-figure-review.html").exists()
    assert (output_dir / "diversification-figure-package.manifest.json").exists()
    assert (output_dir / "figure-reproducibility.manifest.json").exists()


def test_diversification_medusa_cli_reports_explicit_exclusion(capsys) -> None:
    exit_code = main(
        [
            "diversification",
            "medusa",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "diversification_medusa_explicitly_excluded"
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "geiger_medusa_explicitly_excluded_this_round"
    )
    assert payload["errors"][0]["details"]["sampling_metadata_complete"] is True
    assert (
        "descriptive clade diversification outlier scan"
        in payload["errors"][0]["details"]["supported_surfaces"]
    )
    assert (
        "stepwise branch-specific rate-shift search"
        in payload["errors"][0]["details"]["missing_surfaces"]
    )


def test_diversification_bd_ms_cli_reports_explicit_exclusion(capsys) -> None:
    exit_code = main(
        [
            "diversification",
            "bd-ms",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert (
        payload["errors"][0]["code"]
        == "diversification_birth_death_explicitly_excluded"
    )
    assert (
        payload["errors"][0]["details"]["failure_reason"]
        == "geiger_birth_death_explicitly_excluded_this_round"
    )
    assert payload["errors"][0]["details"]["geiger_reference_surface"] == (
        "geiger::bd.ms"
    )
    assert payload["errors"][0]["details"]["geiger_reference_arguments"] == [
        "phy",
        "time",
        "n",
        "missing",
        "crown",
        "epsilon",
    ]
    assert payload["errors"][0]["details"]["sampling_metadata_complete"] is True
