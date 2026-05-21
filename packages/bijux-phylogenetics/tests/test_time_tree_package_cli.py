from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_cli_report_time_tree_package_writes_publication_review_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "time-tree-package"
    exit_code = main(
        [
            "report",
            "time-tree-package",
            str(metadata_fixture("beast2_strict_yule_posterior.trees")),
            "--source-format",
            "beast",
            "--burnin-fraction",
            "0.25",
            "--metadata",
            str(metadata_fixture("example_metadata.tsv")),
            "--label-column",
            "species",
            "--tip-dates",
            str(metadata_fixture("example_tip_dates.tsv")),
            "--title",
            "Rabies time tree",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["retained_tree_count"] > 0
    assert payload["metrics"]["rendered_interval_count"] >= 1
    assert (
        payload["metrics"]["rendered_interval_count"]
        == payload["metrics"]["expected_interval_count"]
    )
    assert payload["metrics"]["ultrametric"] is True
    assert payload["metrics"]["readiness_decision"] == "ready"
    assert payload["data"]["audit"]["publication_ready"] is True
    assert payload["data"]["manifest_path"].endswith("time-tree-package.manifest.json")
    assert (output_dir / "time-tree-review.html").exists()
    assert (output_dir / "node-age-intervals.tsv").exists()
    assert (output_dir / "figure-caption.md").exists()
    assert (output_dir / "figure-reproducibility.manifest.json").exists()
