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


def test_biogeography_model_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    rates_path = tmp_path / "rates.tsv"
    events_path = tmp_path / "events.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--summary-out",
            str(summary_path),
            "--nodes-out",
            str(nodes_path),
            "--rates-out",
            str(rates_path),
            "--events-out",
            str(events_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["transition_rate_row_count"] > 0
    assert payload["metrics"]["changed_branch_count"] >= 0
    assert "root_region" in summary_path.read_text(encoding="utf-8")
    assert "most_likely_region" in nodes_path.read_text(encoding="utf-8")
    assert "source_region" in rates_path.read_text(encoding="utf-8")
    assert "strongly_supported" in events_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_model_cli_accepts_region_vocabulary(capsys) -> None:
    exit_code = main(
        [
            "biogeography",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "sym",
            "--allowed-regions",
            "north,south,island",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "sym"
    assert payload["metrics"]["observed_region_count"] == 3
