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


def test_ecological_niche_transitions_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    rates_path = tmp_path / "rates.tsv"
    branches_path = tmp_path / "branches.tsv"
    counts_path = tmp_path / "counts.tsv"
    clades_path = tmp_path / "clades.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "ecological-niche",
            "transitions",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_ecological_niche.tsv")),
            "--trait",
            "niche",
            "--model",
            "er",
            "--summary-out",
            str(summary_path),
            "--nodes-out",
            str(nodes_path),
            "--rates-out",
            str(rates_path),
            "--branches-out",
            str(branches_path),
            "--counts-out",
            str(counts_path),
            "--clades-out",
            str(clades_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "er"
    assert payload["metrics"]["observed_niche_count"] == 4
    assert payload["metrics"]["changed_branch_count"] >= 2
    assert "root_niche" in summary_path.read_text(encoding="utf-8")
    assert "most_likely_niche" in nodes_path.read_text(encoding="utf-8")
    assert "source_niche" in rates_path.read_text(encoding="utf-8")
    assert "certainty_class" in branches_path.read_text(encoding="utf-8")
    assert "certain_transition_count" in counts_path.read_text(encoding="utf-8")
    assert "shift_burden_score" in clades_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
