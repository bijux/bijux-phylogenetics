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


@pytest.mark.slow
def test_host_association_switches_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    branches_path = tmp_path / "branches.tsv"
    counts_path = tmp_path / "counts.tsv"
    fits_path = tmp_path / "fits.tsv"
    unsupported_path = tmp_path / "unsupported.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "host-association",
            "switches",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_host.tsv")),
            "--trait",
            "host",
            "--model",
            "ard",
            "--constraints",
            str(fixture("example_host_transition_constraints.tsv")),
            "--summary-out",
            str(summary_path),
            "--nodes-out",
            str(nodes_path),
            "--branches-out",
            str(branches_path),
            "--counts-out",
            str(counts_path),
            "--fits-out",
            str(fits_path),
            "--unsupported-out",
            str(unsupported_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["analysis_constraint_mode"] == "constrained"
    assert payload["metrics"]["observed_host_count"] == 3
    assert "preferred_constraint" in summary_path.read_text(encoding="utf-8")
    assert "host_probabilities" in nodes_path.read_text(encoding="utf-8")
    assert "certainty_class" in branches_path.read_text(encoding="utf-8")
    assert "certain_switch_count" in counts_path.read_text(encoding="utf-8")
    assert "constraint_mode" in fits_path.read_text(encoding="utf-8")
    assert "claim_resolved" in unsupported_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
