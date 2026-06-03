from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_dating_relaxed_rate_summary_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "relaxed-rate-summary"

    exit_code = main(
        [
            "phylo",
            "dating",
            "relaxed-rate-summary",
            str(fixture("trees", "relaxed_rate_summary_substitution_tree_4_taxa.nwk")),
            str(fixture("trees", "relaxed_rate_summary_dated_tree_4_taxa.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["branch_count"] == 6
    assert payload["metrics"]["outlier_threshold"] == pytest.approx(2.0, abs=1e-12)
    assert payload["metrics"]["mean_branch_rate"] == pytest.approx(0.2, abs=1e-12)
    assert payload["metrics"]["minimum_branch_rate"] == pytest.approx(0.1, abs=1e-12)
    assert payload["metrics"]["maximum_branch_rate"] == pytest.approx(0.6, abs=1e-12)
    assert payload["metrics"]["outlier_count"] == 1
    assert payload["data"]["outlier_rows"][0]["child_name"] == "A"
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "branch_rates.tsv",
        "outliers.tsv",
        "run.json",
        "summary.tsv",
    ]
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "branch_rates.tsv").is_file()
    assert (out_dir / "outliers.tsv").is_file()
    assert (out_dir / "run.json").is_file()
