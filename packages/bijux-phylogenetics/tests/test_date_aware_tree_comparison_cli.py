from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_compare_clade_ages_reports_date_aware_metrics_and_table(
    tmp_path: Path,
    capsys,
) -> None:
    table_path = tmp_path / "clade-ages.tsv"

    exit_code = main(
        [
            "compare",
            "clade-ages",
            str(fixture("strict_clock_time_tree_4_taxa.nwk")),
            str(fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")),
            "--out",
            str(table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["matched_clades"] == 3
    assert payload["metrics"]["age_rmse"] == pytest.approx(
        3.16227766016838,
        abs=1e-12,
    )
    assert payload["metrics"]["unstable_clades"] == 1
    assert payload["metrics"]["topology_equal"] is True
    assert payload["metrics"]["robinson_foulds_distance"] == 0
    assert payload["data"]["clade_rows"][2]["clade_id"] == "A|B|C|D"
    assert payload["data"]["clade_rows"][2]["unstable_age"] is True
    assert payload["outputs"] == [str(table_path)]
    assert table_path.read_text(encoding="utf-8").startswith(
        "clade_id\tnode_kind\ttaxon_count\tdescendant_taxa\tleft_age\t"
    )
