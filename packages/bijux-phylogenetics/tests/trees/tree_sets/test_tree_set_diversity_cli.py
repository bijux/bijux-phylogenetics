from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_tree_set_diversity_reports_rf_distribution_metrics(
    tmp_path: Path,
    capsys,
) -> None:
    output_path = tmp_path / "rf-distribution.tsv"

    exit_code = main(
        [
            "tree-set",
            "diversity",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["rooted_topology_count"] == 2
    assert payload["metrics"]["pair_count"] == 3
    assert payload["metrics"]["rf_bucket_count"] >= 1
    assert payload["metrics"]["runtime_seconds"] >= 0.0
    assert payload["metrics"]["peak_memory_bytes"] >= 0
    assert payload["metrics"]["skipped_malformed_tree_count"] == 0
    assert output_path.read_text(encoding="utf-8").startswith(
        "robinson_foulds_distance\tnormalized_robinson_foulds\tpair_count\tfrequency\n"
    )
