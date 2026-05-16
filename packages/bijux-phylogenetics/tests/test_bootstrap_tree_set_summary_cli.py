from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_tree_set_bootstrap_summary_writes_artifact_bundle(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "bootstrap-summary"
    exit_code = main(
        [
            "tree-set",
            "bootstrap-summary",
            str(fixture("example_tree_set_left.nwk")),
            "--out-dir",
            str(out_dir),
            "--prefix",
            "bootstrap-review",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["rooted_topology_count"] == 2
    assert payload["metrics"]["unstable_branch_count"] == 2
    assert payload["metrics"]["runtime_seconds"] >= 0.0
    assert payload["metrics"]["peak_memory_bytes"] >= 0
    assert payload["metrics"]["skipped_malformed_tree_count"] == 0
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "bootstrap-review.clade-frequencies.tsv",
        "bootstrap-review.consensus.nwk",
        "bootstrap-review.distance-matrix.tsv",
        "bootstrap-review.rf-distribution.tsv",
        "bootstrap-review.summary.tsv",
        "bootstrap-review.topology-clusters.tsv",
        "bootstrap-review.unstable-branches.tsv",
        "bootstrap-review.unstable-clades.tsv",
    ]
    assert out_dir.joinpath("bootstrap-review.unstable-branches.tsv").exists()
