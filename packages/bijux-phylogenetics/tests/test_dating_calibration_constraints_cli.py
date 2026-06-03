from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_dating_calibration_constraints_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "calibration-constraints"

    exit_code = main(
        [
            "phylo",
            "dating",
            "calibration-constraints",
            str(
                fixture(
                    "trees",
                    "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
                )
            ),
            str(
                fixture(
                    "metadata",
                    "dating_calibration_constraints_contradictory_5_taxa.tsv",
                )
            ),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["calibration_count"] == 3
    assert payload["metrics"]["valid_calibration_count"] == 3
    assert payload["metrics"]["contradictory_calibration_count"] == 3
    assert payload["metrics"]["contradictory_node_count"] == 3
    assert payload["metrics"]["feasible"] is False
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "constraints.tsv").is_file()
    assert (out_dir / "node_windows.tsv").is_file()
    assert (out_dir / "issues.tsv").is_file()
    assert (out_dir / "run.json").is_file()
