from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_dating_least_squares_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "least-squares-dating"

    exit_code = main(
        [
            "phylo",
            "dating",
            "least-squares",
            str(fixture("trees", "least_squares_dating_substitution_tree_4_taxa.nwk")),
            str(fixture("metadata", "least_squares_dating_tip_dates_4_taxa.tsv")),
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
    assert math.isclose(
        payload["metrics"]["estimated_clock_rate"],
        0.2500000000709406,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        payload["metrics"]["root_date"],
        1999.9999999931315,
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    assert math.isclose(
        payload["metrics"]["residual_sum_squares"],
        1.1053895766954377e-17,
        rel_tol=0.0,
        abs_tol=1e-24,
    )
    assert payload["metrics"]["exact_fit"] is True
    assert payload["metrics"]["converged"] is True
    assert (out_dir / "dated_tree.nwk").is_file()
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "node_dates.tsv").is_file()
    assert (out_dir / "branch_residuals.tsv").is_file()
    assert (out_dir / "run.json").is_file()
