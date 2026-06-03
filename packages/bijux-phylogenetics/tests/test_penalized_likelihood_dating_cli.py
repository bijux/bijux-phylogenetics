from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
DATE_ABS_TOLERANCE = 1e-5
SCORE_REL_TOLERANCE = 1e-6
SCORE_ABS_TOLERANCE = 1e-12


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_dating_penalized_likelihood_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "penalized-likelihood-dating"

    exit_code = main(
        [
            "phylo",
            "dating",
            "penalized-likelihood",
            str(
                fixture(
                    "trees",
                    "penalized_likelihood_dating_substitution_tree_4_taxa.nwk",
                )
            ),
            str(
                fixture(
                    "metadata",
                    "penalized_likelihood_dating_tip_dates_4_taxa.tsv",
                )
            ),
            "--smoothing-parameter",
            "0.01",
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
        payload["metrics"]["smoothing_parameter"],
        0.01,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        payload["metrics"]["data_score"],
        2.9289583718815336e-06,
        rel_tol=SCORE_REL_TOLERANCE,
        abs_tol=SCORE_ABS_TOLERANCE,
    )
    assert math.isclose(
        payload["metrics"]["penalty_score"],
        0.00013482416705344673,
        rel_tol=SCORE_REL_TOLERANCE,
        abs_tol=SCORE_ABS_TOLERANCE,
    )
    assert math.isclose(
        payload["metrics"]["total_score"],
        0.00013775312542532825,
        rel_tol=SCORE_REL_TOLERANCE,
        abs_tol=SCORE_ABS_TOLERANCE,
    )
    assert math.isclose(
        payload["metrics"]["root_date"],
        1985.738765803845,
        rel_tol=0.0,
        abs_tol=DATE_ABS_TOLERANCE,
    )
    assert payload["metrics"]["optimization_pass_count"] == 5
    assert payload["metrics"]["function_evaluation_count"] == 771
    assert payload["metrics"]["converged"] is True
    assert (out_dir / "dated_tree.nwk").is_file()
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "node_dates.tsv").is_file()
    assert (out_dir / "branch_rates.tsv").is_file()
    assert (out_dir / "run.json").is_file()
