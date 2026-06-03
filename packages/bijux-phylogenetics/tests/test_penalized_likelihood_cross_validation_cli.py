from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
CV_ERROR_ABS_TOLERANCE = 5e-8


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_phylo_dating_penalized_likelihood_cross_validation_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "penalized-likelihood-cross-validation"

    exit_code = main(
        [
            "phylo",
            "dating",
            "penalized-likelihood-cross-validation",
            str(
                fixture(
                    "trees",
                    "penalized_likelihood_cross_validation_substitution_tree_5_taxa.nwk",
                )
            ),
            str(
                fixture(
                    "metadata",
                    "penalized_likelihood_cross_validation_tip_dates_5_taxa.tsv",
                )
            ),
            str(
                fixture(
                    "metadata",
                    "penalized_likelihood_cross_validation_calibrations_5_taxa.tsv",
                )
            ),
            "--smoothing-parameters",
            "0.01",
            "0.1",
            "1.0",
            "10.0",
            "100.0",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["usable_calibration_count"] == 4
    assert payload["metrics"]["candidate_count"] == 5
    assert math.isclose(
        payload["metrics"]["selected_smoothing_parameter"],
        0.01,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        payload["metrics"]["selected_root_mean_squared_error"],
        0.28564007040784534,
        rel_tol=0.0,
        abs_tol=CV_ERROR_ABS_TOLERANCE,
    )
    assert payload["metrics"]["final_converged"] is True
    assert (out_dir / "dated_tree.nwk").is_file()
    assert (out_dir / "summary.tsv").is_file()
    assert (out_dir / "candidate_scores.tsv").is_file()
    assert (out_dir / "prediction_errors.tsv").is_file()
    assert (out_dir / "calibrations.tsv").is_file()
    assert (out_dir / "node_dates.tsv").is_file()
    assert (out_dir / "branch_rates.tsv").is_file()
    assert (out_dir / "run.json").is_file()
