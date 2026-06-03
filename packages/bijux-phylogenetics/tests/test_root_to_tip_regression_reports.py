from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.root_to_tip import (
    diagnose_root_to_tip_regression,
    write_root_to_tip_regression_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_write_root_to_tip_regression_artifacts_emit_summary_residuals_and_outliers(
    tmp_path: Path,
) -> None:
    report = diagnose_root_to_tip_regression(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
    )

    artifact_paths = write_root_to_tip_regression_artifacts(tmp_path, report)

    assert sorted(artifact_paths) == ["outliers", "residuals", "run_json", "summary"]
    summary_lines = (tmp_path / "summary.tsv").read_text(encoding="utf-8").splitlines()
    residual_lines = (
        (tmp_path / "residuals.tsv").read_text(encoding="utf-8").splitlines()
    )
    outlier_lines = (tmp_path / "outliers.tsv").read_text(encoding="utf-8").splitlines()
    payload = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))

    assert summary_lines[0] == (
        "tree_path\tmetadata_path\tsource_format\ttaxon_column\tdate_column\t"
        "tip_count\tslope\tintercept\tr_squared\tresidual_mean_square\t"
        "outlier_threshold\toutlier_count\toutlier_tips\tsampling_time_min\t"
        "sampling_time_max\troot_to_tip_min\troot_to_tip_max"
    )
    assert summary_lines[1].endswith("\t1\tG\t0\t6\t0.5\t19.5")
    assert residual_lines[0] == (
        "tip\tsampling_time\troot_to_tip_distance\tfitted_distance\tresidual\t"
        "leverage\tstudentized_residual\toutlier"
    )
    assert residual_lines[-1].endswith("\t2.23606797749979\ttrue")
    assert outlier_lines == [
        "rank\ttip\tsampling_time\troot_to_tip_distance\tresidual\tstudentized_residual\tleverage",
        "1\tG\t6\t19.5\t6.96428571428572\t2.23606797749979\t0.464285714285714",
    ]
    assert payload["outliers"] == [
        {
            "leverage": 0.4642857142857143,
            "rank": 1,
            "residual": 6.964285714285715,
            "root_to_tip_distance": 19.5,
            "sampling_time": 6.0,
            "studentized_residual": 2.23606797749979,
            "tip": "G",
        }
    ]
