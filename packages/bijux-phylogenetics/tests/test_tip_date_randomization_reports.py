from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.root_to_tip import (
    diagnose_tip_date_randomization,
    write_tip_date_randomization_artifacts,
    write_tip_date_randomization_permutations_tsv,
    write_tip_date_randomization_run_json,
    write_tip_date_randomization_summary_tsv,
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


def build_report():
    return diagnose_tip_date_randomization(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
        permutations=19,
        seed=17,
    )


def test_write_tip_date_randomization_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "summary.tsv"

    write_tip_date_randomization_summary_tsv(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "tree_path\tmetadata_path\tsource_format\ttaxon_column\tdate_column\t"
        "tip_count\tobserved_slope\tobserved_intercept\tobserved_r_squared\t"
        "observed_residual_mean_square\tobserved_outlier_count\toutlier_threshold\t"
        "permutations\tseed\tpermuted_r_squared_at_or_above_observed\t"
        "tip_date_randomization_p_value\tnull_distribution_minimum\t"
        "null_distribution_mean\tnull_distribution_maximum"
    )
    assert lines[1].endswith(
        "\t7\t2.39285714285714\t-1.82142857142857\t0.639094533029613\t"
        "18.1071428571429\t1\t2\t19\t17\t0\t0.05\t0.000142369020501021\t"
        "0.152267413979139\t0.547266514806378"
    )


def test_write_tip_date_randomization_permutations_tsv_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "permutations.tsv"

    write_tip_date_randomization_permutations_tsv(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "permutation_index\tpermuted_slope\tpermuted_intercept\t"
        "permuted_r_squared\tpermuted_residual_mean_square\tat_or_above_observed"
    )
    assert lines[1] == (
        "1\t0.607142857142857\t3.53571428571429\t0.041144646924829\t"
        "48.1071428571429\tfalse"
    )
    assert len(lines) == 20


def test_write_tip_date_randomization_run_json_serializes_report_fields(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_tip_date_randomization_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["tip_count"] == 7
    assert payload["permutations"] == 19
    assert payload["seed"] == 17
    assert payload["p_value"] == 0.05
    assert payload["permuted_r_squared_at_or_above_observed"] == 0
    assert payload["observed_regression"]["outliers"][0]["tip"] == "G"
    assert len(payload["permutation_rows"]) == 19


def test_write_tip_date_randomization_artifacts_emit_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    artifact_paths = write_tip_date_randomization_artifacts(tmp_path, report)

    assert sorted(artifact_paths) == [
        "observed_outliers",
        "observed_residuals",
        "permutations",
        "run_json",
        "summary",
    ]
    assert (
        (tmp_path / "summary.tsv")
        .read_text(encoding="utf-8")
        .startswith("tree_path\tmetadata_path\tsource_format\t")
    )
    assert (
        (tmp_path / "permutations.tsv")
        .read_text(encoding="utf-8")
        .startswith("permutation_index\tpermuted_slope\t")
    )
    assert (
        (tmp_path / "observed_residuals.tsv")
        .read_text(encoding="utf-8")
        .startswith("tip\tsampling_time\troot_to_tip_distance\t")
    )
    assert (
        (tmp_path / "observed_outliers.tsv")
        .read_text(encoding="utf-8")
        .startswith("rank\ttip\tsampling_time\troot_to_tip_distance\t")
    )
    payload = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert payload["observed_regression"]["outliers"][0]["tip"] == "G"
