from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.dating import (
    summarize_relaxed_rate_branches_from_paths,
    write_relaxed_rate_branch_summary_artifacts,
    write_relaxed_rate_branch_summary_tsv,
    write_relaxed_rate_branch_table,
    write_relaxed_rate_outliers_tsv,
    write_relaxed_rate_run_json,
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
    return summarize_relaxed_rate_branches_from_paths(
        fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk"),
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk"),
    )


def test_write_relaxed_rate_branch_summary_tsv_writes_expected_row(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "summary.tsv"

    write_relaxed_rate_branch_summary_tsv(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [rows[0]]
    row = rows[0]
    assert row["substitution_tree_path"] == str(
        fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk")
    )
    assert row["dated_tree_path"] == str(
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")
    )
    assert row["tip_count"] == "4"
    assert row["internal_node_count"] == "3"
    assert row["branch_count"] == "6"
    assert float(row["mean_branch_rate"]) == pytest.approx(0.2, abs=1e-12)
    assert float(row["maximum_branch_rate"]) == pytest.approx(0.6, abs=1e-12)
    assert row["outlier_count"] == "1"
    assert "taxon:A" in row["outlier_branch_ids"]


def test_write_relaxed_rate_branch_table_writes_expected_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "branch_rates.tsv"

    write_relaxed_rate_branch_table(output_path, report)

    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 6
    rows_by_descendant_taxa = {
        tuple(row["descendant_taxa"].split("|")): row for row in rows
    }
    assert float(rows_by_descendant_taxa[("A",)]["branch_rate"]) == pytest.approx(
        0.6,
        abs=1e-12,
    )
    assert float(rows_by_descendant_taxa[("A",)]["rate_z_score"]) == pytest.approx(
        2.1908902300206643,
        abs=1e-14,
    )
    assert rows_by_descendant_taxa[("A",)]["outlier"] == "true"
    assert rows_by_descendant_taxa[("A", "B")]["child_name"] == ""


def test_write_relaxed_rate_outliers_tsv_writes_ranked_outlier_rows(
    tmp_path: Path,
) -> None:
    report = build_report()
    output_path = tmp_path / "outliers.tsv"

    write_relaxed_rate_outliers_tsv(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines == [
        "rank\tbranch_id\tchild_name\tdescendant_taxa\tsubstitution_branch_length\tdated_time_duration\tbranch_rate\trate_z_score",
        "1\troot:clade:A|B|C|D/clade:A|B|C/clade:A|B/taxon:A\tA\tA\t1.2\t2\t0.6\t2.19089023002066",
    ]


def test_write_relaxed_rate_run_json_serializes_report_fields(tmp_path: Path) -> None:
    report = build_report()
    output_path = tmp_path / "run.json"

    write_relaxed_rate_run_json(output_path, report)

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["tip_count"] == 4
    assert payload["internal_node_count"] == 3
    assert payload["branch_count"] == 6
    assert payload["outlier_count"] == 1
    assert len(payload["branch_rows"]) == 6
    assert payload["outlier_rows"][0]["child_name"] == "A"


def test_write_relaxed_rate_branch_summary_artifacts_materializes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = build_report()

    outputs = write_relaxed_rate_branch_summary_artifacts(tmp_path, report)

    assert sorted(outputs) == [
        "branch_rates_path",
        "outliers_path",
        "run_json_path",
        "summary_path",
    ]
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith("substitution_tree_path\tdated_tree_path\t")
    )
    assert (
        outputs["branch_rates_path"]
        .read_text(encoding="utf-8")
        .startswith("branch_id\tchild_name\tdescendant_taxa\t")
    )
    assert (
        outputs["outliers_path"]
        .read_text(encoding="utf-8")
        .startswith("rank\tbranch_id\tchild_name\tdescendant_taxa\t")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["outlier_rows"][0]["child_name"] == "A"
