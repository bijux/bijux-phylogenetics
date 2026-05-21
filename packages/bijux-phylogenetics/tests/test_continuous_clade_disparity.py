from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.comparative.disparity import (
    summarize_continuous_clade_disparity,
    write_continuous_clade_disparity_summary_table,
    write_continuous_clade_disparity_table,
    write_disparity_through_time_exclusion_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_summarize_continuous_clade_disparity_reports_method_formula_and_ranges() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"),
        trait_columns=["ou_truth", "early_burst_truth"],
    )

    assert report.distance_metric == "avg-squared-euclidean"
    assert "squared Euclidean distance" in report.method_formula
    assert report.root_disparity == report.clade_rows[0].disparity
    assert report.minimum_clade_disparity <= report.root_disparity
    assert report.maximum_clade_disparity >= report.root_disparity


def test_summarize_continuous_clade_disparity_matches_known_simple_clade_partitions() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait_columns=["response"],
    )

    assert [row.ape_node_id for row in report.clade_rows] == [5, 6, 7]
    assert [row.descendant_taxa for row in report.clade_rows] == [
        ["A", "B", "C", "D"],
        ["A", "B"],
        ["C", "D"],
    ]
    observed = [row.disparity for row in report.clade_rows]
    expected = [2.16666666666667, 2.25, 2.25]
    for value, target in zip(observed, expected, strict=True):
        assert math.isclose(value, target, rel_tol=1e-12, abs_tol=1e-12)


def test_summarize_continuous_clade_disparity_supports_known_multivariate_clades() -> (
    None
):
    report = summarize_continuous_clade_disparity(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait_columns=["response", "predictor_one"],
    )

    observed = [row.disparity for row in report.clade_rows]
    expected = [5.5, 3.25, 3.25]
    for value, target in zip(observed, expected, strict=True):
        assert math.isclose(value, target, rel_tol=1e-12, abs_tol=1e-12)


def test_continuous_clade_disparity_writers_emit_summary_and_long_tables(
    tmp_path: Path,
) -> None:
    report = summarize_continuous_clade_disparity(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait_columns=["response", "predictor_one"],
    )
    summary_path = tmp_path / "disparity-summary.tsv"
    clades_path = tmp_path / "disparity-clades.tsv"
    excluded_path = tmp_path / "disparity-excluded.tsv"

    write_continuous_clade_disparity_summary_table(summary_path, report)
    write_continuous_clade_disparity_table(clades_path, report)
    write_disparity_through_time_exclusion_table(excluded_path, report)

    with summary_path.open(encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    with clades_path.open(encoding="utf-8", newline="") as handle:
        clade_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(summary_rows) == 1
    assert summary_rows[0]["method_formula"]
    assert len(clade_rows) == 3
    assert excluded_path.exists()
