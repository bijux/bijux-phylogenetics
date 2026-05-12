from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.biogeography import (
    TimeBinDefinition,
    summarize_time_stratified_geographic_transitions,
    write_time_stratified_branch_table,
    write_time_stratified_exclusion_table,
    write_time_stratified_transition_matrix_table,
    write_time_stratified_transition_summary_table,
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


def test_summarize_time_stratified_geographic_transitions_reports_interval_matrix() -> None:
    report = summarize_time_stratified_geographic_transitions(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="ard",
        time_bins=[
            TimeBinDefinition(label="early", start_depth=0.0, end_depth=0.1),
            TimeBinDefinition(label="late", start_depth=0.1, end_depth=0.3),
        ],
    )

    assert report.summary.time_bin_count == 2
    assert report.summary.matrix_row_count == len(report.matrix_rows)
    assert report.branch_rows
    assert report.matrix_rows


def test_summarize_time_stratified_geographic_transitions_validates_simulated_history(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "simulated-geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tnorth\nC\tisland\nD\tsouth\n",
        encoding="utf-8",
    )
    report = summarize_time_stratified_geographic_transitions(
        fixture("example_tree.nwk"),
        traits_path,
        trait="region",
        model="ard",
        time_bins=[
            TimeBinDefinition(label="early", start_depth=0.0, end_depth=0.1),
            TimeBinDefinition(label="late", start_depth=0.1, end_depth=0.3),
        ],
    )

    early_north_to_south = next(
        row
        for row in report.matrix_rows
        if row.time_bin_label == "early"
        and row.source_region == "north"
        and row.target_region == "south"
    )
    late_north_to_south = next(
        row
        for row in report.matrix_rows
        if row.time_bin_label == "late"
        and row.source_region == "north"
        and row.target_region == "south"
    )
    early_south_to_island = next(
        row
        for row in report.matrix_rows
        if row.time_bin_label == "early"
        and row.source_region == "south"
        and row.target_region == "island"
    )
    late_south_to_island = next(
        row
        for row in report.matrix_rows
        if row.time_bin_label == "late"
        and row.source_region == "south"
        and row.target_region == "island"
    )

    assert early_north_to_south.allocated_transition_weight > 0.0
    assert late_north_to_south.allocated_transition_weight == 0.0
    assert early_south_to_island.allocated_transition_weight == 0.0
    assert late_south_to_island.allocated_transition_weight > 0.0


def test_write_time_stratified_geographic_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_time_stratified_geographic_transitions(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="er",
        time_bins=[
            TimeBinDefinition(label="early", start_depth=0.0, end_depth=0.1),
            TimeBinDefinition(label="late", start_depth=0.1, end_depth=0.3),
        ],
    )

    summary_path = write_time_stratified_transition_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    matrix_path = write_time_stratified_transition_matrix_table(
        tmp_path / "matrix.tsv",
        report,
    )
    branches_path = write_time_stratified_branch_table(
        tmp_path / "branches.tsv",
        report,
    )
    exclusions_path = write_time_stratified_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "tree_depth" in summary_path.read_text(encoding="utf-8")
    assert "time_stratified_rate" in matrix_path.read_text(encoding="utf-8")
    assert "allocated_transition_weight" in branches_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_summarize_time_stratified_geographic_transitions_rejects_overlapping_bins() -> None:
    with pytest.raises(ValueError, match="must not overlap"):
        summarize_time_stratified_geographic_transitions(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model="er",
            time_bins=[
                TimeBinDefinition(label="one", start_depth=0.0, end_depth=0.2),
                TimeBinDefinition(label="two", start_depth=0.1, end_depth=0.3),
            ],
        )


def test_summarize_time_stratified_geographic_transitions_warns_on_partial_coverage() -> None:
    report = summarize_time_stratified_geographic_transitions(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="er",
        time_bins=[
            TimeBinDefinition(label="shallow", start_depth=0.0, end_depth=0.2),
        ],
    )

    assert report.summary.warning_count >= 1
    assert any(
        "do not cover the full tree depth" in warning for warning in report.warnings
    )
