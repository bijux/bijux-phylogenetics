from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.diversification import (
    compare_diversification_models,
    compute_diversification_gamma_statistic,
    compute_lineage_through_time_curve,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
    inspect_diversification_time_tree,
    render_diversification_report,
    run_trait_dependent_diversification_analysis,
    validate_time_tree_for_diversification,
    write_clade_diversification_table,
    write_diversification_gamma_statistic_table,
    write_lineage_through_time_table,
    write_trait_dependent_diversification_table,
)
from bijux_phylogenetics.errors import DiversificationAnalysisError, UnrootedTreeError

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


def test_validate_time_tree_for_diversification_reports_root_age() -> None:
    report = validate_time_tree_for_diversification(fixture("example_tree.nwk"))

    assert report.rooted is True
    assert report.ultrametric is True
    assert report.branch_length_status == "complete"
    assert report.tip_count == 4
    assert report.root_age == 0.3


def test_validate_time_tree_for_diversification_accepts_near_ultrametric_tree() -> None:
    report = validate_time_tree_for_diversification(
        fixture("example_tree_near_ultrametric.nwk")
    )

    assert report.rooted is True
    assert report.ultrametric is True
    assert report.root_age == pytest.approx(0.300000001, abs=1e-12)


def test_compute_lineage_through_time_curve_tracks_lineage_increases() -> None:
    report = compute_lineage_through_time_curve(fixture("example_tree.nwk"))

    assert [
        (point.time_before_present, point.lineage_count) for point in report.points
    ] == [
        (0.3, 2),
        (0.2, 3),
        (0.1, 4),
        (0.0, 4),
    ]


def test_compute_diversification_gamma_statistic_matches_governed_ultrametric_examples() -> (
    None
):
    balanced = compute_diversification_gamma_statistic(fixture("example_tree.nwk"))
    larger = compute_diversification_gamma_statistic(fixture("example_tree_eight_taxa.nwk"))
    zero_internal = compute_diversification_gamma_statistic(
        fixture("example_tree_ultrametric_zero_internal.nwk")
    )

    assert balanced.tip_count == 4
    assert balanced.bifurcating is True
    assert balanced.gamma_statistic == pytest.approx(-0.544331053951818, abs=1e-15)
    assert larger.gamma_statistic == pytest.approx(-1.414213562373095, abs=1e-14)
    assert zero_internal.gamma_statistic == pytest.approx(
        -0.979795897113271,
        abs=1e-15,
    )


def test_compute_diversification_gamma_statistic_warns_on_incomplete_sampling() -> None:
    report = compute_diversification_gamma_statistic(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
    )

    assert report.sampling_fraction == 0.75
    assert any("assumes complete taxon sampling" in warning for warning in report.warnings)


def test_compute_diversification_gamma_statistic_rejects_small_or_non_bifurcating_trees(
    tmp_path: Path,
) -> None:
    small_tree_path = tmp_path / "two-tip-ultrametric.nwk"
    small_tree_path.write_text("(A:0.3,B:0.3);", encoding="utf-8")
    singleton_tree_path = tmp_path / "singleton-root-ultrametric.nwk"
    singleton_tree_path.write_text(
        "(((A:0.1,B:0.1):0.2):0.0,C:0.3);",
        encoding="utf-8",
    )

    with pytest.raises(DiversificationAnalysisError) as small_error:
        compute_diversification_gamma_statistic(small_tree_path)
    assert (
        small_error.value.code
        == "diversification_gamma_statistic_requires_three_or_more_tips"
    )

    with pytest.raises(DiversificationAnalysisError) as bifurcating_error:
        compute_diversification_gamma_statistic(singleton_tree_path)
    assert (
        bifurcating_error.value.code
        == "diversification_gamma_statistic_requires_bifurcating_tree"
    )


def test_write_lineage_through_time_table_exports_curve(tmp_path: Path) -> None:
    output_path = tmp_path / "ltt.tsv"
    report = compute_lineage_through_time_curve(fixture("example_tree.nwk"))

    write_lineage_through_time_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "node\ttime_before_present\tlineage_count\tevent",
        "A|B|C|D\t0.3\t2\troot",
        "C|D\t0.2\t3\tspeciation",
        "A|B\t0.1\t4\tspeciation",
        "present\t0\t4\tpresent",
    ]


def test_inspect_diversification_time_tree_rejects_invalid_time_tree() -> None:
    with pytest.raises(DiversificationAnalysisError):
        inspect_diversification_time_tree(fixture("example_tree_no_lengths.nwk"))

    with pytest.raises(UnrootedTreeError):
        validate_time_tree_for_diversification(fixture("example_tree_unrooted.nwk"))


def test_detect_incomplete_taxon_sampling_metadata_reports_missing_and_invalid_rows() -> (
    None
):
    report = detect_incomplete_taxon_sampling_metadata(
        fixture("example_tree.nwk"),
        fixture("example_sampling_fractions_incomplete.tsv"),
    )

    assert report.complete is False
    assert report.sampling_column == "sampling_fraction"
    assert report.missing_taxa == ["D"]
    assert [issue.code for issue in report.invalid_rows] == [
        "missing-sampling-fraction",
        "out-of-range-sampling-fraction",
    ]


def test_estimate_diversification_rate_applies_sampling_correction() -> None:
    report = estimate_diversification_rate(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
        model="birth-death",
    )

    assert report.model == "birth-death"
    assert report.crown_age == 0.3
    assert report.observed_tip_count == 4
    assert report.sampling_fraction == 0.75
    assert report.corrected_tip_count == 5.33333333333333
    assert report.birth_rate >= report.net_diversification_rate
    assert report.aic > 0.0


def test_compare_diversification_models_returns_aic_ranked_rows() -> None:
    report = compare_diversification_models(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
    )

    assert report.better_model in {"yule", "birth-death"}
    assert [row.model for row in report.rows] == ["yule", "birth-death"]
    assert all(row.aic > 0.0 for row in report.rows)


def test_detect_diversification_outlier_clades_flags_high_and_low_clades() -> None:
    report = detect_diversification_outlier_clades(fixture("example_tree.nwk"))

    assert report.global_rate > 0.0
    assert {row.classification for row in report.observations} == {
        "baseline",
        "high",
        "low",
    }
    assert [row.node for row in report.high_diversification_clades] == ["A|B"]
    assert [row.node for row in report.low_diversification_clades] == ["C|D"]


def test_run_trait_dependent_diversification_analysis_reports_monophyly_and_polyphyly() -> (
    None
):
    monophyletic = run_trait_dependent_diversification_analysis(
        fixture("example_tree.nwk"),
        fixture("example_traits_diversification.tsv"),
        trait="habitat",
    )
    polyphyletic = run_trait_dependent_diversification_analysis(
        fixture("example_tree.nwk"),
        fixture("example_traits_diversification_polyphyletic.tsv"),
        trait="habitat",
    )

    assert [row.monophyletic for row in monophyletic.states] == [True, True]
    assert [row.diversification_rate for row in monophyletic.states] == [
        6.93147180559945,
        3.46573590279973,
    ]
    assert any("not monophyletic" in warning for warning in polyphyletic.warnings)


def test_write_diversification_tables_and_report_outputs_files(tmp_path: Path) -> None:
    clade_path = tmp_path / "clades.tsv"
    gamma_path = tmp_path / "gamma-statistic.tsv"
    trait_path = tmp_path / "trait-dependent.tsv"
    html_path = tmp_path / "diversification-report.html"

    clades = detect_diversification_outlier_clades(fixture("example_tree.nwk"))
    gamma = compute_diversification_gamma_statistic(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
    )
    trait_report = run_trait_dependent_diversification_analysis(
        fixture("example_tree.nwk"),
        fixture("example_traits_diversification.tsv"),
        trait="habitat",
    )
    report = render_diversification_report(
        tree_path=fixture("example_tree.nwk"),
        out_path=html_path,
        metadata_path=fixture("example_sampling_fractions.tsv"),
        traits_path=fixture("example_traits_diversification.tsv"),
        trait="habitat",
    )

    write_clade_diversification_table(clade_path, clades)
    write_diversification_gamma_statistic_table(gamma_path, gamma)
    write_trait_dependent_diversification_table(trait_path, trait_report)

    assert "classification" in clade_path.read_text(encoding="utf-8")
    assert "gamma_statistic" in gamma_path.read_text(encoding="utf-8")
    assert "monophyletic" in trait_path.read_text(encoding="utf-8")
    assert "diversification-model-comparison" in html_path.read_text(encoding="utf-8")
    assert "diversification-gamma-statistic" in html_path.read_text(encoding="utf-8")
    assert report.report_kind == "diversification"
