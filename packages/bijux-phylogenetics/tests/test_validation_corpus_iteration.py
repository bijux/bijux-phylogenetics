from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.benchmark import (
    benchmark_alignment_site_scaling,
    benchmark_tree_set_consensus,
)
from bijux_phylogenetics.validation_corpus import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_memory_benchmark_dashboard,
    build_method_accuracy_dashboard,
    build_method_limitation_registry,
    build_messy_benchmark_corpus,
    build_regression_dataset_corpus,
    build_runtime_benchmark_dashboard,
    build_scientific_validation_report,
    validate_simulation_reproducibility,
)


FIXTURES = Path(__file__).parent / "fixtures"


def test_build_clean_benchmark_corpus_keeps_core_dataset_ready() -> None:
    report = build_clean_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 242
    assert report.passed is True
    assert report.passed_case_count == 1
    case = report.cases[0]
    assert case.readiness_decision == "ready_with_warnings"
    assert "comparative" in case.allowed_analyses
    assert "publication" in case.allowed_analyses
    assert case.blockers == []


def test_build_broken_benchmark_corpus_preserves_expected_failure_signatures() -> None:
    report = build_broken_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 243
    assert report.passed is True
    observed = {case.name: case for case in report.cases}
    assert observed["duplicate_tip_tree"].observed_code == "duplicate_taxon_error"
    assert observed["invalid_alignment_lengths"].observed_code == "invalid_alignment_error"
    assert observed["dataset_missing_metadata_taxon"].readiness_decision == "blocked"
    assert "metadata table is missing one or more tree taxa" in observed["dataset_missing_metadata_taxon"].blockers


def test_build_messy_benchmark_corpus_captures_multi_surface_warning_cases() -> None:
    report = build_messy_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 244
    assert report.passed is True
    observed = {case.name: case for case in report.cases}
    first = observed["reordered_alignment_extra_taxa_invalid_dates_and_calibrations"]
    assert "calibration table contains invalid fossil calibration targets or ages" in first.blockers
    assert "one or more dataset surfaces silently reorder shared taxa relative to the tree" in first.warnings
    second = observed["low_information_alignment_with_trait_mismatch"]
    assert "alignment is not currently safe for core inference workflows" in second.blockers
    assert "alignment contains near-duplicate sequences" in second.warnings


def test_build_regression_dataset_corpus_matches_checked_in_summaries() -> None:
    report = build_regression_dataset_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 245
    assert report.passed is True
    assert report.passed_case_count == 2
    observed = {case.name: case for case in report.cases}
    assert observed["core_inference_ready_dataset"].observed["risk_level"] == "low"
    assert observed["warning_rich_dataset"].observed["warning_count"] == 12


def test_benchmark_alignment_site_scaling_reports_site_axis_observations() -> None:
    report = benchmark_alignment_site_scaling(replicates=1, site_counts=[24, 48], sequence_count=4)

    assert report.sequence_count == 4
    assert [row.item_count for row in report.observations] == [24, 48]


def test_benchmark_tree_set_consensus_reports_tree_count_scaling() -> None:
    report = benchmark_tree_set_consensus(replicates=1, tree_counts=[4, 8], tip_count=4)

    assert report.tip_count == 4
    assert [row.item_count for row in report.observations] == [4, 8]


def test_build_method_accuracy_dashboard_summarizes_fixture_and_corpus_pass_rates() -> None:
    report = build_method_accuracy_dashboard(fixtures_root=FIXTURES)

    assert report.goal_id == 246
    surfaces = {row.surface: row for row in report.rows}
    assert surfaces["level1-reference-validation"].coverage_count > 0
    assert surfaces["clean-benchmark-corpus"].accuracy == 1.0
    assert surfaces["regression-dataset-corpus"].failed_count == 0


def test_build_runtime_and_memory_dashboards_cover_sites_and_posterior_samples() -> None:
    runtime = build_runtime_benchmark_dashboard(replicates=1)
    memory = build_memory_benchmark_dashboard(replicates=1)

    assert runtime.goal_id == 247
    assert memory.goal_id == 248
    assert {row.scaling_axis for row in runtime.rows} >= {"sites", "posterior_samples", "taxa"}
    assert {row.scaling_axis for row in memory.rows} >= {"sites", "posterior_samples", "taxa"}


def test_build_method_limitation_registry_marks_experimental_and_validated_surfaces() -> None:
    report = build_method_limitation_registry()

    assert report.goal_id == 250
    statuses = {entry.method: entry.status for entry in report.entries}
    assert statuses["tree-validation"] == "validated"
    assert statuses["bayesian-time-tree"] == "experimental"


def test_build_scientific_validation_report_separates_claim_statuses() -> None:
    report = build_scientific_validation_report(fixtures_root=FIXTURES)

    assert report.goal_id == 249
    statuses = {claim.status for claim in report.claims}
    assert statuses >= {"validated", "experimental", "unvalidated", "unsafe"}


def test_validate_simulation_reproducibility_confirms_same_seed_repeatability() -> None:
    report = validate_simulation_reproducibility(fixtures_root=FIXTURES)

    assert report.goal_id == 251
    assert report.passed is True
    assert len(report.cases) == 6
    assert all(case.digest for case in report.cases)
