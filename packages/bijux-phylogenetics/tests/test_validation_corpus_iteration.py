from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.validation_corpus import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
    build_regression_dataset_corpus,
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
