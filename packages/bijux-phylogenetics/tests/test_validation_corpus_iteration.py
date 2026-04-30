from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.validation_corpus import build_clean_benchmark_corpus


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
