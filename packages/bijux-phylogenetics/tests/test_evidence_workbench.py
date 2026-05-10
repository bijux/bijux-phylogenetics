from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    build_evidence_book_selection,
    list_registered_evidence_studies,
    refresh_evidence_book,
    rerun_evidence_book_selection,
)


REPO_ROOT = Path(__file__).resolve().parents[3]


def test_list_registered_evidence_studies_reports_partial_rerun_support() -> None:
    studies = list_registered_evidence_studies(REPO_ROOT)

    assert [study.study_id for study in studies] == [
        "comparative-trust-boundaries",
        "primate-longevity-signal",
        "primate-pgls-and-signal",
        "taxon-trust",
    ]
    assert [study.supports_partial_rerun for study in studies] == [
        True,
        True,
        True,
        False,
    ]


def test_refresh_evidence_book_writes_docs_and_reviewer_outputs() -> None:
    report = refresh_evidence_book(REPO_ROOT)

    assert report.reviewer_summary_count == 23
    assert REPO_ROOT / DOCS_EVIDENCE_OVERVIEW in report.updated_paths
    assert (
        REPO_ROOT
        / "evidence-book"
        / "studies"
        / "primate-pgls-and-signal"
        / "evidence-002"
        / "reviewer-summary.json"
    ) in report.updated_paths


def test_rerun_evidence_book_selection_updates_selected_bundle_and_refreshes_book() -> (
    None
):
    report = rerun_evidence_book_selection(
        REPO_ROOT,
        "primate-pgls-and-signal",
        ["evidence-002"],
    )

    assert report.rerun_report.study_id == "primate-pgls-and-signal"
    assert report.rerun_report.selected_evidence_ids == ["evidence-002"]
    assert (
        "evidence-book/studies/primate-pgls-and-signal/evidence-002/manifest.json"
        in report.rerun_report.updated_paths
    )
    assert REPO_ROOT / DOCS_EVIDENCE_OVERVIEW in report.refresh_report.updated_paths


def test_build_evidence_book_selection_rebuilds_selected_evidence_only() -> None:
    report = build_evidence_book_selection(
        REPO_ROOT,
        "primate-pgls-and-signal",
        ["evidence-002"],
    )

    assert report.study_id == "primate-pgls-and-signal"
    assert report.selected_evidence_ids == ["evidence-002"]
    assert (
        "evidence-book/studies/primate-pgls-and-signal/evidence-002/manifest.json"
        in report.updated_paths
    )
    assert REPO_ROOT / DOCS_EVIDENCE_OVERVIEW in report.refresh_report.updated_paths
