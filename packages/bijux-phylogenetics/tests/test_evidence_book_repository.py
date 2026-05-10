from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.evidence.book import validate_evidence_book


REPO_ROOT = Path(__file__).resolve().parents[3]
STUDIES_ROOT = REPO_ROOT / "evidence-book" / "studies"


def test_repository_evidence_book_passes_validation() -> None:
    report = validate_evidence_book(REPO_ROOT)

    assert report.valid is True, [
        f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
    ]


def test_repository_evidence_book_only_keeps_pcm_studies() -> None:
    assert sorted(path.name for path in STUDIES_ROOT.iterdir() if path.is_dir()) == [
        "primate-longevity-signal",
        "primate-pgls-and-signal",
    ]


def test_repository_evidence_book_study_roots_are_minimal() -> None:
    assert {path.name for path in (STUDIES_ROOT / "primate-longevity-signal").iterdir()} == {
        "README.md",
        "datasets",
        "evidence-001",
        "evidence-002",
        "evidence-003",
        "evidence-004",
        "evidence-005",
        "evidence-006",
        "evidence-007",
        "evidence-008",
        "evidence-009",
        "provenance",
        "reference",
    }
    assert {path.name for path in (STUDIES_ROOT / "primate-pgls-and-signal").iterdir()} == {
        "README.md",
        "datasets",
        "evidence-001",
        "evidence-002",
        "evidence-003",
        "evidence-004",
        "evidence-005",
        "evidence-006",
        "evidence-007",
        "evidence-008",
        "evidence-009",
        "evidence-010",
        "provenance",
        "reference",
    }
