from __future__ import annotations

from pathlib import Path

from ..support.paths import REPOSITORY_ROOT

import pytest

from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    list_registered_evidence_studies,
    refresh_evidence_book,
)

REPO_ROOT = REPOSITORY_ROOT
pytestmark = pytest.mark.slow


def test_list_registered_evidence_studies_cover_all_pcm_studies() -> None:
    studies = list_registered_evidence_studies(REPO_ROOT)

    assert [study.study_id for study in studies] == [
        "cichlid-bayesian-linear-and-phylogenetic-mixed-models",
        "cichlid-phylogenetic-generalized-linear-mixed-models",
        "primate-continuous-and-categorical-regression",
        "primate-longevity-signal",
        "primate-pgls-and-signal",
    ]
    assert all(study.supports_partial_rerun for study in studies)


def test_refresh_evidence_book_writes_repo_index_outputs() -> None:
    report = refresh_evidence_book(REPO_ROOT)

    assert REPO_ROOT / DOCS_EVIDENCE_OVERVIEW in report.updated_paths
    assert (
        REPO_ROOT
        / "evidence-book"
        / "studies"
        / "primate-pgls-and-signal"
        / "evidence-002"
        / "results"
        / "reviewer-summary.json"
    ) in report.updated_paths
