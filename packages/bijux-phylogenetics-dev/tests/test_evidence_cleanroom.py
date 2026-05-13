from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics_dev.quality.evidence_cleanroom import (
    build_evidence_cleanroom_report,
    build_selected_evidence_cleanroom_reports,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.slow
def test_repository_cleanroom_rerun_keeps_primate_longevity_selection_clean() -> None:
    report = build_evidence_cleanroom_report(
        REPO_ROOT,
        study_id="primate-longevity-signal",
        evidence_ids=["evidence-002"],
    )

    assert report.validation_issue_count == 0
    assert report.artifact_issue_count == 0
    assert report.worktree_clean is True
    assert report.updated_path_count > 0


def test_cleanroom_rerun_requires_an_explicit_evidence_selection() -> None:
    with pytest.raises(ValueError, match="at least one evidence id is required"):
        build_evidence_cleanroom_report(
            REPO_ROOT,
            study_id="primate-longevity-signal",
            evidence_ids=[],
        )


@pytest.mark.slow
def test_repository_selected_cleanroom_reruns_keep_governed_selections_clean() -> None:
    report = build_selected_evidence_cleanroom_reports(REPO_ROOT)

    assert report.selection_count >= 1
    assert all(entry.validation_issue_count == 0 for entry in report.reports)
    assert all(entry.artifact_issue_count == 0 for entry in report.reports)
    assert all(entry.worktree_clean for entry in report.reports)
