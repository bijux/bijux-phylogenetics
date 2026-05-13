from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def test_repository_exposes_contributor_onboarding_files() -> None:
    code_of_conduct = (REPO_ROOT / "CODE_OF_CONDUCT.md").read_text(encoding="utf-8")
    contributing = (REPO_ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")

    assert "## Reporting" in code_of_conduct
    assert "mailto:bijan@bijux.io" in code_of_conduct
    assert "## Validation Expectations" in contributing
    assert "`make help`" in contributing
    assert "`artifacts/`" in contributing
