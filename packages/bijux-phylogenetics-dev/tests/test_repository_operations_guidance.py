from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def test_root_readme_documents_local_artifact_contract() -> None:
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")

    assert "transient local outputs belong under `artifacts/`" in readme
    assert "`artifacts/root/check-venv/`" in readme
    assert "`artifacts/root/docs/site/`" in readme
