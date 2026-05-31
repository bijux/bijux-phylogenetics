"""Test-path bootstrap for repository-owned maintainer package tests."""

from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]

for source_root in (
    REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "src",
    REPO_ROOT / "packages" / "bijux-phylogenetics" / "src",
    REPO_ROOT / "packages" / "phylogenetic" / "src",
):
    source_root_text = str(source_root)
    if source_root_text not in sys.path:
        sys.path.insert(0, source_root_text)
