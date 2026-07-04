from __future__ import annotations

import os
from pathlib import Path

_PACKAGE_MARKER = Path("src") / "bijux_phylogenetics"
_REPOSITORY_PACKAGE_MARKER = Path("packages") / "bijux-phylogenetics"
_FIXTURE_MARKER = Path("tests") / "fixtures" / "metadata"


def package_root(*, start: Path | None = None) -> Path:
    """Locate the repository package root for local governed assets."""
    explicit_root = _explicit_package_root()
    if explicit_root is not None:
        return explicit_root
    for candidate in _search_candidates(start=start):
        if _is_package_root(candidate):
            return candidate
        nested_package_root = candidate / _REPOSITORY_PACKAGE_MARKER
        if _is_package_root(nested_package_root):
            return nested_package_root
    raise FileNotFoundError(
        "could not locate the bijux-phylogenetics package root from the active environment"
    )


def repository_root(*, start: Path | None = None) -> Path:
    """Locate the repository root that owns the phylogenetics package checkout."""
    explicit_root = _explicit_repository_root()
    if explicit_root is not None:
        return explicit_root
    discovered_package_root = package_root(start=start)
    if discovered_package_root.parent.name == "packages":
        return discovered_package_root.parents[1]
    return discovered_package_root.parent


def governed_fixture_root(*, start: Path | None = None) -> Path:
    """Locate the repository-owned governed fixture tree."""
    explicit_root = _explicit_fixture_root()
    if explicit_root is not None:
        return explicit_root
    root = package_root(start=start) / "tests" / "fixtures"
    if _is_fixture_root(root):
        return root
    raise FileNotFoundError(
        "could not locate the bijux-phylogenetics governed fixture tree"
    )


def _search_candidates(*, start: Path | None) -> list[Path]:
    candidates: list[Path] = []
    seen: set[Path] = set()
    start_points = [Path.cwd()]
    if start is not None:
        start_points.append(start)
    else:
        start_points.append(Path(__file__))
    for start_point in start_points:
        resolved = start_point.resolve()
        lineage = [resolved]
        lineage.extend(resolved.parents)
        for candidate in lineage:
            if candidate in seen:
                continue
            seen.add(candidate)
            candidates.append(candidate)
    return candidates


def _explicit_repository_root() -> Path | None:
    value = os.environ.get("BIJUX_PHYLOGENETICS_REPOSITORY_ROOT")
    if not value:
        return None
    candidate = Path(value).expanduser().resolve()
    if _is_repository_root(candidate):
        return candidate
    package_candidate = candidate / _REPOSITORY_PACKAGE_MARKER
    if _is_package_root(package_candidate):
        return candidate
    return None


def _explicit_package_root() -> Path | None:
    value = os.environ.get("BIJUX_PHYLOGENETICS_PACKAGE_ROOT")
    if not value:
        return None
    candidate = Path(value).expanduser().resolve()
    if _is_package_root(candidate):
        return candidate
    package_candidate = candidate / _REPOSITORY_PACKAGE_MARKER
    if _is_package_root(package_candidate):
        return package_candidate
    return None


def _explicit_fixture_root() -> Path | None:
    value = os.environ.get("BIJUX_PHYLOGENETICS_FIXTURE_ROOT")
    if not value:
        return None
    candidate = Path(value).expanduser().resolve()
    if _is_fixture_root(candidate):
        return candidate
    return None


def _is_package_root(candidate: Path) -> bool:
    return candidate.is_dir() and (candidate / _PACKAGE_MARKER).is_dir()


def _is_repository_root(candidate: Path) -> bool:
    return candidate.is_dir() and _is_package_root(candidate / _REPOSITORY_PACKAGE_MARKER)


def _is_fixture_root(candidate: Path) -> bool:
    return candidate.is_dir() and (candidate / "metadata").is_dir()
