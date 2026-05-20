from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedTreeSetFixture:
    """One governed tree-set fixture shared across Bijux and external references."""

    fixture_id: str
    relative_path: str
    source_format: str
    parse_expectation: str
    tree_count: int
    shared_taxa: tuple[str, ...]
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _fixtures_root() / self.relative_path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_tree_set_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_tree_set_fixtures() -> list[SharedTreeSetFixture]:
    """Return the governed shared tree-set fixture corpus."""
    catalog = _load_catalog()
    return [
        SharedTreeSetFixture(
            fixture_id=entry["fixture_id"],
            relative_path=entry["relative_path"],
            source_format=entry["source_format"],
            parse_expectation=entry["parse_expectation"],
            tree_count=entry["tree_count"],
            shared_taxa=tuple(entry["shared_taxa"]),
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_tree_set_fixture(fixture_id: str) -> SharedTreeSetFixture:
    """Resolve one governed shared tree-set fixture by durable fixture id."""
    for fixture in list_shared_tree_set_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_tree_set_fixtures())
    )
    raise ValueError(
        f"unsupported shared tree-set fixture '{fixture_id}'; expected one of: {supported}"
    )
