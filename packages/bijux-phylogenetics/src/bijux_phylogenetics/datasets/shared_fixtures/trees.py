from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedTreeFixture:
    """One governed tree fixture shared across Bijux and external references."""

    fixture_id: str
    relative_path: str
    parse_expectation: str
    validation_expectation: str
    ape_read_tree_expectation: str
    tip_count: int | None
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _fixtures_root() / self.relative_path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_tree_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_tree_fixtures() -> list[SharedTreeFixture]:
    """Return the governed shared tree fixture corpus."""
    catalog = _load_catalog()
    return [
        SharedTreeFixture(
            fixture_id=entry["fixture_id"],
            relative_path=entry["relative_path"],
            parse_expectation=entry["parse_expectation"],
            validation_expectation=entry["validation_expectation"],
            ape_read_tree_expectation=entry["ape_read_tree_expectation"],
            tip_count=entry.get("tip_count"),
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_tree_fixture(fixture_id: str) -> SharedTreeFixture:
    """Resolve one governed shared tree fixture by durable fixture id."""
    for fixture in list_shared_tree_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_tree_fixtures())
    )
    raise ValueError(
        f"unsupported shared tree fixture '{fixture_id}'; expected one of: {supported}"
    )
