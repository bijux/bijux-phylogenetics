from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedDistanceMatrixFixture:
    """One governed distance-matrix fixture shared across Bijux and `ape`."""

    fixture_id: str
    relative_path: str
    taxa: tuple[str, ...]
    pair_count: int
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _fixtures_root() / self.relative_path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_distance_matrix_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_distance_matrix_fixtures() -> list[SharedDistanceMatrixFixture]:
    """Return the governed shared distance-matrix fixture corpus."""
    catalog = _load_catalog()
    return [
        SharedDistanceMatrixFixture(
            fixture_id=entry["fixture_id"],
            relative_path=entry["relative_path"],
            taxa=tuple(entry["taxa"]),
            pair_count=entry["pair_count"],
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_distance_matrix_fixture(fixture_id: str) -> SharedDistanceMatrixFixture:
    """Resolve one governed shared distance-matrix fixture by durable fixture id."""
    for fixture in list_shared_distance_matrix_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_distance_matrix_fixtures())
    )
    raise ValueError(
        "unsupported shared distance-matrix fixture "
        f"{fixture_id!r}; expected one of: {supported}"
    )
