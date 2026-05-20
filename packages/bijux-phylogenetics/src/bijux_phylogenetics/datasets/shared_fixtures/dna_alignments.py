from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedDnaAlignmentFixture:
    """One governed DNA alignment fixture shared across Bijux and external references."""

    fixture_id: str
    relative_path: str
    load_expectation: str
    ape_read_dna_expectation: str
    translation_expectation: str
    sequence_count: int
    alignment_length: int | None
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _fixtures_root() / self.relative_path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_dna_alignment_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_dna_alignment_fixtures() -> list[SharedDnaAlignmentFixture]:
    """Return the governed shared DNA alignment fixture corpus."""
    catalog = _load_catalog()
    return [
        SharedDnaAlignmentFixture(
            fixture_id=entry["fixture_id"],
            relative_path=entry["relative_path"],
            load_expectation=entry["load_expectation"],
            ape_read_dna_expectation=entry["ape_read_dna_expectation"],
            translation_expectation=entry["translation_expectation"],
            sequence_count=entry["sequence_count"],
            alignment_length=entry.get("alignment_length"),
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_dna_alignment_fixture(fixture_id: str) -> SharedDnaAlignmentFixture:
    """Resolve one governed shared DNA alignment fixture by durable fixture id."""
    for fixture in list_shared_dna_alignment_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_dna_alignment_fixtures())
    )
    raise ValueError(
        "unsupported shared DNA alignment fixture "
        f"'{fixture_id}'; expected one of: {supported}"
    )
