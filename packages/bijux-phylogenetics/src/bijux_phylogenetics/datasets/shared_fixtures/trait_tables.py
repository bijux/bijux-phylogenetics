from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedTraitTableFixture:
    """One governed trait-table fixture shared across Bijux and live reference checks."""

    fixture_id: str
    relative_path: str
    tree_fixture_id: str
    taxon_column: str
    primary_trait_columns: tuple[str, ...]
    standard_error_columns: tuple[str, ...]
    validation_expectation: str
    ape_reference_expectation: str
    row_count: int
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _fixtures_root() / self.relative_path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_trait_table_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_trait_table_fixtures() -> list[SharedTraitTableFixture]:
    """Return the governed shared trait-table fixture corpus."""
    catalog = _load_catalog()
    return [
        SharedTraitTableFixture(
            fixture_id=entry["fixture_id"],
            relative_path=entry["relative_path"],
            tree_fixture_id=entry["tree_fixture_id"],
            taxon_column=entry["taxon_column"],
            primary_trait_columns=tuple(entry["primary_trait_columns"]),
            standard_error_columns=tuple(entry.get("standard_error_columns", [])),
            validation_expectation=entry["validation_expectation"],
            ape_reference_expectation=entry["ape_reference_expectation"],
            row_count=entry["row_count"],
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_trait_table_fixture(fixture_id: str) -> SharedTraitTableFixture:
    """Resolve one governed shared trait-table fixture by durable fixture id."""
    for fixture in list_shared_trait_table_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_trait_table_fixtures())
    )
    raise ValueError(
        f"unsupported shared trait-table fixture '{fixture_id}'; expected one of: {supported}"
    )
