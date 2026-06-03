from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root
from .trait_tables import SharedTraitTableFixture, get_shared_trait_table_fixture
from .trees import SharedTreeFixture, get_shared_tree_fixture


@dataclass(frozen=True, slots=True)
class SharedGeigerContinuousFixture:
    """One governed continuous-trait fixture for live `geiger::fitContinuous` parity."""

    fixture_id: str
    trait_table_fixture_id: str
    trait_name: str
    standard_error_trait_name: str | None
    taxon_column: str
    supported_model_names: tuple[str, ...]
    validation_expectation: str
    geiger_reference_expectation: str
    feature_tags: tuple[str, ...]
    simulation_metadata: dict[str, object] | None
    notes: str

    @property
    def trait_table_fixture(self) -> SharedTraitTableFixture:
        return get_shared_trait_table_fixture(self.trait_table_fixture_id)

    @property
    def tree_fixture(self) -> SharedTreeFixture:
        return get_shared_tree_fixture(self.trait_table_fixture.tree_fixture_id)

    @property
    def tree_path(self) -> Path:
        return self.tree_fixture.path

    @property
    def traits_path(self) -> Path:
        return self.trait_table_fixture.path


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return (
        _fixtures_root() / "metadata" / "shared_geiger_continuous_fixture_catalog.json"
    )


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_geiger_continuous_fixtures() -> list[SharedGeigerContinuousFixture]:
    """Return the governed continuous-trait corpus for live `geiger` parity."""
    catalog = _load_catalog()
    return [
        SharedGeigerContinuousFixture(
            fixture_id=entry["fixture_id"],
            trait_table_fixture_id=entry["trait_table_fixture_id"],
            trait_name=entry["trait_name"],
            standard_error_trait_name=entry.get("standard_error_trait_name"),
            taxon_column=entry["taxon_column"],
            supported_model_names=tuple(entry["supported_model_names"]),
            validation_expectation=entry["validation_expectation"],
            geiger_reference_expectation=entry["geiger_reference_expectation"],
            feature_tags=tuple(entry["feature_tags"]),
            simulation_metadata=entry.get("simulation_metadata"),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_geiger_continuous_fixture(
    fixture_id: str,
) -> SharedGeigerContinuousFixture:
    """Resolve one governed `geiger` continuous fixture by durable id."""
    for fixture in list_shared_geiger_continuous_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(
            fixture.fixture_id for fixture in list_shared_geiger_continuous_fixtures()
        )
    )
    raise ValueError(
        f"unsupported shared geiger continuous fixture '{fixture_id}'; expected one of: {supported}"
    )
