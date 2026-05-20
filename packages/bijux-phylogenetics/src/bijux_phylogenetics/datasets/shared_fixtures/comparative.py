from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root
from .trait_tables import (
    SharedTraitTableFixture,
    get_shared_trait_table_fixture,
)
from .trees import (
    SharedTreeFixture,
    get_shared_tree_fixture,
)


@dataclass(frozen=True, slots=True)
class SharedPhytoolsComparativeFixture:
    """One governed comparative fixture pairing for live `phytools` parity."""

    fixture_id: str
    tree_fixture_id: str
    trait_table_fixture_id: str
    trait_name: str
    taxon_column: str
    trait_kind: str
    validation_expectation: str
    phytools_reference_expectation: str
    feature_tags: tuple[str, ...]
    simulation_metadata: dict[str, object] | None
    notes: str

    @property
    def tree_fixture(self) -> SharedTreeFixture:
        return get_shared_tree_fixture(self.tree_fixture_id)

    @property
    def trait_table_fixture(self) -> SharedTraitTableFixture:
        return get_shared_trait_table_fixture(self.trait_table_fixture_id)

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
        _fixtures_root()
        / "metadata"
        / "shared_phytools_comparative_fixture_catalog.json"
    )


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_phytools_comparative_fixtures() -> list[
    SharedPhytoolsComparativeFixture
]:
    """Return the governed comparative fixture corpus used by live `phytools` parity."""
    catalog = _load_catalog()
    return [
        SharedPhytoolsComparativeFixture(
            fixture_id=entry["fixture_id"],
            tree_fixture_id=entry["tree_fixture_id"],
            trait_table_fixture_id=entry["trait_table_fixture_id"],
            trait_name=entry["trait_name"],
            taxon_column=entry["taxon_column"],
            trait_kind=entry["trait_kind"],
            validation_expectation=entry["validation_expectation"],
            phytools_reference_expectation=entry["phytools_reference_expectation"],
            feature_tags=tuple(entry["feature_tags"]),
            simulation_metadata=entry.get("simulation_metadata"),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_phytools_comparative_fixture(
    fixture_id: str,
) -> SharedPhytoolsComparativeFixture:
    """Resolve one governed `phytools` comparative fixture by durable fixture id."""
    for fixture in list_shared_phytools_comparative_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(
            fixture.fixture_id
            for fixture in list_shared_phytools_comparative_fixtures()
        )
    )
    raise ValueError(
        f"unsupported shared phytools comparative fixture '{fixture_id}'; expected one of: {supported}"
    )
