from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedTreeSimulationFixture:
    """One governed random-tree simulation case shared across Bijux and `ape`."""

    fixture_id: str
    simulation_model: str
    reference_function: str
    replicate_count: int
    tip_count: int
    seed: int
    branch_length_model: str | None
    population_size: float | None
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def path(self) -> Path:
        return _catalog_path()


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_tree_simulation_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def list_shared_tree_simulation_fixtures() -> list[SharedTreeSimulationFixture]:
    """Return the governed random-tree simulation cases."""
    catalog = _load_catalog()
    return [
        SharedTreeSimulationFixture(
            fixture_id=entry["fixture_id"],
            simulation_model=entry["simulation_model"],
            reference_function=entry["reference_function"],
            replicate_count=entry["replicate_count"],
            tip_count=entry["tip_count"],
            seed=entry["seed"],
            branch_length_model=entry.get("branch_length_model"),
            population_size=entry.get("population_size"),
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_tree_simulation_fixture(
    fixture_id: str,
) -> SharedTreeSimulationFixture:
    """Resolve one governed simulation case by durable fixture id."""
    for fixture in list_shared_tree_simulation_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_tree_simulation_fixtures())
    )
    raise ValueError(
        "unsupported shared tree simulation fixture "
        f"{fixture_id!r}; expected one of: {supported}"
    )
