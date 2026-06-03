from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_distance_matrix_fixture,
    get_shared_dna_alignment_fixture,
    get_shared_trait_table_fixture,
    get_shared_tree_fixture,
    get_shared_tree_set_fixture,
    get_shared_tree_simulation_fixture,
)


def package_root() -> Path:
    return Path(__file__).resolve().parents[5]


def default_fixtures_root() -> Path:
    return package_root() / "tests" / "fixtures"


@dataclass(frozen=True, slots=True)
class ApeParityFixtureResolver:
    """Resolve governed fixture paths for the live `ape` parity suite."""

    fixtures_root: Path | None = None

    def fixture_path(self, fixture_kind: str, fixture_id: str) -> Path:
        if fixture_kind == "tree":
            fixture = get_shared_tree_fixture(fixture_id)
        elif fixture_kind == "tree-set":
            fixture = get_shared_tree_set_fixture(fixture_id)
        elif fixture_kind == "dna-alignment":
            fixture = get_shared_dna_alignment_fixture(fixture_id)
        elif fixture_kind == "distance-matrix":
            fixture = get_shared_distance_matrix_fixture(fixture_id)
        elif fixture_kind == "simulation":
            fixture = get_shared_tree_simulation_fixture(fixture_id)
        else:
            raise ValueError(f"unsupported ape parity fixture kind '{fixture_kind}'")
        if self.fixtures_root is None:
            return fixture.path
        if fixture_kind == "simulation":
            return (
                self.fixtures_root
                / "metadata"
                / "shared_tree_simulation_fixture_catalog.json"
            )
        return self.fixtures_root / fixture.relative_path

    def trait_path(self, fixture_id: str) -> Path:
        fixture = get_shared_trait_table_fixture(fixture_id)
        if self.fixtures_root is None:
            return fixture.path
        return self.fixtures_root / fixture.relative_path
