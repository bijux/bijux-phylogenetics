from __future__ import annotations

import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_ou_traits,
    simulate_protein_alignment,
)

from .contracts import SimulationReproducibilityCase, SimulationReproducibilityReport
from .dataset_corpora import default_fixtures_root, fixture
from .presentation import normalize_jsonable


def validate_simulation_reproducibility(
    *, fixtures_root: Path | None = None
) -> SimulationReproducibilityReport:
    """Verify that repeated simulations with the same seed produce identical structured results."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = fixture(root, "trees", "example_tree.nwk")

    def digest(payload: object) -> tuple[str, str]:
        normalized = normalize_jsonable(payload)
        encoded = json.dumps(normalized, sort_keys=True).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest(), encoded.decode("utf-8")

    cases: list[tuple[str, object]] = [
        (
            "birth-death-tree-set",
            lambda: (
                simulate_birth_death_trees(tree_count=3, tip_count=4, seed=7)[1].records
            ),
        ),
        (
            "brownian-traits",
            lambda: simulate_brownian_traits(
                tree_path, root_state=1.0, sigma=0.5, seed=7
            ),
        ),
        (
            "ou-traits",
            lambda: simulate_ou_traits(
                tree_path, root_state=1.0, sigma=0.5, alpha=0.7, theta=0.2, seed=7
            ),
        ),
        (
            "discrete-traits",
            lambda: simulate_discrete_traits(
                tree_path, states=["north", "south"], transition_rate=0.8, seed=7
            ),
        ),
        (
            "dna-alignment",
            lambda: simulate_dna_alignment(
                tree_path, sequence_length=16, substitution_rate=0.9, seed=7
            ),
        ),
        (
            "protein-alignment",
            lambda: simulate_protein_alignment(
                tree_path, sequence_length=12, substitution_rate=0.6, seed=7
            ),
        ),
    ]

    results: list[SimulationReproducibilityCase] = []
    for surface, callback in cases:
        first = callback()
        second = callback()
        first_digest, first_payload = digest(first)
        second_digest, second_payload = digest(second)
        notes: list[str] = []
        if first_digest != second_digest or first_payload != second_payload:
            notes.append("same-seed simulation output drifted between repeated runs")
        results.append(
            SimulationReproducibilityCase(
                surface=surface,
                passed=not notes,
                digest=first_digest,
                notes=notes,
            )
        )
    return SimulationReproducibilityReport(
        goal_id=251,
        passed=all(case.passed for case in results),
        cases=results,
        limitations=[
            "the reproducibility check covers deterministic seeded library surfaces and does not yet assert cross-environment bit-for-bit stability",
        ],
    )
