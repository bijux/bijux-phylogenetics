from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.validation import (
    validate_continuous_ancestral_reconstruction,
    validate_discrete_ancestral_reconstruction,
)


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_validate_continuous_ancestral_reconstruction_reports_error_against_truth() -> None:
    report = validate_continuous_ancestral_reconstruction(
        fixture("example_tree.nwk"),
        simulation_model="brownian",
        reconstruction_model="brownian",
        replicates=3,
        seed=11,
    )
    assert report.replicates == 3
    assert report.node_count == 9
    assert report.mean_absolute_error >= 0.0
    assert report.root_mean_squared_error >= 0.0


def test_validate_discrete_ancestral_reconstruction_reports_accuracy_and_calibration() -> None:
    report = validate_discrete_ancestral_reconstruction(
        fixture("example_tree.nwk"),
        reconstruction_model="fitch",
        replicates=3,
        seed=11,
        states=["forest", "tundra", "desert"],
    )
    assert report.replicates == 3
    assert report.node_count == 9
    assert 0.0 <= report.accuracy <= 1.0
    assert 0.0 <= report.mean_true_state_probability <= 1.0
    assert 0.0 <= report.uncertainty_calibration_gap <= 1.0
