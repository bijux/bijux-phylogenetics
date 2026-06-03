from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..continuous import (
    simulate_brownian_trait_collection,
    simulate_speciational_trait_collection,
)
from ..discrete import simulate_discrete_histories


@dataclass(frozen=True, slots=True)
class _GeigerDiscreteRateRow:
    source_state: str
    target_state: str
    rate: float


_CONTINUOUS_MEAN_TOLERANCE = 0.12
_CONTINUOUS_STANDARD_DEVIATION_TOLERANCE = 0.08
_CONTINUOUS_COVARIANCE_TOLERANCE = 0.08
_CONTINUOUS_CORRELATION_TOLERANCE = 0.12
_DISCRETE_TIP_STATE_FREQUENCY_TOLERANCE = 0.05

_REFERENCE_TREE_FIXTURE = "example_tree_internal_long_branch.nwk"
_CONTINUOUS_REFERENCE_ROOT_STATE = 1.5
_CONTINUOUS_REFERENCE_SIGMA_SQUARED = 0.25
_CONTINUOUS_REFERENCE_REPLICATES = 256
_DISCRETE_REFERENCE_REPLICATES = 1024
_REFERENCE_SEED = 11
_DISCRETE_REFERENCE_RATE_ROWS = (
    _GeigerDiscreteRateRow(source_state="0", target_state="1", rate=0.4),
    _GeigerDiscreteRateRow(source_state="1", target_state="0", rate=0.1),
)

_SELECTED_COVARIANCE_LABELS = ("A|B", "A|C", "C|D")


@dataclass(slots=True)
class GeigerSimCharReferenceObservation:
    """One governed simulation-envelope comparison against local geiger::sim.char."""

    case_id: str
    model: str
    reference_source: str
    input_fixtures: list[Path]
    true_parameters: dict[str, object]
    expected_metrics: dict[str, object]
    observed_metrics: dict[str, object]
    passed: bool
    notes: list[str]


@dataclass(slots=True)
class GeigerSimCharReferenceValidationReport:
    """Integrated governed envelope validation across local geiger::sim.char cases."""

    case_count: int
    all_passed: bool
    observations: list[GeigerSimCharReferenceObservation]


GEIGER_SIM_CHAR_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "brownian-internal-long-branch-envelope": {
        "expected_metrics": {
            "tip_distribution": {
                "A": {"mean_value": 1.459827, "standard_deviation": 0.5256817},
                "B": {"mean_value": 1.460553, "standard_deviation": 0.5319871},
                "C": {"mean_value": 1.510404, "standard_deviation": 0.2266562},
                "D": {"mean_value": 1.511171, "standard_deviation": 0.2242928},
            },
            "tip_covariance": {
                "A|B": {"covariance": 0.2560237, "correlation": 0.9154954},
                "A|C": {"covariance": -0.004060494, "correlation": -0.03407912},
                "C|D": {"covariance": 0.02569813, "correlation": 0.505497},
            },
        }
    },
    "speciational-internal-long-branch-envelope": {
        "expected_metrics": {
            "tip_distribution": {
                "A": {"mean_value": 1.467032, "standard_deviation": 0.686249},
                "B": {"mean_value": 1.469327, "standard_deviation": 0.7034484},
                "C": {"mean_value": 1.532902, "standard_deviation": 0.7167499},
                "D": {"mean_value": 1.535325, "standard_deviation": 0.7092762},
            },
            "tip_covariance": {
                "A|B": {"covariance": 0.2463678, "correlation": 0.5103523},
                "A|C": {"covariance": -0.01455987, "correlation": -0.02960112},
                "C|D": {"covariance": 0.2569813, "correlation": 0.505497},
            },
        }
    },
    "discrete-rate-matrix-internal-long-branch-envelope": {
        "expected_metrics": {
            "tip_state_frequency": {
                "A:0": 0.6513672,
                "A:1": 0.3486328,
                "B:0": 0.65625,
                "B:1": 0.34375,
                "C:0": 0.9189453,
                "C:1": 0.08105469,
                "D:0": 0.9277344,
                "D:1": 0.07226562,
            }
        }
    },
}


def validate_geiger_sim_char_reference_examples() -> (
    GeigerSimCharReferenceValidationReport
):
    """Validate governed simulation envelopes against local geiger::sim.char behavior."""
    observations = [
        _validate_brownian_reference_case(),
        _validate_speciational_reference_case(),
        _validate_discrete_rate_matrix_reference_case(),
    ]
    return GeigerSimCharReferenceValidationReport(
        case_count=len(observations),
        all_passed=all(observation.passed for observation in observations),
        observations=observations,
    )


def _validate_brownian_reference_case() -> GeigerSimCharReferenceObservation:
    tree_path = _tree_fixture_path(_REFERENCE_TREE_FIXTURE)
    report = simulate_brownian_trait_collection(
        tree_path,
        root_state=_CONTINUOUS_REFERENCE_ROOT_STATE,
        sigma_squared=_CONTINUOUS_REFERENCE_SIGMA_SQUARED,
        replicates=_CONTINUOUS_REFERENCE_REPLICATES,
        seed=_REFERENCE_SEED,
    )
    payload = GEIGER_SIM_CHAR_REFERENCE_PAYLOADS[
        "brownian-internal-long-branch-envelope"
    ]
    observed_metrics = _continuous_metrics(report)
    passed = _continuous_metrics_match(
        expected=payload["expected_metrics"],
        observed=observed_metrics,
    )
    return GeigerSimCharReferenceObservation(
        case_id="brownian-internal-long-branch-envelope",
        model="BM",
        reference_source="stored local geiger::sim.char Brownian summary envelope",
        input_fixtures=[tree_path],
        true_parameters=_continuous_true_parameters(),
        expected_metrics=payload["expected_metrics"],
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "Brownian simulation summary drifted outside the governed geiger envelope"
        ],
    )


def _validate_speciational_reference_case() -> GeigerSimCharReferenceObservation:
    tree_path = _tree_fixture_path(_REFERENCE_TREE_FIXTURE)
    report = simulate_speciational_trait_collection(
        tree_path,
        root_state=_CONTINUOUS_REFERENCE_ROOT_STATE,
        sigma_squared=_CONTINUOUS_REFERENCE_SIGMA_SQUARED,
        replicates=_CONTINUOUS_REFERENCE_REPLICATES,
        seed=_REFERENCE_SEED,
    )
    payload = GEIGER_SIM_CHAR_REFERENCE_PAYLOADS[
        "speciational-internal-long-branch-envelope"
    ]
    observed_metrics = _continuous_metrics(report)
    passed = _continuous_metrics_match(
        expected=payload["expected_metrics"],
        observed=observed_metrics,
    )
    return GeigerSimCharReferenceObservation(
        case_id="speciational-internal-long-branch-envelope",
        model="speciational",
        reference_source="stored local geiger::sim.char speciational summary envelope",
        input_fixtures=[tree_path],
        true_parameters=_continuous_true_parameters(),
        expected_metrics=payload["expected_metrics"],
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "Speciational simulation summary drifted outside the governed geiger envelope"
        ],
    )


def _validate_discrete_rate_matrix_reference_case() -> (
    GeigerSimCharReferenceObservation
):
    tree_path = _tree_fixture_path(_REFERENCE_TREE_FIXTURE)
    report = simulate_discrete_histories(
        tree_path,
        states=["0", "1"],
        rate_rows=list(_DISCRETE_REFERENCE_RATE_ROWS),
        root_state="0",
        replicates=_DISCRETE_REFERENCE_REPLICATES,
        seed=_REFERENCE_SEED,
    )
    payload = GEIGER_SIM_CHAR_REFERENCE_PAYLOADS[
        "discrete-rate-matrix-internal-long-branch-envelope"
    ]
    observed_metrics = _discrete_tip_state_frequency_metrics(report)
    passed = _discrete_tip_state_frequency_metrics_match(
        expected=payload["expected_metrics"]["tip_state_frequency"],
        observed=observed_metrics["tip_state_frequency"],
    )
    return GeigerSimCharReferenceObservation(
        case_id="discrete-rate-matrix-internal-long-branch-envelope",
        model="discrete",
        reference_source="stored local geiger::sim.char discrete tip-state envelope",
        input_fixtures=[tree_path],
        true_parameters={
            "states": ["0", "1"],
            "root_state": "0",
            "rate_rows": [
                {
                    "source_state": row.source_state,
                    "target_state": row.target_state,
                    "rate": row.rate,
                }
                for row in _DISCRETE_REFERENCE_RATE_ROWS
            ],
            "replicate_count": _DISCRETE_REFERENCE_REPLICATES,
            "seed": _REFERENCE_SEED,
        },
        expected_metrics=payload["expected_metrics"],
        observed_metrics=observed_metrics,
        passed=passed,
        notes=[]
        if passed
        else [
            "Discrete tip-state frequencies drifted outside the governed geiger envelope"
        ],
    )


def _continuous_true_parameters() -> dict[str, object]:
    return {
        "root_state": _CONTINUOUS_REFERENCE_ROOT_STATE,
        "sigma_squared": _CONTINUOUS_REFERENCE_SIGMA_SQUARED,
        "replicate_count": _CONTINUOUS_REFERENCE_REPLICATES,
        "seed": _REFERENCE_SEED,
    }


def _continuous_metrics(report) -> dict[str, object]:
    return {
        "tip_distribution": {
            row.label: {
                "mean_value": row.mean_value,
                "standard_deviation": row.standard_deviation,
            }
            for row in report.rows
            if row.row_kind == "tip_distribution"
        },
        "tip_covariance": {
            row.label: {
                "covariance": row.covariance,
                "correlation": row.correlation,
            }
            for row in report.rows
            if row.row_kind == "tip_covariance"
            and row.label in _SELECTED_COVARIANCE_LABELS
        },
    }


def _continuous_metrics_match(
    *,
    expected: dict[str, object],
    observed: dict[str, object],
) -> bool:
    expected_tip_distribution = expected["tip_distribution"]
    observed_tip_distribution = observed["tip_distribution"]
    for taxon, expected_row in expected_tip_distribution.items():
        observed_row = observed_tip_distribution[taxon]
        if (
            abs(observed_row["mean_value"] - expected_row["mean_value"])
            > _CONTINUOUS_MEAN_TOLERANCE
        ):
            return False
        if (
            abs(observed_row["standard_deviation"] - expected_row["standard_deviation"])
            > _CONTINUOUS_STANDARD_DEVIATION_TOLERANCE
        ):
            return False
    expected_tip_covariance = expected["tip_covariance"]
    observed_tip_covariance = observed["tip_covariance"]
    for label, expected_row in expected_tip_covariance.items():
        observed_row = observed_tip_covariance[label]
        if (
            abs(observed_row["covariance"] - expected_row["covariance"])
            > _CONTINUOUS_COVARIANCE_TOLERANCE
        ):
            return False
        if (
            abs(observed_row["correlation"] - expected_row["correlation"])
            > _CONTINUOUS_CORRELATION_TOLERANCE
        ):
            return False
    return True


def _discrete_tip_state_frequency_metrics(report) -> dict[str, object]:
    return {
        "tip_state_frequency": {
            row.label: row.mean_value
            for row in report.rows
            if row.row_kind == "tip_state_frequency"
        }
    }


def _discrete_tip_state_frequency_metrics_match(
    *,
    expected: dict[str, float],
    observed: dict[str, float],
) -> bool:
    for label, expected_value in expected.items():
        if (
            abs(observed[label] - expected_value)
            > _DISCRETE_TIP_STATE_FREQUENCY_TOLERANCE
        ):
            return False
    return True


def _tree_fixture_path(name: str) -> Path:
    return Path(__file__).resolve().parents[4] / "tests" / "fixtures" / "trees" / name
