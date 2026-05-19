from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.parity.geiger.boundary_warning_registry import (
    GeigerBoundaryWarningRow,
)
from bijux_phylogenetics.parity.geiger.likelihood_policy import (
    GeigerLikelihoodPolicyRow,
)
from bijux_phylogenetics.parity.geiger.model_confidence import GeigerModelConfidenceRow
from bijux_phylogenetics.parity.geiger.optimizer_triage import GeigerOptimizerTriageRow
from bijux_phylogenetics.parity.geiger.parameterization_registry import (
    GeigerParameterizationRegistryRow,
)


@dataclass(frozen=True, slots=True)
class GeigerParityObservation:
    """One live parity comparison between Bijux and `geiger`."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixtures: tuple[Path, ...]
    model_name: str
    optimizer_settings: dict[str, object] | None
    tolerance: float
    r_version: str | None
    geiger_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None
    reference_rows: list[dict[str, object]] | None
    bijux_rows: list[dict[str, object]] | None
    reference_error: dict[str, object] | None
    bijux_error: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class GeigerParitySummaryRow:
    """One function-level summary across governed `geiger` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class GeigerParityReport:
    """Aggregate report for governed live `geiger` parity cases."""

    observations: list[GeigerParityObservation]
    optimizer_triage_rows: list[GeigerOptimizerTriageRow]
    boundary_warning_rows: list[GeigerBoundaryWarningRow]
    likelihood_policy_rows: list[GeigerLikelihoodPolicyRow]
    model_confidence_rows: list[GeigerModelConfidenceRow]
    parameterization_registry_rows: list[GeigerParameterizationRegistryRow]
    summary_rows: list[GeigerParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]
