from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class PhytoolsParityObservation:
    """One live parity comparison between Bijux and `phytools`."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    r_version: str | None
    phytools_version: str | None
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
class PhytoolsParitySummaryRow:
    """One function-level summary across governed `phytools` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class PhytoolsParityReport:
    """Aggregate report for governed live `phytools` parity cases."""

    observations: list[PhytoolsParityObservation]
    summary_rows: list[PhytoolsParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]
