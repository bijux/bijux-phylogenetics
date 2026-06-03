from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment


@dataclass(slots=True)
class BayesianPosteriorReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    posterior_tree_path: Path
    trace_path: Path
    kept_tree_count: int
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class CalibrationAuditReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    calibration_path: Path
    tip_dates_path: Path | None
    invalid_calibration_count: int
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class BayesianRunComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    left_tree_set_path: Path
    right_tree_set_path: Path
    left_trace_path: Path
    right_trace_path: Path
    trace_kind: str
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class BayesianDiagnosticsReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    posterior_tree_path: Path
    primary_log_path: Path
    chain_count: int
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class BayesianMlComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    ml_tree_path: Path
    posterior_tree_path: Path
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TimeTreeReadinessReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    warning_count: int
    method_tier: MethodTierAssessment
    machine_manifest: dict[str, object]
