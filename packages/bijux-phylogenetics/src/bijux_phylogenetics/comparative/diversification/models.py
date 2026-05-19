from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class TimeTreeValidationReport:
    tree_path: Path
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    tip_count: int
    root_age: float
    warnings: list[str]


@dataclass(slots=True)
class LineageThroughTimePoint:
    node: str
    time_before_present: float
    lineage_count: int
    event: str


@dataclass(slots=True)
class LineageThroughTimeReport:
    tree_path: Path
    tip_count: int
    root_age: float
    points: list[LineageThroughTimePoint]


@dataclass(slots=True)
class SamplingFractionIssue:
    taxon: str
    code: str
    raw_value: str
    message: str


@dataclass(slots=True)
class SamplingFractionReport:
    tree_path: Path
    metadata_path: Path
    taxon_column: str
    sampling_column: str | None
    complete: bool
    matched_taxa: list[str]
    missing_taxa: list[str]
    invalid_rows: list[SamplingFractionIssue]
    sampling_fraction: float | None
    heterogeneous_values: bool
    warnings: list[str]


@dataclass(slots=True)
class DiversificationRateReport:
    tree_path: Path
    model: str
    crown_age: float
    observed_tip_count: int
    corrected_tip_count: float
    sampling_fraction: float
    birth_rate: float
    death_rate: float
    net_diversification_rate: float
    relative_extinction: float
    likelihood_kind: str
    log_likelihood: float
    aic: float
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationGammaStatisticReport:
    tree_path: Path
    tip_count: int
    rooted: bool
    ultrametric: bool
    bifurcating: bool
    root_age: float
    branching_time_count: int
    interval_count: int
    minimum_branching_time: float
    maximum_branching_time: float
    gamma_statistic: float
    sampling_fraction: float | None
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationModelComparisonRow:
    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    sampling_fraction: float
    net_diversification_rate: float
    relative_extinction: float


@dataclass(slots=True)
class DiversificationModelComparisonReport:
    tree_path: Path
    metadata_path: Path | None
    better_model: str
    rows: list[DiversificationModelComparisonRow]


@dataclass(slots=True)
class CladeDiversificationObservation:
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    tip_count: int
    crown_age: float
    diversification_rate: float
    z_score: float
    classification: str


@dataclass(slots=True)
class CladeDiversificationScanReport:
    tree_path: Path
    model: str
    global_rate: float
    observations: list[CladeDiversificationObservation]
    high_diversification_clades: list[CladeDiversificationObservation]
    low_diversification_clades: list[CladeDiversificationObservation]
    warnings: list[str]


@dataclass(slots=True)
class TraitDependentDiversificationState:
    state: str
    taxon_count: int
    taxa: list[str]
    monophyletic: bool
    crown_age: float | None
    diversification_rate: float | None
    warnings: list[str]


@dataclass(slots=True)
class TraitDependentDiversificationReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    observed_states: list[str]
    states: list[TraitDependentDiversificationState]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    machine_manifest: dict[str, object]
    methods_summary_text: str
    methods_summary_warning_count: int
    methods_summary_path: Path | None
    report: DiversificationMethodReport


@dataclass(slots=True)
class DiversificationMethodReport:
    tree_path: Path
    metadata_path: Path | None
    traits_path: Path | None
    taxon_column: str | None
    sampling_column: str | None
    estimate_model: str
    clade_model: str
    clade_min_tip_count: int
    validation: TimeTreeValidationReport
    lineage: LineageThroughTimeReport
    gamma_statistic: DiversificationGammaStatisticReport
    primary_estimate: DiversificationRateReport
    model_comparison: DiversificationModelComparisonReport
    clade_scan: CladeDiversificationScanReport
    sampling_report: SamplingFractionReport | None
    trait_report: TraitDependentDiversificationReport | None


@dataclass(slots=True)
class DiversificationMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    better_model: str
    sampling_metadata_complete: bool | None
    clade_observation_count: int
    text: str
    report: DiversificationMethodReport


@dataclass(slots=True)
class MedusaExclusionReport:
    tree_path: Path
    metadata_path: Path | None
    validation: TimeTreeValidationReport
    sampling_report: SamplingFractionReport | None
    supported_surfaces: list[str]
    missing_surfaces: list[str]
    exclusion_code: str
    exclusion_reason: str
    warnings: list[str]


@dataclass(slots=True)
class GeigerBirthDeathExclusionReport:
    tree_path: Path
    metadata_path: Path | None
    validation: TimeTreeValidationReport
    sampling_report: SamplingFractionReport | None
    geiger_reference_surface: str
    geiger_reference_arguments: list[str]
    owned_surface: str
    exclusion_code: str
    exclusion_reason: str
    warnings: list[str]
