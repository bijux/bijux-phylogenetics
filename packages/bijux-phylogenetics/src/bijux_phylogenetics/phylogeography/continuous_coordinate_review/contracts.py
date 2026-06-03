from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class CoordinateEstimateRow:
    """One reviewed coordinate estimate for a tip or internal node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    latitude: float
    longitude: float
    latitude_standard_error: float
    longitude_standard_error: float
    radial_standard_error_km: float
    lower_95_latitude: float
    upper_95_latitude: float
    lower_95_longitude: float
    upper_95_longitude: float
    confidence: float
    unstable: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class CoordinateMovementBranchRow:
    """One branchwise continuous movement review row."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_latitude: float
    parent_longitude: float
    child_latitude: float
    child_longitude: float
    great_circle_km: float
    branch_rate_km_per_unit: float | None
    support: float
    impossible_jump: bool
    outlier_jump: bool
    flag_codes: list[str]


@dataclass(frozen=True, slots=True)
class CoordinateMovementOutlierRow:
    """One flagged branchwise movement row."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    great_circle_km: float
    branch_rate_km_per_unit: float | None
    median_distance_km: float
    distance_threshold_km: float
    median_rate_km_per_unit: float | None
    rate_threshold_km: float | None
    impossible_jump: bool
    outlier_jump: bool
    flag_codes: list[str]


@dataclass(frozen=True, slots=True)
class CoordinateMovementExclusionRow:
    """One excluded coordinate row from continuous phylogeography."""

    taxon: str
    raw_latitude: str
    raw_longitude: str
    reason: str
    note: str


@dataclass(frozen=True, slots=True)
class CoordinateMovementSummary:
    """Reviewer-facing summary for one continuous phylogeographic reconstruction."""

    taxon_column: str
    latitude_column: str
    longitude_column: str
    model: str
    alpha: float
    analyzed_taxon_count: int
    excluded_taxon_count: int
    internal_node_count: int
    weak_support_node_count: int
    outlier_jump_count: int
    impossible_jump_count: int
    flagged_branch_count: int
    maximum_jump_km: float
    root_latitude: float
    root_longitude: float
    root_radial_standard_error_km: float
    warning_count: int


@dataclass(slots=True)
class PhylogeographicCoordinateReport:
    """Owned continuous-coordinate phylogeographic review surface."""

    tree_path: Path
    table_path: Path
    taxon_column: str
    latitude_column: str
    longitude_column: str
    model: str
    alpha: float
    summary: CoordinateMovementSummary
    estimate_rows: list[CoordinateEstimateRow]
    branch_rows: list[CoordinateMovementBranchRow]
    outlier_rows: list[CoordinateMovementOutlierRow]
    exclusion_rows: list[CoordinateMovementExclusionRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class CoordinateMovementVisualization:
    """One rendered coordinate-space movement visualization artifact."""

    output_path: Path
    format: str
    width: int
    height: int
    highlighted_branch_count: int
