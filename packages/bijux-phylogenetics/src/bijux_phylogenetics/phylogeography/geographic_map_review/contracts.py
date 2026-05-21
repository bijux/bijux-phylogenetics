from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class GeographicMapMarkerRow:
    """One geographic marker plotted on the owned map surface."""

    marker_id: str
    label: str
    marker_kind: str
    latitude: float
    longitude: float
    state_label: str
    descendant_taxa: list[str]
    confidence: float
    uncertainty_km: float | None
    is_tip: bool
    is_root: bool
    active_line_count: int


@dataclass(frozen=True, slots=True)
class GeographicMapLineRow:
    """One geographic line plotted on the owned map surface."""

    line_id: str
    line_kind: str
    source_label: str
    target_label: str
    source_latitude: float
    source_longitude: float
    target_latitude: float
    target_longitude: float
    child_descendant_taxa: list[str]
    support: float
    midpoint_depth: float | None
    branch_length: float | None
    distance_km: float
    state_transition: str
    flag_codes: list[str]
    visible: bool


@dataclass(frozen=True, slots=True)
class GeographicMapExclusionRow:
    """One record omitted from the owned map review surface."""

    subject_id: str
    subject_kind: str
    raw_left: str
    raw_right: str
    reason: str
    note: str


@dataclass(frozen=True, slots=True)
class GeographicMapSummary:
    """Reviewer-facing summary for one geographic map review."""

    mode: str
    model: str
    analyzed_taxon_count: int
    excluded_record_count: int
    tip_marker_count: int
    internal_marker_count: int
    root_marker_count: int
    line_count: int
    visible_line_count: int
    tree_depth: float
    time_filter_applied: bool
    minimum_midpoint_depth_filter: float | None
    maximum_midpoint_depth_filter: float | None
    earliest_visible_midpoint_depth: float | None
    latest_visible_midpoint_depth: float | None
    warning_count: int


@dataclass(slots=True)
class GeographicMapReport:
    """Owned HTML map review surface for geographic reconstruction outputs."""

    mode: str
    tree_path: Path
    table_path: Path
    summary: GeographicMapSummary
    marker_rows: list[GeographicMapMarkerRow]
    line_rows: list[GeographicMapLineRow]
    exclusion_rows: list[GeographicMapExclusionRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class GeographicMapArtifact:
    """One rendered geographic HTML map artifact."""

    output_path: Path
    format: str
    width: int
    height: int
    visible_line_count: int
