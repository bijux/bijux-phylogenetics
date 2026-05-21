from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..likelihood import GeographicExcludedTaxonRow


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasSummary:
    """One summary row for region-sampling bias review on one tree."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    weighting_mode: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_region_count: int
    region_dominated: bool
    dominant_region: str
    dominant_region_fraction: float
    weighted_region_dominated: bool
    weighted_dominant_region: str
    weighted_dominant_region_fraction: float
    root_region_unweighted: str
    root_region_weighted: str
    root_region_changed: bool
    compared_internal_node_count: int
    changed_internal_node_count: int
    compared_transition_count: int
    changed_transition_count: int
    warning_count: int


@dataclass(frozen=True, slots=True)
class GeographicSamplingCountRow:
    """One observed-region sample-count and weight row."""

    region: str
    sample_count: int
    sample_fraction: float
    applied_weight: float
    weighted_sample_count: float
    weighted_sample_fraction: float
    dominant_unweighted: bool
    dominant_weighted: bool


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasNodeRow:
    """One weighted-versus-unweighted node region comparison row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    is_root: bool
    unweighted_region: str
    weighted_region: str
    unweighted_confidence: float
    weighted_confidence: float
    confidence_delta: float
    changed: bool
    unweighted_region_probabilities: dict[str, float]
    weighted_region_probabilities: dict[str, float]


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasTransitionRow:
    """One weighted-versus-unweighted branch transition comparison row."""

    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    unweighted_source_region: str
    unweighted_target_region: str
    weighted_source_region: str
    weighted_target_region: str
    unweighted_transition: str
    weighted_transition: str
    unweighted_changed: bool
    weighted_changed: bool
    changed_by_weighting: bool
    unweighted_support: float
    weighted_support: float


@dataclass(slots=True)
class GeographicSamplingBiasReport:
    """Owned review surface for region-sampling bias correction on one tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    weighting_mode: str
    summary: GeographicSamplingBiasSummary
    count_rows: list[GeographicSamplingCountRow]
    node_rows: list[GeographicSamplingBiasNodeRow]
    transition_rows: list[GeographicSamplingBiasTransitionRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]
