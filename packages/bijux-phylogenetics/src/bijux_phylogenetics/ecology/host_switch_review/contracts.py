from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class HostStateNodeRow:
    """One internal-node host-state reconstruction row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_host: str
    host_probabilities: dict[str, float]
    confidence: float
    ambiguous: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class HostSwitchBranchRow:
    """One branchwise host-switch review row."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_host: str
    child_most_likely_host: str
    parent_host_set: list[str]
    child_host_set: list[str]
    overlapping_hosts: list[str]
    changed: bool
    transition: str
    certainty_class: str
    parent_confidence: float
    child_confidence: float
    transition_allowed: bool


@dataclass(frozen=True, slots=True)
class HostSwitchCountRow:
    """One aggregated host-switch count row."""

    transition: str
    source_host: str
    target_host: str
    transition_allowed: bool
    certain_switch_count: int
    uncertain_switch_count: int
    total_switch_count: int


@dataclass(frozen=True, slots=True)
class HostSwitchFitRow:
    """One fitted host-transition regime row."""

    constraint_mode: str
    model: str
    analyzed_taxon_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    root_host: str
    root_confidence: float


@dataclass(frozen=True, slots=True)
class UnsupportedHostSwitchClaimRow:
    """One unconstrained host-switch claim forbidden by the supplied transition policy."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    unconstrained_source_host: str
    unconstrained_target_host: str
    unconstrained_certainty_class: str
    constrained_source_host: str
    constrained_target_host: str
    constrained_certainty_class: str
    claim_resolved: bool


@dataclass(frozen=True, slots=True)
class HostSwitchExclusionRow:
    """One excluded row from the host-switching workflow."""

    taxon: str
    raw_host: str
    normalized_host: str | None
    reason: str
    note: str


@dataclass(frozen=True, slots=True)
class HostSwitchSummary:
    """Reviewer-facing summary for host-switching analysis."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analysis_constraint_mode: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_host_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    host_switch_count: int
    certain_host_switch_count: int
    uncertain_host_switch_count: int
    allowed_transition_count: int
    forbidden_transition_count: int
    constrained_log_likelihood: float | None
    unconstrained_log_likelihood: float
    constrained_aic: float | None
    unconstrained_aic: float
    preferred_constraint: str
    unsupported_switch_claim_count: int
    root_host: str
    root_confidence: float
    warning_count: int


@dataclass(slots=True)
class HostSwitchingReport:
    """Owned host-switching review surface on one parasite or pathogen tree."""

    tree_path: Path
    traits_path: Path
    constraint_path: Path | None
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analysis_constraint_mode: str
    summary: HostSwitchSummary
    node_rows: list[HostStateNodeRow]
    branch_rows: list[HostSwitchBranchRow]
    count_rows: list[HostSwitchCountRow]
    fit_rows: list[HostSwitchFitRow]
    unsupported_claim_rows: list[UnsupportedHostSwitchClaimRow]
    exclusion_rows: list[HostSwitchExclusionRow]
    warnings: list[str]
