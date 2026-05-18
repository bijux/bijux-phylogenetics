from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class BootstrapSupportNode:
    node: str
    descendant_taxa: list[str]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class BootstrapSupportSummaryReport:
    tree_path: Path
    internal_node_count: int
    supported_node_count: int
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    support_histogram: dict[str, int]
    nodes: list[BootstrapSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class FastTreeSupportNode:
    node: str
    descendant_taxa: list[str]
    local_support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class FastTreeSupportSummaryReport:
    tree_path: Path
    internal_node_count: int
    annotated_node_count: int
    minimum_local_support: float | None
    maximum_local_support: float | None
    median_local_support: float | None
    weakly_supported_clade_count: int
    support_histogram: dict[str, int]
    approximate_method: bool
    support_label_kind: str
    support_scale: str
    nodes: list[FastTreeSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class WeakBackboneReport:
    tree_path: Path
    threshold: float
    evaluated_backbone_node_count: int
    weak_backbone_node_count: int
    weak_nodes: list[BootstrapSupportNode]
    warnings: list[str]


@dataclass(slots=True)
class ShAlrtSupportNode:
    node: str
    descendant_taxa: list[str]
    sh_alrt_support: float | None
    sh_alrt_support_fraction: float | None
    ufboot_support: float | None
    ufboot_support_fraction: float | None
    is_backbone: bool
    sh_alrt_strong: bool
    ufboot_strong: bool
    conflicting_support_signal: bool
    support_agreement: str


@dataclass(slots=True)
class ShAlrtSupportSummaryReport:
    tree_path: Path
    internal_node_count: int
    annotated_node_count: int
    fully_scored_node_count: int
    minimum_sh_alrt_support: float | None
    maximum_sh_alrt_support: float | None
    minimum_ufboot_support: float | None
    maximum_ufboot_support: float | None
    weak_sh_alrt_clade_count: int
    weak_ufboot_clade_count: int
    conflicting_support_signal_count: int
    nodes: list[ShAlrtSupportNode]
    warnings: list[str]
