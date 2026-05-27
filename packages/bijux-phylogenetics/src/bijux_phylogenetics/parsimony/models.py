from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class FitchCharacterMatrix:
    """Taxon-by-character matrix normalized for unordered Fitch scoring."""

    matrix_path: Path | None
    taxon_column: str
    character_ids: list[str]
    states_by_taxon: dict[str, dict[str, str]]

    @property
    def taxon_count(self) -> int:
        return len(self.states_by_taxon)

    @property
    def character_count(self) -> int:
        return len(self.character_ids)


@dataclass(frozen=True, slots=True)
class FitchCharacterScore:
    """Per-character unordered Fitch tree-length row."""

    character_id: str
    step_count: int
    observed_states: list[str]


@dataclass(frozen=True, slots=True)
class FitchNodeStateSet:
    """Per-node unordered Fitch candidate-state row for one character."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state_set: list[str]
    is_tip: bool
    observed_state: str | None


@dataclass(frozen=True, slots=True)
class FitchScoreReport:
    """Complete unordered Fitch scoring report over one tree and matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_steps: int
    step_rows: list[FitchCharacterScore]
    node_state_rows: list[FitchNodeStateSet]


ParsimonyCharacterMatrix = FitchCharacterMatrix


@dataclass(frozen=True, slots=True)
class WagnerCharacterScore:
    """Per-character ordered Wagner weighted tree-length row."""

    character_id: str
    weighted_step_count: int
    observed_states: list[str]
    state_order: list[str]
    optimal_root_states: list[str]


@dataclass(frozen=True, slots=True)
class WagnerNodeCost:
    """Per-node ordered Wagner cost row for one candidate state."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state: str
    cost: int
    is_optimal_state: bool


@dataclass(frozen=True, slots=True)
class WagnerScoreReport:
    """Complete ordered Wagner scoring report over one tree and matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_cost: int
    step_rows: list[WagnerCharacterScore]
    node_cost_rows: list[WagnerNodeCost]


@dataclass(frozen=True, slots=True)
class SankoffCostMatrix:
    """Validated state-to-state Sankoff transition-cost matrix."""

    matrix_path: Path | None
    states: list[str]
    costs: dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class SankoffCharacterScore:
    """Per-character Sankoff minimum-cost row."""

    character_id: str
    minimum_cost: float
    observed_states: list[str]
    matrix_states: list[str]


@dataclass(frozen=True, slots=True)
class SankoffNodeCost:
    """Per-node Sankoff cost-vector row for one candidate state."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state: str
    cost: float
    is_optimal_state: bool


@dataclass(frozen=True, slots=True)
class SankoffNodeSelection:
    """Per-node Sankoff optimal-state selection row."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    optimal_states: list[str]
    tie_states: list[str]


@dataclass(frozen=True, slots=True)
class SankoffScoreReport:
    """Complete Sankoff scoring report over one tree, matrix, and cost matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_cost: float
    step_rows: list[SankoffCharacterScore]
    node_cost_rows: list[SankoffNodeCost]
    node_selection_rows: list[SankoffNodeSelection]


@dataclass(frozen=True, slots=True)
class DolloCharacterScore:
    """Per-character Dollo summary row."""

    character_id: str
    derived_taxon_count: int
    gain_node: str | None
    gain_node_name: str | None
    gain_descendant_taxa: list[str]
    total_losses: int
    impossible_state_warning: str | None


@dataclass(frozen=True, slots=True)
class DolloBranchChange:
    """Per-branch Dollo change row."""

    character_id: str
    change_kind: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]


@dataclass(frozen=True, slots=True)
class DolloScoreReport:
    """Complete Dollo scoring report over one tree and binary matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_gains: int
    total_losses: int
    step_rows: list[DolloCharacterScore]
    branch_change_rows: list[DolloBranchChange]
