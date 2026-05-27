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
