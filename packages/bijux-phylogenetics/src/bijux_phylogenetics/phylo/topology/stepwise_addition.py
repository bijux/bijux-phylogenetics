from __future__ import annotations

from dataclasses import dataclass

from .tree import TreeNode, descendant_taxa

STEPWISE_ADDITION_ROOT_BRANCH_ID = "root"
_SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS = frozenset({"minimize", "maximize"})


@dataclass(frozen=True, slots=True)
class StepwiseAdditionEdgeCandidate:
    """One branch on which a new taxon may be inserted."""

    target_node_id: str | None
    branch_id: str
    descendant_taxa: tuple[str, ...]


def validate_stepwise_addition_taxa(taxa: list[str]) -> list[str]:
    """Require at least two distinct non-empty taxa while preserving insertion order."""
    if len(taxa) < 2:
        raise ValueError("stepwise addition requires at least two taxa")
    seen: set[str] = set()
    duplicates: list[str] = []
    for taxon in taxa:
        if not taxon.strip():
            raise ValueError("stepwise addition does not allow blank taxon labels")
        if taxon in seen and taxon not in duplicates:
            duplicates.append(taxon)
        seen.add(taxon)
    if duplicates:
        raise ValueError(
            "stepwise addition requires distinct taxa; duplicates: "
            + ", ".join(duplicates)
        )
    return list(taxa)


def validate_stepwise_objective_direction(objective_direction: str) -> str:
    """Validate whether the score objective is minimized or maximized."""
    normalized_direction = objective_direction.strip().lower()
    if normalized_direction not in _SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS:
        raise ValueError(
            "objective_direction must be one of "
            + ", ".join(sorted(_SUPPORTED_STEPWISE_OBJECTIVE_DIRECTIONS))
        )
    return normalized_direction


def stepwise_addition_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort insertion edges deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)
