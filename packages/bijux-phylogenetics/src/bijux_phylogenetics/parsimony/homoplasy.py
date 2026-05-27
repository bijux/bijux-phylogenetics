from __future__ import annotations

from collections import Counter

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import FitchCharacterMatrix
from .wagner import resolve_wagner_character_state_order


def character_kind(
    matrix: FitchCharacterMatrix,
    character_id: str,
    observed_states: list[str],
) -> str:
    """Classify one character by its parsimony-information content."""
    if len(observed_states) == 1:
        return "constant"
    counts = Counter(
        matrix.states_by_taxon[taxon][character_id] for taxon in matrix.states_by_taxon
    )
    repeated_state_count = sum(1 for count in counts.values() if count >= 2)
    if repeated_state_count >= 2:
        return "parsimony-informative"
    return "parsimony-uninformative"


def minimum_possible_steps(
    method: str,
    observed_states: list[str],
    *,
    character_id: str,
    state_order: list[str] | None,
) -> float:
    """Compute the minimum theoretical step count for one supported method."""
    if method in {"fitch", "acctran", "deltran"}:
        return float(max(len(observed_states) - 1, 0))
    if method == "wagner":
        resolved_state_order = resolve_wagner_character_state_order(
            character_id=character_id,
            observed_states=observed_states,
            state_order=state_order,
        )
        coordinates = resolved_state_order.coordinates
        return float(
            max(coordinates[state] for state in observed_states)
            - min(coordinates[state] for state in observed_states)
        )
    if method in {"dollo", "camin-sokal"}:
        return 0.0 if observed_states == ["0"] else 1.0
    raise ParsimonyAnalysisError(
        "minimum-step calculation requires one supported parsimony method",
        code="parsimony_homoplasy_method_unknown",
        details={"method": method},
    )


def unordered_maximum_possible_steps(
    matrix: FitchCharacterMatrix,
    character_id: str,
) -> float:
    """Compute the unordered maximum possible step count for one character."""
    counts = Counter(
        matrix.states_by_taxon[taxon][character_id] for taxon in matrix.states_by_taxon
    )
    if not counts:
        raise ParsimonyAnalysisError(
            "maximum-step calculation requires at least one observed state",
            code="parsimony_homoplasy_empty_character",
            details={"character_id": character_id},
        )
    return float(matrix.taxon_count - max(counts.values()))
