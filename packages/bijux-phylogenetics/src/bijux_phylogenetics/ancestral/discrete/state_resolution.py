from __future__ import annotations

from collections.abc import Mapping, Sequence


def clade_observed_state_counts(
    clade_taxa: Sequence[str],
    candidate_states: Sequence[str],
    observed_states_by_taxon: Mapping[str, str],
) -> dict[str, int]:
    """Count observed tip states that support a candidate clade state."""
    counts = dict.fromkeys(candidate_states, 0)
    for taxon in clade_taxa:
        observed_state = observed_states_by_taxon.get(taxon)
        if observed_state in counts:
            counts[observed_state] += 1
    return counts


def resolve_clade_consensus_state(
    *,
    clade_taxa: Sequence[str],
    candidate_states: Sequence[str],
    observed_states_by_taxon: Mapping[str, str],
    fallback_state: str,
) -> str:
    """Resolve one clade state without inventing conflict from ambiguous ties."""
    if len(candidate_states) == 1:
        return candidate_states[0]
    state_counts = clade_observed_state_counts(
        clade_taxa,
        candidate_states,
        observed_states_by_taxon,
    )
    best_count = max(state_counts.values(), default=0)
    if best_count <= 0:
        return fallback_state
    tied_states = sorted(
        state for state, count in state_counts.items() if count == best_count
    )
    return tied_states[0]
