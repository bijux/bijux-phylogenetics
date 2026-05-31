from __future__ import annotations

from collections import defaultdict

from .dna import DNA_STATE_ORDER
from .models import (
    MarginalAncestralSequenceProbabilityReport,
    MarginalAncestralSiteSummaryRow,
)


def summarize_marginal_ancestral_sites(
    report: MarginalAncestralSequenceProbabilityReport,
) -> list[MarginalAncestralSiteSummaryRow]:
    """Collapse full state rows into one posterior summary per node and site."""
    grouped_rows: dict[
        tuple[str, int],
        dict[str, object],
    ] = defaultdict(dict)
    for row in report.posterior_rows:
        key = (row.node_id, row.site_position)
        bucket = grouped_rows[key]
        bucket.setdefault("node_id", row.node_id)
        bucket.setdefault("node_name", row.node_name)
        bucket.setdefault("descendant_taxa", row.descendant_taxa)
        bucket.setdefault("pattern_id", row.pattern_id)
        bucket.setdefault("site_position", row.site_position)
        posterior_by_state = bucket.setdefault("posterior_by_state", {})
        if not isinstance(posterior_by_state, dict):
            raise TypeError("posterior_by_state grouping bucket must be a dictionary")
        posterior_by_state[row.state] = row.posterior_probability

    summary_rows: list[MarginalAncestralSiteSummaryRow] = []
    for key in sorted(grouped_rows, key=lambda item: (item[0], item[1])):
        bucket = grouped_rows[key]
        posterior_by_state = bucket["posterior_by_state"]
        if not isinstance(posterior_by_state, dict):
            raise TypeError("posterior_by_state grouping bucket must be a dictionary")
        _validate_complete_state_set(posterior_by_state)
        most_likely_state, max_posterior_probability = max(
            ((state, float(posterior_by_state[state])) for state in DNA_STATE_ORDER),
            key=lambda item: (item[1], -DNA_STATE_ORDER.index(item[0])),
        )
        summary_rows.append(
            MarginalAncestralSiteSummaryRow(
                node_id=str(bucket["node_id"]),
                node_name=(
                    None if bucket["node_name"] is None else str(bucket["node_name"])
                ),
                descendant_taxa=list(bucket["descendant_taxa"]),
                pattern_id=str(bucket["pattern_id"]),
                site_position=int(bucket["site_position"]),
                most_likely_state=most_likely_state,
                max_posterior_probability=max_posterior_probability,
                posterior_probability_a=float(posterior_by_state["A"]),
                posterior_probability_c=float(posterior_by_state["C"]),
                posterior_probability_g=float(posterior_by_state["G"]),
                posterior_probability_t=float(posterior_by_state["T"]),
            )
        )
    return summary_rows


def _validate_complete_state_set(posterior_by_state: dict[str, float]) -> None:
    missing_states = [
        state for state in DNA_STATE_ORDER if state not in posterior_by_state
    ]
    if missing_states:
        raise ValueError(
            "marginal ancestral posterior grouping requires complete DNA state rows; "
            f"missing {', '.join(missing_states)}"
        )
