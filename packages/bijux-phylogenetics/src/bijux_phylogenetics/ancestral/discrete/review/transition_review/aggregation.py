from __future__ import annotations

import statistics

from .contracts import (
    AncestralTransitionBranchRow,
    AncestralTransitionCountRow,
    AncestralTransitionTreeRow,
    AncestralTransitionTreeSetCountRow,
)


def _summarize_transition_rows(
    branch_rows: list[AncestralTransitionBranchRow],
) -> list[AncestralTransitionCountRow]:
    grouped: dict[str, list[AncestralTransitionBranchRow]] = {}
    for row in branch_rows:
        if row.changed:
            grouped.setdefault(row.transition, []).append(row)
    transition_rows: list[AncestralTransitionCountRow] = []
    for transition in sorted(grouped):
        transition_branch_rows = grouped[transition]
        source_state, target_state = transition.split("->", maxsplit=1)
        certain_change_count = sum(
            row.certainty_class == "certain_change" for row in transition_branch_rows
        )
        uncertain_change_count = sum(
            row.certainty_class == "uncertain_change" for row in transition_branch_rows
        )
        transition_rows.append(
            AncestralTransitionCountRow(
                transition=transition,
                source_state=source_state,
                target_state=target_state,
                certain_change_count=certain_change_count,
                uncertain_change_count=uncertain_change_count,
                total_change_count=len(transition_branch_rows),
            )
        )
    return transition_rows


def _summarize_transition_rows_across_trees(
    *,
    tree_rows: list[AncestralTransitionTreeRow],
    per_tree_transition_counts: dict[str, list[AncestralTransitionCountRow]],
) -> list[AncestralTransitionTreeSetCountRow]:
    counts_by_transition: dict[str, dict[int, AncestralTransitionCountRow]] = {}
    tree_count = len(tree_rows)
    for tree_row in tree_rows:
        transition_rows = per_tree_transition_counts[str(tree_row.source_tree_index)]
        for transition_row in transition_rows:
            counts_by_transition.setdefault(transition_row.transition, {})[
                tree_row.source_tree_index
            ] = transition_row
    summary_rows: list[AncestralTransitionTreeSetCountRow] = []
    for transition in sorted(counts_by_transition):
        present_rows = counts_by_transition[transition]
        source_state, target_state = transition.split("->", maxsplit=1)
        certain_counts = [
            float(row.certain_change_count) for row in present_rows.values()
        ]
        uncertain_counts = [
            float(row.uncertain_change_count) for row in present_rows.values()
        ]
        total_counts = [float(row.total_change_count) for row in present_rows.values()]
        tree_presence_fraction = len(present_rows) / tree_count
        mean_uncertain_change_count = statistics.fmean(uncertain_counts)
        if tree_presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif mean_uncertain_change_count > 0.0:
            stability_class = "uncertainty_sensitive"
        else:
            stability_class = "stable"
        summary_rows.append(
            AncestralTransitionTreeSetCountRow(
                transition=transition,
                source_state=source_state,
                target_state=target_state,
                tree_presence_count=len(present_rows),
                tree_presence_fraction=tree_presence_fraction,
                mean_certain_change_count=statistics.fmean(certain_counts),
                mean_uncertain_change_count=mean_uncertain_change_count,
                mean_total_change_count=statistics.fmean(total_counts),
                minimum_total_change_count=int(min(total_counts)),
                maximum_total_change_count=int(max(total_counts)),
                lower_95_total_change_count=_empirical_quantile(total_counts, 0.025),
                upper_95_total_change_count=_empirical_quantile(total_counts, 0.975),
                stability_class=stability_class,
            )
        )
    return summary_rows


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")


def _empirical_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * probability
    lower_index = int(index)
    upper_index = min(lower_index + 1, len(ordered) - 1)
    if lower_index == upper_index:
        return ordered[lower_index]
    fraction = index - lower_index
    return ordered[lower_index] + (
        (ordered[upper_index] - ordered[lower_index]) * fraction
    )
