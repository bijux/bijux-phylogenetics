from __future__ import annotations

from ..numeric import _quantile
from .models import (
    StochasticMapBranchOccupancyRow,
    StochasticMapBranchTransitionCountRow,
    StochasticMapCollectionReport,
    StochasticMapReplicate,
    StochasticMapStateTimeRow,
    StochasticMapSummaryReport,
    StochasticMapSummaryRow,
    StochasticMapTransitionCountMatrixRow,
    StochasticMapTransitionCountReport,
)


def _summarize_stochastic_map_replicates(
    replicates: list[StochasticMapReplicate],
    *,
    simulation_failure_count: int,
    expected_transitions: list[str] | None = None,
) -> StochasticMapSummaryReport:
    total_counts = sorted(
        float(replicate.total_transition_count) for replicate in replicates
    )
    transition_names = sorted(
        set(expected_transitions or [])
        | {
            transition
            for replicate in replicates
            for transition in replicate.transition_counts
        }
    )
    rows: list[StochasticMapSummaryRow] = []
    for transition in transition_names:
        values = [
            replicate.transition_counts.get(transition, 0) for replicate in replicates
        ]
        sorted_values = sorted(float(value) for value in values)
        rows.append(
            StochasticMapSummaryRow(
                transition=transition,
                mean_count=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_count=min(values, default=0),
                maximum_count=max(values, default=0),
                presence_fraction=float(
                    format(
                        sum(1 for value in values if value > 0) / max(len(values), 1),
                        ".15g",
                    )
                ),
            )
        )
    state_names = sorted(
        {state for replicate in replicates for state in replicate.state_time_totals}
    )
    state_time_rows: list[StochasticMapStateTimeRow] = []
    for state in state_names:
        values = [
            float(replicate.state_time_totals.get(state, 0.0))
            for replicate in replicates
        ]
        sorted_values = sorted(values)
        state_time_rows.append(
            StochasticMapStateTimeRow(
                state=state,
                mean_time=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_time=min(values, default=0.0),
                maximum_time=max(values, default=0.0),
            )
        )
    branch_lookup: dict[tuple[int, str, str, float], list[dict[str, float]]] = {}
    for replicate in replicates:
        for history in replicate.branch_histories:
            key = (
                history.branch_index,
                history.parent_node,
                history.child_node,
                float(history.branch_length),
            )
            branch_lookup.setdefault(key, [])
            state_times = dict.fromkeys(state_names, 0.0)
            for segment in history.segments:
                state_times[segment.state] = state_times.get(
                    segment.state, 0.0
                ) + float(segment.duration)
            branch_lookup[key].append(state_times)
    branch_occupancy_rows: list[StochasticMapBranchOccupancyRow] = []
    for (
        branch_index,
        parent_node,
        child_node,
        branch_length,
    ), state_times in sorted(
        branch_lookup.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][2],
        ),
    ):
        for state in state_names:
            values = [
                replicate_state_times.get(state, 0.0)
                for replicate_state_times in state_times
            ]
            sorted_values = sorted(float(value) for value in values)
            mean_time = float(format(sum(values) / max(len(values), 1), ".15g"))
            mean_fraction = 0.0
            if branch_length > 0.0:
                mean_fraction = float(format(mean_time / branch_length, ".15g"))
            branch_occupancy_rows.append(
                StochasticMapBranchOccupancyRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    state=state,
                    branch_length=branch_length,
                    mean_time=mean_time,
                    lower_95_interval=_quantile(sorted_values, 0.025),
                    upper_95_interval=_quantile(sorted_values, 0.975),
                    minimum_time=min(values, default=0.0),
                    maximum_time=max(values, default=0.0),
                    mean_fraction=mean_fraction,
                    presence_fraction=float(
                        format(
                            sum(1 for value in values if value > 0.0)
                            / max(len(values), 1),
                            ".15g",
                        )
                    ),
                )
            )
    warnings: list[str] = []
    if simulation_failure_count > 0:
        warnings.append(
            "one or more stochastic-map replicates failed to sample a branch history consistent with the conditioned endpoint states"
        )
    return StochasticMapSummaryReport(
        replicate_count=len(replicates),
        mean_total_transition_count=float(
            format(sum(total_counts) / max(len(total_counts), 1), ".15g")
        ),
        lower_95_total_transition_count=_quantile(total_counts, 0.025),
        upper_95_total_transition_count=_quantile(total_counts, 0.975),
        rows=rows,
        state_time_rows=state_time_rows,
        branch_occupancy_rows=branch_occupancy_rows,
        simulation_failure_count=simulation_failure_count,
        warnings=warnings,
    )


def _stochastic_map_warning_union(*warning_lists: list[str]) -> list[str]:
    merged: list[str] = []
    for warnings in warning_lists:
        for warning in warnings:
            if warning not in merged:
                merged.append(warning)
    return merged


def summarize_discrete_stochastic_maps(
    report: StochasticMapCollectionReport,
) -> StochasticMapSummaryReport:
    """Summarize one stochastic-map collection without regenerating maps."""
    return _summarize_stochastic_map_replicates(
        report.maps,
        simulation_failure_count=len(report.failures),
        expected_transitions=report.fit_audit.allowed_transitions,
    )


def count_discrete_stochastic_map_transitions(
    report: StochasticMapCollectionReport,
) -> StochasticMapTransitionCountReport:
    """Count directional transitions across one stochastic-map collection."""
    summary = summarize_discrete_stochastic_maps(report)
    state_order = list(report.fit_audit.state_order)
    if not state_order:
        state_order = sorted(
            {
                state
                for replicate in report.maps
                for state in replicate.state_time_totals
            }
        )
    transition_order = [
        f"{source_state}->{target_state}"
        for source_state in state_order
        for target_state in state_order
    ]
    aggregate_lookup = {row.transition: row for row in summary.rows}
    aggregate_rows = [
        aggregate_lookup.get(
            transition,
            StochasticMapSummaryRow(
                transition=transition,
                mean_count=0.0,
                lower_95_interval=0.0,
                upper_95_interval=0.0,
                minimum_count=0,
                maximum_count=0,
                presence_fraction=0.0,
            ),
        )
        for transition in transition_order
    ]
    matrix_rows = [
        StochasticMapTransitionCountMatrixRow(
            replicate_index=replicate.replicate_index,
            total_transition_count=replicate.total_transition_count,
            transition_counts={
                transition: int(replicate.transition_counts.get(transition, 0))
                for transition in transition_order
            },
        )
        for replicate in report.maps
    ]
    branch_keys = [
        (
            history.branch_index,
            history.parent_node,
            history.child_node,
        )
        for history in report.maps[0].branch_histories
    ]
    branch_transition_values: dict[tuple[int, str, str, str], list[int]] = {
        (*branch_key, transition): []
        for branch_key in branch_keys
        for transition in transition_order
    }
    for replicate in report.maps:
        replicate_branch_counts: dict[tuple[int, str, str], dict[str, int]] = {}
        for history in replicate.branch_histories:
            transition_counts = dict.fromkeys(transition_order, 0)
            inferred_transitions = [
                f"{event.source_state}->{event.target_state}"
                for event in history.events
            ]
            if not inferred_transitions and len(history.segments) > 1:
                inferred_transitions = [
                    f"{left.state}->{right.state}"
                    for left, right in zip(
                        history.segments,
                        history.segments[1:],
                        strict=False,
                    )
                    if left.state != right.state
                ]
            for transition in inferred_transitions:
                transition_counts[transition] = transition_counts.get(transition, 0) + 1
            replicate_branch_counts[
                (
                    history.branch_index,
                    history.parent_node,
                    history.child_node,
                )
            ] = transition_counts
        for branch_key in branch_keys:
            transition_counts = replicate_branch_counts.get(
                branch_key,
                dict.fromkeys(transition_order, 0),
            )
            for transition in transition_order:
                branch_transition_values[(*branch_key, transition)].append(
                    int(transition_counts.get(transition, 0))
                )
    branch_rows: list[StochasticMapBranchTransitionCountRow] = []
    for (
        branch_index,
        parent_node,
        child_node,
        transition,
    ), values in sorted(
        branch_transition_values.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][2],
            item[0][3],
        ),
    ):
        sorted_values = sorted(float(value) for value in values)
        branch_rows.append(
            StochasticMapBranchTransitionCountRow(
                branch_index=branch_index,
                parent_node=parent_node,
                child_node=child_node,
                transition=transition,
                mean_count=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_count=min(values, default=0),
                maximum_count=max(values, default=0),
                presence_fraction=float(
                    format(
                        sum(1 for value in values if value > 0) / max(len(values), 1),
                        ".15g",
                    )
                ),
            )
        )
    return StochasticMapTransitionCountReport(
        replicate_count=summary.replicate_count,
        mean_total_transition_count=summary.mean_total_transition_count,
        lower_95_total_transition_count=summary.lower_95_total_transition_count,
        upper_95_total_transition_count=summary.upper_95_total_transition_count,
        transition_order=transition_order,
        matrix_rows=matrix_rows,
        aggregate_rows=aggregate_rows,
        branch_rows=branch_rows,
        warnings=list(summary.warnings),
    )
