from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .models import (
    StochasticMapBranchHistory,
    StochasticMapBranchOccupancyRow,
    StochasticMapCollectionReport,
    StochasticMapDensityReport,
    StochasticMapModelFitAudit,
    StochasticMapReplicate,
    StochasticMapSimulationFailure,
    StochasticMapStateSegment,
    StochasticMapStateTimeRow,
    StochasticMapSummaryReport,
    StochasticMapSummaryRow,
    StochasticMapTransitionCountReport,
    StochasticMapTransitionEvent,
)


def write_stochastic_map_summary_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one transition-by-transition stochastic-map uncertainty table."""
    rows = [
        {
            "transition": row.transition,
            "mean_count": row.mean_count,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_count": row.minimum_count,
            "maximum_count": row.maximum_count,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "mean_count",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_count",
            "maximum_count",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_transition_count_matrix(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one countSimmap-style transition matrix with one row per replicate."""
    columns = ["replicate_index", "total_transition_count", *report.transition_order]
    rows = [
        {
            "replicate_index": row.replicate_index,
            "total_transition_count": row.total_transition_count,
            **{
                transition: row.transition_counts.get(transition, 0)
                for transition in report.transition_order
            },
        }
        for row in report.matrix_rows
    ]
    return write_taxon_rows(path, columns=columns, rows=rows)


def write_stochastic_map_aggregate_transition_matrix(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one mean transition matrix aggregated over a stochastic-map collection."""
    source_states = sorted(
        {
            transition.split("->", 1)[0]
            for transition in report.transition_order
            if "->" in transition
        }
    )
    target_states = sorted(
        {
            transition.split("->", 1)[1]
            for transition in report.transition_order
            if "->" in transition
        }
    )
    mean_lookup = {row.transition: row.mean_count for row in report.aggregate_rows}
    rows = [
        {
            "source_state": source_state,
            **{
                target_state: mean_lookup.get(f"{source_state}->{target_state}", 0.0)
                for target_state in target_states
            },
        }
        for source_state in source_states
    ]
    return write_taxon_rows(
        path,
        columns=["source_state", *target_states],
        rows=rows,
    )


def write_stochastic_map_branch_transition_count_table(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one per-branch transition-count summary table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "transition": row.transition,
            "mean_count": row.mean_count,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_count": row.minimum_count,
            "maximum_count": row.maximum_count,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "transition",
            "mean_count",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_count",
            "maximum_count",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_state_time_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one per-state time-in-state summary table for a stochastic-map collection."""
    rows = [
        {
            "state": row.state,
            "mean_time": row.mean_time,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_time": row.minimum_time,
            "maximum_time": row.maximum_time,
        }
        for row in report.state_time_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "state",
            "mean_time",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_time",
            "maximum_time",
        ],
        rows=rows,
    )


def write_stochastic_map_branch_occupancy_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one per-branch state-occupancy summary table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "state": row.state,
            "branch_length": row.branch_length,
            "mean_time": row.mean_time,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_time": row.minimum_time,
            "maximum_time": row.maximum_time,
            "mean_fraction": row.mean_fraction,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_occupancy_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "branch_length",
            "mean_time",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_time",
            "maximum_time",
            "mean_fraction",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_branch_probability_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one per-branch state-probability table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "state": row.state,
            "branch_length": row.branch_length,
            "mean_probability": row.mean_probability,
            "lower_95_probability": row.lower_95_probability,
            "upper_95_probability": row.upper_95_probability,
            "minimum_probability": row.minimum_probability,
            "maximum_probability": row.maximum_probability,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_state_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "branch_length",
            "mean_probability",
            "lower_95_probability",
            "upper_95_probability",
            "minimum_probability",
            "maximum_probability",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_density_branch_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one per-branch focal-state density summary table."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "branch_length": row.branch_length,
            "focal_state": row.focal_state,
            "baseline_state": row.baseline_state or "",
            "mean_posterior_probability": row.mean_posterior_probability,
            "minimum_posterior_probability": row.minimum_posterior_probability,
            "maximum_posterior_probability": row.maximum_posterior_probability,
            "uncertainty": row.uncertainty,
            "slice_count": row.slice_count,
        }
        for row in report.branch_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "branch_length",
            "focal_state",
            "baseline_state",
            "mean_posterior_probability",
            "minimum_posterior_probability",
            "maximum_posterior_probability",
            "uncertainty",
            "slice_count",
        ],
        rows=rows,
    )


def write_stochastic_map_density_slice_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one flat branch-slice density table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "branch_length": row.branch_length,
            "slice_index": row.slice_index,
            "start_depth": row.start_depth,
            "end_depth": row.end_depth,
            "start_time_fraction": row.start_time_fraction,
            "end_time_fraction": row.end_time_fraction,
            "posterior_probability": row.posterior_probability,
            "posterior_uncertainty": row.posterior_uncertainty,
        }
        for row in report.density_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "branch_length",
            "slice_index",
            "start_depth",
            "end_depth",
            "start_time_fraction",
            "end_time_fraction",
            "posterior_probability",
            "posterior_uncertainty",
        ],
        rows=rows,
    )


def write_stochastic_map_segment_table(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Export one flat branch-state segment table for a stochastic-map collection."""
    rows = [
        {
            "replicate_index": replicate.replicate_index,
            "branch_index": segment.branch_index,
            "parent_node": segment.parent_node,
            "child_node": segment.child_node,
            "state": segment.state,
            "start_time_fraction": segment.start_time_fraction,
            "end_time_fraction": segment.end_time_fraction,
            "duration": segment.duration,
        }
        for replicate in report.maps
        for history in replicate.branch_histories
        for segment in history.segments
    ]
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "start_time_fraction",
            "end_time_fraction",
            "duration",
        ],
        rows=rows,
    )


def write_stochastic_map_event_table(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Export one flat transition-event table for a stochastic-map collection."""
    rows = [
        {
            "replicate_index": replicate.replicate_index,
            "branch_index": history.branch_index,
            "parent_node": history.parent_node,
            "child_node": history.child_node,
            "event_index": event_index,
            "source_state": event.source_state,
            "target_state": event.target_state,
            "branch_length": history.branch_length,
            "event_time_fraction": event.event_time_fraction,
            "event_time": float(
                format(history.branch_length * event.event_time_fraction, ".15g")
            ),
        }
        for replicate in report.maps
        for history in replicate.branch_histories
        for event_index, event in enumerate(history.events)
    ]
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "branch_index",
            "parent_node",
            "child_node",
            "event_index",
            "source_state",
            "target_state",
            "branch_length",
            "event_time_fraction",
            "event_time",
        ],
        rows=rows,
    )


def write_stochastic_map_collection(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Write one stochastic-map collection as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def load_stochastic_map_collection(path: Path) -> StochasticMapCollectionReport:
    """Load one stochastic-map collection from JSON."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    maps = [
        StochasticMapReplicate(
            replicate_index=replicate["replicate_index"],
            root_state=replicate["root_state"],
            total_transition_count=replicate["total_transition_count"],
            transition_counts=replicate["transition_counts"],
            branch_histories=[
                StochasticMapBranchHistory(
                    branch_index=history["branch_index"],
                    parent_node=history["parent_node"],
                    child_node=history["child_node"],
                    branch_length=history["branch_length"],
                    start_state=history["start_state"],
                    end_state=history["end_state"],
                    event_count=history["event_count"],
                    events=[
                        StochasticMapTransitionEvent(
                            branch_index=event["branch_index"],
                            parent_node=event["parent_node"],
                            child_node=event["child_node"],
                            source_state=event["source_state"],
                            target_state=event["target_state"],
                            event_time_fraction=event["event_time_fraction"],
                        )
                        for event in history["events"]
                    ],
                    segments=[
                        StochasticMapStateSegment(
                            branch_index=segment["branch_index"],
                            parent_node=segment["parent_node"],
                            child_node=segment["child_node"],
                            state=segment["state"],
                            start_time_fraction=segment["start_time_fraction"],
                            end_time_fraction=segment["end_time_fraction"],
                            duration=segment["duration"],
                        )
                        for segment in history.get("segments", [])
                    ],
                )
                for history in replicate["branch_histories"]
            ],
            state_time_totals=replicate.get("state_time_totals", {}),
        )
        for replicate in payload["maps"]
    ]
    summary = StochasticMapSummaryReport(
        replicate_count=payload["summary"]["replicate_count"],
        mean_total_transition_count=payload["summary"]["mean_total_transition_count"],
        lower_95_total_transition_count=payload["summary"][
            "lower_95_total_transition_count"
        ],
        upper_95_total_transition_count=payload["summary"][
            "upper_95_total_transition_count"
        ],
        rows=[
            StochasticMapSummaryRow(
                transition=row["transition"],
                mean_count=row["mean_count"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_count=row["minimum_count"],
                maximum_count=row["maximum_count"],
                presence_fraction=row["presence_fraction"],
            )
            for row in payload["summary"]["rows"]
        ],
        state_time_rows=[
            StochasticMapStateTimeRow(
                state=row["state"],
                mean_time=row["mean_time"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_time=row["minimum_time"],
                maximum_time=row["maximum_time"],
            )
            for row in payload["summary"].get("state_time_rows", [])
        ],
        branch_occupancy_rows=[
            StochasticMapBranchOccupancyRow(
                branch_index=row["branch_index"],
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                state=row["state"],
                branch_length=row["branch_length"],
                mean_time=row["mean_time"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_time=row["minimum_time"],
                maximum_time=row["maximum_time"],
                mean_fraction=row.get("mean_fraction", 0.0),
                presence_fraction=row.get("presence_fraction", 1.0),
            )
            for row in payload["summary"].get("branch_occupancy_rows", [])
        ],
        simulation_failure_count=payload["summary"].get("simulation_failure_count", 0),
        warnings=payload["summary"]["warnings"],
    )
    return StochasticMapCollectionReport(
        tree_path=Path(payload["tree_path"]),
        traits_path=Path(payload["traits_path"]),
        taxon_column=payload["taxon_column"],
        trait=payload["trait"],
        model=payload["model"],
        state_ordering=payload["state_ordering"],
        ordered_states=payload["ordered_states"],
        replicates=payload["replicates"],
        seed=payload["seed"],
        conditioned_on_node_estimates=payload.get(
            "conditioned_on_node_estimates", False
        ),
        fit_audit=StochasticMapModelFitAudit(
            state_order=payload.get("fit_audit", {}).get("state_order", []),
            allowed_transitions=payload.get("fit_audit", {}).get(
                "allowed_transitions", []
            ),
            parameter_count=payload.get("fit_audit", {}).get("parameter_count", 0),
            log_likelihood=payload.get("fit_audit", {}).get("log_likelihood", 0.0),
            aic=payload.get("fit_audit", {}).get("aic", 0.0),
            aicc=payload.get("fit_audit", {}).get("aicc", 0.0),
            overparameterized=payload.get("fit_audit", {}).get(
                "overparameterized", False
            ),
            optimizer_converged=payload.get("fit_audit", {}).get(
                "optimizer_converged", True
            ),
            optimizer_iteration_count=payload.get("fit_audit", {}).get(
                "optimizer_iteration_count", 0
            ),
            optimizer_function_evaluation_count=payload.get("fit_audit", {}).get(
                "optimizer_function_evaluation_count", 0
            ),
            optimizer_hit_lower_parameter_bound=payload.get("fit_audit", {}).get(
                "optimizer_hit_lower_parameter_bound", False
            ),
            optimizer_hit_upper_parameter_bound=payload.get("fit_audit", {}).get(
                "optimizer_hit_upper_parameter_bound", False
            ),
            baseline_model=payload.get("fit_audit", {}).get("baseline_model"),
            baseline_aic=payload.get("fit_audit", {}).get("baseline_aic"),
            baseline_delta_aic=payload.get("fit_audit", {}).get("baseline_delta_aic"),
            preferred_model_by_aic=payload.get("fit_audit", {}).get(
                "preferred_model_by_aic"
            ),
            warnings=payload.get("fit_audit", {}).get("warnings", []),
        ),
        warnings=payload.get("warnings", []),
        maps=maps,
        failures=[
            StochasticMapSimulationFailure(
                replicate_index=row["replicate_index"],
                branch_index=row["branch_index"],
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                source_state=row["source_state"],
                target_state=row["target_state"],
                branch_length=row["branch_length"],
                attempt_count=row["attempt_count"],
                reason=row["reason"],
            )
            for row in payload.get("failures", [])
        ],
        summary=summary,
    )
