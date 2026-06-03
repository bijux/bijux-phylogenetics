from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_discrete_dataset
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    count_discrete_stochastic_map_transitions,
    simulate_discrete_stochastic_maps,
    summarize_discrete_stochastic_map_density,
    summarize_discrete_stochastic_maps,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.simulation import simulate_discrete_histories

from ..registry import PhytoolsParityCase


def _stochastic_map_parity_rows(
    summary,
    *,
    include_branch_occupancy: bool,
) -> list[dict[str, object]]:
    rows = [
        {
            "row_kind": "transition_count",
            "label": row.transition,
            "mean_value": row.mean_count,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "presence_fraction": row.presence_fraction,
        }
        for row in summary.rows
    ] + [
        {
            "row_kind": "state_time",
            "label": row.state,
            "mean_value": row.mean_time,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "presence_fraction": 1.0,
        }
        for row in summary.state_time_rows
    ]
    if include_branch_occupancy:
        rows.extend(
            {
                "row_kind": "branch_state_occupancy",
                "label": f"{row.parent_node}->{row.child_node}:{row.state}",
                "mean_value": row.mean_time,
                "lower_95_interval": row.lower_95_interval,
                "upper_95_interval": row.upper_95_interval,
                "presence_fraction": row.presence_fraction,
            }
            for row in summary.branch_occupancy_rows
        )
    return sorted(rows, key=lambda row: (str(row["row_kind"]), str(row["label"])))


def _stochastic_map_count_parity_rows(report) -> list[dict[str, object]]:
    return sorted(
        [
            {
                "row_kind": "transition_count",
                "label": row.transition,
                "mean_value": row.mean_count,
                "lower_95_interval": row.lower_95_interval,
                "upper_95_interval": row.upper_95_interval,
                "presence_fraction": row.presence_fraction,
            }
            for row in report.aggregate_rows
        ],
        key=lambda row: (str(row["row_kind"]), str(row["label"])),
    )


def _stochastic_map_density_parity_rows(report) -> list[dict[str, object]]:
    return sorted(
        [
            {
                "label": f"{row.parent_node}->{row.child_node}",
                "mean_posterior_probability": row.mean_posterior_probability,
                "minimum_posterior_probability": row.minimum_posterior_probability,
                "maximum_posterior_probability": row.maximum_posterior_probability,
                "uncertainty": row.uncertainty,
                "slice_count": row.slice_count,
            }
            for row in report.branch_rows
        ],
        key=lambda row: str(row["label"]),
    )


def _discrete_history_parity_rows(report) -> list[dict[str, object]]:
    return [
        {
            "row_kind": row.row_kind,
            "label": row.label,
            "mean_value": row.mean_value,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.rows
    ]


def build_discrete_case_payload(
    case: PhytoolsParityCase,
    *,
    tree_path: Path,
    traits_path: Path | None,
) -> tuple[dict[str, object], list[dict[str, object]] | None] | None:
    if case.operation == "discrete-fit-mk":
        report = fit_discrete_mk_model(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
        )
        rows = sorted(
            [
                {
                    "source_state": row.source_state,
                    "target_state": row.target_state,
                    "transition_allowed": row.transition_allowed,
                    "step_distance": row.step_distance,
                    "rate": row.rate,
                }
                for row in report.transition_rate_rows
            ],
            key=lambda row: (str(row["source_state"]), str(row["target_state"])),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "excluded_taxon_count": len(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "excluded_taxa": list(report.input_audit.pruned_missing_value_taxa),
                "model": report.model,
                "state_count": len(report.input_audit.observed_states),
                "parameter_count": report.parameter_count,
                "log_likelihood": report.log_likelihood,
                "aic": report.aic,
                "aicc": report.aicc,
                "overparameterized": report.overparameterized,
                "baseline_model": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.baseline_model
                ),
                "preferred_model_by_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.preferred_model_by_aic
                ),
            },
            rows,
        )
    if case.operation == "discrete-stochastic-map":
        dataset = load_discrete_dataset(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        report = simulate_discrete_stochastic_maps(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
            replicates=case.stochastic_map_replicate_count or 128,
            seed=case.stochastic_map_seed or 1,
        )
        return (
            {
                "taxon_count": len(dataset.taxa),
                "trait_name": report.trait,
                "excluded_taxon_count": len(dataset.dropped_missing_taxa),
                "excluded_taxa": list(dataset.dropped_missing_taxa),
                "model": report.model,
                "state_count": len(report.fit_audit.state_order),
                "parameter_count": report.fit_audit.parameter_count,
                "log_likelihood": report.fit_audit.log_likelihood,
                "aic": report.fit_audit.aic,
                "aicc": report.fit_audit.aicc,
                "overparameterized": report.fit_audit.overparameterized,
                "baseline_model": report.fit_audit.baseline_model,
                "preferred_model_by_aic": report.fit_audit.preferred_model_by_aic,
                "requested_replicate_count": report.replicates,
                "successful_replicate_count": report.summary.replicate_count,
                "simulation_failure_count": report.summary.simulation_failure_count,
                "conditioned_on_node_estimates": report.conditioned_on_node_estimates,
                "seed": report.seed,
                "mean_total_transition_count": report.summary.mean_total_transition_count,
                "lower_95_total_transition_count": report.summary.lower_95_total_transition_count,
                "upper_95_total_transition_count": report.summary.upper_95_total_transition_count,
            },
            _stochastic_map_parity_rows(
                report.summary,
                include_branch_occupancy=False,
            ),
        )
    if case.operation == "discrete-stochastic-map-description":
        dataset = load_discrete_dataset(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        collection = simulate_discrete_stochastic_maps(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
            replicates=case.stochastic_map_replicate_count or 128,
            seed=case.stochastic_map_seed or 1,
        )
        report = summarize_discrete_stochastic_maps(collection)
        return (
            {
                "taxon_count": len(dataset.taxa),
                "trait_name": collection.trait,
                "excluded_taxon_count": len(dataset.dropped_missing_taxa),
                "excluded_taxa": list(dataset.dropped_missing_taxa),
                "model": collection.model,
                "state_count": len(collection.fit_audit.state_order),
                "parameter_count": collection.fit_audit.parameter_count,
                "log_likelihood": collection.fit_audit.log_likelihood,
                "aic": collection.fit_audit.aic,
                "aicc": collection.fit_audit.aicc,
                "overparameterized": collection.fit_audit.overparameterized,
                "baseline_model": collection.fit_audit.baseline_model,
                "preferred_model_by_aic": collection.fit_audit.preferred_model_by_aic,
                "requested_replicate_count": collection.replicates,
                "successful_replicate_count": report.replicate_count,
                "simulation_failure_count": report.simulation_failure_count,
                "seed": collection.seed,
                "branch_count": len(collection.maps[0].branch_histories),
                "mean_total_transition_count": report.mean_total_transition_count,
                "lower_95_total_transition_count": report.lower_95_total_transition_count,
                "upper_95_total_transition_count": report.upper_95_total_transition_count,
            },
            _stochastic_map_parity_rows(
                report,
                include_branch_occupancy=True,
            ),
        )
    if case.operation == "discrete-stochastic-map-count":
        dataset = load_discrete_dataset(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        collection = simulate_discrete_stochastic_maps(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
            replicates=case.stochastic_map_replicate_count or 128,
            seed=case.stochastic_map_seed or 1,
        )
        report = count_discrete_stochastic_map_transitions(collection)
        return (
            {
                "taxon_count": len(dataset.taxa),
                "trait_name": collection.trait,
                "excluded_taxon_count": len(dataset.dropped_missing_taxa),
                "excluded_taxa": list(dataset.dropped_missing_taxa),
                "model": collection.model,
                "state_count": len(collection.fit_audit.state_order),
                "parameter_count": collection.fit_audit.parameter_count,
                "log_likelihood": collection.fit_audit.log_likelihood,
                "aic": collection.fit_audit.aic,
                "aicc": collection.fit_audit.aicc,
                "overparameterized": collection.fit_audit.overparameterized,
                "baseline_model": collection.fit_audit.baseline_model,
                "preferred_model_by_aic": collection.fit_audit.preferred_model_by_aic,
                "requested_replicate_count": collection.replicates,
                "successful_replicate_count": report.replicate_count,
                "simulation_failure_count": len(collection.failures),
                "seed": collection.seed,
                "mean_total_transition_count": report.mean_total_transition_count,
                "lower_95_total_transition_count": report.lower_95_total_transition_count,
                "upper_95_total_transition_count": report.upper_95_total_transition_count,
            },
            _stochastic_map_count_parity_rows(report),
        )
    if case.operation == "discrete-stochastic-map-density":
        if traits_path is None:
            raise ValueError(
                "discrete-stochastic-map-density requires one traits fixture"
            )
        dataset = load_discrete_dataset(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        collection = simulate_discrete_stochastic_maps(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
            replicates=case.stochastic_map_replicate_count or 128,
            seed=case.stochastic_map_seed or 1,
        )
        report = summarize_discrete_stochastic_map_density(
            collection,
            resolution=case.density_resolution or 100,
            focal_state=case.focal_state,
        )
        return (
            {
                "taxon_count": len(dataset.taxa),
                "trait_name": collection.trait,
                "excluded_taxon_count": len(dataset.dropped_missing_taxa),
                "excluded_taxa": list(dataset.dropped_missing_taxa),
                "model": collection.model,
                "state_count": len(report.state_order),
                "parameter_count": collection.fit_audit.parameter_count,
                "log_likelihood": collection.fit_audit.log_likelihood,
                "aic": collection.fit_audit.aic,
                "aicc": collection.fit_audit.aicc,
                "overparameterized": collection.fit_audit.overparameterized,
                "baseline_model": collection.fit_audit.baseline_model,
                "preferred_model_by_aic": collection.fit_audit.preferred_model_by_aic,
                "requested_replicate_count": collection.replicates,
                "successful_replicate_count": report.replicate_count,
                "simulation_failure_count": len(collection.failures),
                "seed": collection.seed,
                "branch_count": len(report.branch_rows),
                "focal_state": report.focal_state,
                "baseline_state": report.baseline_state,
                "resolution": report.resolution,
                "total_tree_depth": report.total_tree_depth,
            },
            _stochastic_map_density_parity_rows(report),
        )
    if case.operation == "simulate-discrete-history":
        report = simulate_discrete_histories(
            tree_path,
            states=list(case.simulation_states or ()),
            rate_rows=list(case.simulation_rate_rows or ()),
            root_state=case.simulation_root_state,
            root_state_probabilities=case.simulation_root_state_probabilities,
            replicates=case.simulation_replicate_count or 128,
            seed=case.simulation_seed or 1,
        )
        return (
            {
                "taxon_count": report.tip_count,
                "trait_name": case.trait_name,
                "branch_count": report.branch_count,
                "state_count": len(report.states),
                "requested_replicate_count": report.replicate_count,
                "successful_replicate_count": report.replicate_count,
                "fixed_root_state": report.fixed_root_state,
                "root_prior_probabilities": (
                    None
                    if report.fixed_root_state is not None
                    else report.root_state_probabilities
                ),
                "seed": report.seed,
                "mean_total_transition_count": report.mean_total_transition_count,
                "lower_95_total_transition_count": report.lower_95_total_transition_count,
                "upper_95_total_transition_count": report.upper_95_total_transition_count,
            },
            _discrete_history_parity_rows(report),
        )
    if case.operation == "discrete-ancestral-rerooting":
        if traits_path is None:
            raise ValueError("discrete-ancestral-rerooting requires one traits fixture")
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
            root_prior_mode=case.root_prior_mode,
        )
        rows = sorted(
            [
                {
                    "node": estimate.node,
                    "state": state,
                    "probability": probability,
                }
                for estimate in report.estimates
                if not estimate.is_tip
                for state, probability in estimate.state_probabilities.items()
            ],
            key=lambda row: (str(row["node"]), str(row["state"])),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "excluded_taxon_count": len(report.dropped_missing_taxa),
                "excluded_taxa": list(report.dropped_missing_taxa),
                "model": report.model,
                "state_count": len(report.observed_states),
                "internal_node_count": sum(
                    1 for estimate in report.estimates if not estimate.is_tip
                ),
                "root_prior_mode": report.root_prior_mode,
                "phytools_rerooting_method_comparable": (
                    report.rerooting_method_compatibility.comparable
                ),
            },
            rows,
        )
    return None
