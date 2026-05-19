from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.simulation import (
    ContinuousTraitSimulationCollectionReport,
    CorrelatedContinuousTraitSimulationCollectionReport,
    simulate_brownian_trait_collection,
    simulate_correlated_brownian_trait_collection,
)

from ..registry import PhytoolsParityCase


def _continuous_trait_collection_parity_rows(
    report: ContinuousTraitSimulationCollectionReport,
) -> list[dict[str, object]]:
    return [
        {
            "row_kind": row.row_kind,
            "label": row.label,
            "mean_value": "" if row.mean_value is None else row.mean_value,
            "standard_deviation": (
                "" if row.standard_deviation is None else row.standard_deviation
            ),
            "minimum": "" if row.minimum is None else row.minimum,
            "median": "" if row.median is None else row.median,
            "maximum": "" if row.maximum is None else row.maximum,
            "covariance": "" if row.covariance is None else row.covariance,
            "correlation": "" if row.correlation is None else row.correlation,
        }
        for row in report.rows
    ]


def _correlated_continuous_trait_collection_parity_rows(
    report: CorrelatedContinuousTraitSimulationCollectionReport,
) -> list[dict[str, object]]:
    return [
        {
            "row_kind": row.row_kind,
            "label": row.label,
            "mean_value": "" if row.mean_value is None else row.mean_value,
            "standard_deviation": (
                "" if row.standard_deviation is None else row.standard_deviation
            ),
            "minimum": "" if row.minimum is None else row.minimum,
            "median": "" if row.median is None else row.median,
            "maximum": "" if row.maximum is None else row.maximum,
            "covariance": "" if row.covariance is None else row.covariance,
            "correlation": "" if row.correlation is None else row.correlation,
        }
        for row in report.rows
    ]


def build_continuous_case_payload(
    case: PhytoolsParityCase,
    *,
    tree_path: Path,
    traits_path: Path | None,
) -> tuple[dict[str, object], list[dict[str, object]] | None] | None:
    if case.operation == "simulate-continuous-brownian":
        report = simulate_brownian_trait_collection(
            tree_path,
            root_state=case.continuous_root_state or 0.0,
            sigma_squared=case.continuous_sigma_squared,
            replicates=case.continuous_replicate_count or 256,
            seed=case.continuous_seed or 1,
        )
        return (
            {
                "taxon_count": report.tip_count,
                "branch_count": report.branch_count,
                "requested_replicate_count": report.replicate_count,
                "successful_replicate_count": report.replicate_count,
                "seed": report.seed,
                "root_state": report.root_state,
                "sigma_squared": report.sigma_squared,
            },
            _continuous_trait_collection_parity_rows(report),
        )
    if case.operation == "simulate-continuous-correlated-brownian":
        report = simulate_correlated_brownian_trait_collection(
            tree_path,
            trait_names=list(case.continuous_trait_names or ()),
            evolutionary_covariance_matrix=[
                list(row) for row in (case.continuous_covariance_matrix or ())
            ],
            root_states=list(case.continuous_root_states or ()),
            replicates=case.continuous_replicate_count or 256,
            seed=case.continuous_seed or 1,
        )
        return (
            {
                "taxon_count": report.tip_count,
                "branch_count": report.branch_count,
                "trait_count": len(report.trait_names),
                "requested_replicate_count": report.replicate_count,
                "successful_replicate_count": report.replicate_count,
                "seed": report.seed,
            },
            _correlated_continuous_trait_collection_parity_rows(report),
        )
    if case.operation == "continuous-ancestral-fast-anc":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model="brownian",
            estimator="fast-anc",
        )
        rows = sorted(
            [
                {
                    "node": estimate.node,
                    "estimate": estimate.estimate,
                    "standard_error": estimate.standard_error,
                    "lower_95_interval": estimate.lower_95_interval,
                    "upper_95_interval": estimate.upper_95_interval,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
            key=lambda row: str(row["node"]),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "internal_node_count": len(rows),
                "excluded_taxon_count": len(report.dropped_missing_taxa)
                + len(report.dropped_non_numeric_taxa),
                "excluded_taxa": sorted(
                    report.dropped_missing_taxa + report.dropped_non_numeric_taxa
                ),
                "tree_is_ultrametric": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.tree_is_ultrametric
                ),
                "covariance_condition_number": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.covariance_condition_number
                ),
                "log_likelihood": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.log_likelihood
                ),
                "warning_count": len(report.warnings),
            },
            rows,
        )
    if case.operation == "continuous-ancestral-anc-ml":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model="brownian",
            estimator="anc-ml",
        )
        rows = sorted(
            [
                {
                    "node": estimate.node,
                    "estimate": estimate.estimate,
                    "standard_error": estimate.standard_error,
                    "lower_95_interval": estimate.lower_95_interval,
                    "upper_95_interval": estimate.upper_95_interval,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
            key=lambda row: str(row["node"]),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "internal_node_count": len(rows),
                "excluded_taxon_count": len(report.dropped_missing_taxa)
                + len(report.dropped_non_numeric_taxa),
                "excluded_taxa": sorted(
                    report.dropped_missing_taxa + report.dropped_non_numeric_taxa
                ),
                "tree_is_ultrametric": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.tree_is_ultrametric
                ),
                "sigma_squared": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.residual_sigma_squared
                ),
                "log_likelihood": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.log_likelihood
                ),
                "warning_count": len(report.warnings),
            },
            rows,
        )
    return None
