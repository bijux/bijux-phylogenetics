from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.ancestral import (
    summarize_continuous_ancestral_report,
    summarize_discrete_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.comparative import (
    write_brownian_trait_evolution_exclusion_table,
    write_brownian_trait_evolution_summary_table,
    write_clade_trait_clade_table,
    write_clade_trait_exclusion_table,
    write_clade_trait_summary_table,
    write_ou_trait_evolution_exclusion_table,
    write_ou_trait_evolution_summary_table,
    write_pgls_lambda_profile_table,
    write_phylogenetic_signal_permutation_table,
    write_phylogenetic_signal_summary_table,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .models import (
    CentralEuropeanSeashoreFloraWorkflowBundle,
    CentralEuropeanSeashoreFloraWorkflowReport,
)


def write_central_european_seashore_flora_workflow_bundle(
    output_root: Path,
    report: CentralEuropeanSeashoreFloraWorkflowReport,
) -> CentralEuropeanSeashoreFloraWorkflowBundle:
    """Write the governed comparative workflow outputs for the packaged plant dataset."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    pgls_lambda_profile_path = write_pgls_lambda_profile_table(
        output_root / "pgls-lambda-profile.tsv",
        report.pgls.lambda_fit,
    )
    brownian_summary_path = write_brownian_trait_evolution_summary_table(
        output_root / "brownian-summary.tsv",
        report.brownian,
    )
    brownian_exclusion_path = write_brownian_trait_evolution_exclusion_table(
        output_root / "brownian-excluded.tsv",
        report.brownian,
    )
    ou_summary_path = write_ou_trait_evolution_summary_table(
        output_root / "ou-summary.tsv",
        report.ou,
    )
    ou_exclusion_path = write_ou_trait_evolution_exclusion_table(
        output_root / "ou-excluded.tsv",
        report.ou,
    )
    signal_summary_path = write_phylogenetic_signal_summary_table(
        output_root / "signal-summary.tsv",
        report.signal,
    )
    signal_permutations_path = write_phylogenetic_signal_permutation_table(
        output_root / "signal-permutations.tsv",
        report.signal,
    )
    continuous_ancestral_summary_path = write_continuous_ancestral_summary_table(
        output_root / "continuous-ancestral-summary.tsv",
        report.continuous_ancestral,
    )
    continuous_ancestral_uncertainty_path = (
        write_continuous_ancestral_uncertainty_table(
            output_root / "continuous-ancestral-uncertainty.tsv",
            report.continuous_ancestral,
        )
    )
    continuous_ancestral_exclusion_path = write_continuous_ancestral_exclusion_table(
        output_root / "continuous-ancestral-excluded.tsv",
        report.continuous_ancestral,
    )
    discrete_ancestral_summary_path = write_discrete_ancestral_summary_table(
        output_root / "discrete-ancestral-summary.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_probability_path = write_discrete_ancestral_probability_table(
        output_root / "discrete-ancestral-probabilities.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_exclusion_path = write_discrete_ancestral_exclusion_table(
        output_root / "discrete-ancestral-excluded.tsv",
        report.discrete_ancestral,
    )
    clade_summary_path = write_clade_trait_summary_table(
        output_root / "clade-trait-summary.tsv",
        report.clade_traits,
    )
    clade_rows_path = write_clade_trait_clade_table(
        output_root / "clade-trait-clades.tsv",
        report.clade_traits,
    )
    clade_exclusion_path = write_clade_trait_exclusion_table(
        output_root / "clade-trait-excluded.tsv",
        report.clade_traits,
    )
    return CentralEuropeanSeashoreFloraWorkflowBundle(
        output_root=output_root,
        summary_path=summary_path,
        pgls_lambda_profile_path=pgls_lambda_profile_path,
        brownian_summary_path=brownian_summary_path,
        brownian_exclusion_path=brownian_exclusion_path,
        ou_summary_path=ou_summary_path,
        ou_exclusion_path=ou_exclusion_path,
        signal_summary_path=signal_summary_path,
        signal_permutations_path=signal_permutations_path,
        continuous_ancestral_summary_path=continuous_ancestral_summary_path,
        continuous_ancestral_uncertainty_path=continuous_ancestral_uncertainty_path,
        continuous_ancestral_exclusion_path=continuous_ancestral_exclusion_path,
        discrete_ancestral_summary_path=discrete_ancestral_summary_path,
        discrete_ancestral_probability_path=discrete_ancestral_probability_path,
        discrete_ancestral_exclusion_path=discrete_ancestral_exclusion_path,
        clade_summary_path=clade_summary_path,
        clade_rows_path=clade_rows_path,
        clade_exclusion_path=clade_exclusion_path,
    )


def _write_workflow_summary_table(
    path: Path,
    report: CentralEuropeanSeashoreFloraWorkflowReport,
) -> Path:
    slope = next(
        coefficient
        for coefficient in report.pgls.coefficients
        if coefficient.name == report.dataset.workflow_pgls_predictor
    )
    continuous_summary = summarize_continuous_ancestral_report(
        report.continuous_ancestral
    )
    discrete_summary = summarize_discrete_ancestral_report(report.discrete_ancestral)
    continuous_root = next(
        estimate
        for estimate in report.continuous_ancestral.estimates
        if not estimate.is_tip and estimate.node == continuous_summary.root_node
    )
    discrete_root = next(
        estimate
        for estimate in report.discrete_ancestral.estimates
        if not estimate.is_tip and estimate.node == discrete_summary.root_node
    )
    return write_taxon_rows(
        path,
        columns=[
            "dataset_id",
            "taxon_count",
            "taxon_column",
            "continuous_trait",
            "pgls_predictor",
            "pgls_lambda",
            "pgls_predictor_estimate",
            "pgls_predictor_p_value",
            "brownian_sigma_squared",
            "brownian_aicc",
            "ou_alpha",
            "ou_theta",
            "ou_sigma_squared",
            "ou_aicc",
            "signal_blombergs_k",
            "signal_pagels_lambda",
            "signal_permutation_p_value",
            "continuous_root_node",
            "continuous_root_estimate",
            "discrete_trait",
            "discrete_root_node",
            "discrete_root_state",
            "discrete_root_confidence",
            "clade_trait",
            "top_exceptional_clade",
            "top_exceptionality_score",
            "exceptional_clade_count",
        ],
        rows=[
            {
                "dataset_id": report.dataset.dataset_id,
                "taxon_count": str(report.dataset.taxon_count),
                "taxon_column": report.dataset.taxon_column,
                "continuous_trait": report.dataset.workflow_continuous_trait,
                "pgls_predictor": report.dataset.workflow_pgls_predictor,
                "pgls_lambda": format(report.pgls.lambda_value, ".15g"),
                "pgls_predictor_estimate": format(slope.estimate, ".15g"),
                "pgls_predictor_p_value": format(slope.p_value, ".15g"),
                "brownian_sigma_squared": format(report.brownian.sigma_squared, ".15g"),
                "brownian_aicc": format(report.brownian.aicc, ".15g"),
                "ou_alpha": format(report.ou.alpha, ".15g"),
                "ou_theta": format(report.ou.theta, ".15g"),
                "ou_sigma_squared": format(report.ou.sigma_squared, ".15g"),
                "ou_aicc": format(report.ou.aicc, ".15g"),
                "signal_blombergs_k": format(report.signal.blombergs_k.k, ".15g"),
                "signal_pagels_lambda": format(
                    report.signal.pagels_lambda.lambda_value,
                    ".15g",
                ),
                "signal_permutation_p_value": format(
                    report.signal.signal_test.p_value,
                    ".15g",
                ),
                "continuous_root_node": continuous_root.node,
                "continuous_root_estimate": format(continuous_root.estimate, ".15g"),
                "discrete_trait": report.dataset.workflow_discrete_trait,
                "discrete_root_node": discrete_root.node,
                "discrete_root_state": discrete_root.most_likely_state,
                "discrete_root_confidence": format(discrete_root.confidence, ".15g"),
                "clade_trait": report.dataset.workflow_clade_trait,
                "top_exceptional_clade": report.clade_traits.top_exceptional_clade
                or "",
                "top_exceptionality_score": format(
                    report.clade_traits.top_exceptionality_score or 0.0,
                    ".15g",
                ),
                "exceptional_clade_count": len(report.clade_traits.exceptional_clades),
            }
        ],
    )
