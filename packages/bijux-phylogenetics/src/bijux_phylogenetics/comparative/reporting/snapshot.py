from __future__ import annotations

from bijux_phylogenetics.comparative.assessment import (
    assess_comparative_method_maturity,
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous import (
    BrownianMotionFitReport,
    OUTraitEvolutionSummaryReport,
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    summarize_ou_trait_evolution,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSInputReport,
    PGLSResult,
    inspect_pgls_inputs,
    run_pgls,
)
from bijux_phylogenetics.comparative.reporting.contracts import (
    ComparativeAuditRow,
    ComparativeModelSnapshot,
)
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_independent_contrasts,
    estimate_pagels_lambda,
)


def build_comparative_model_snapshot(
    tree_path,
    traits_path,
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
    taxon_column: str | None,
    lambda_value: float | str,
) -> ComparativeModelSnapshot:
    pgls_inputs = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    summary = summarize_numeric_trait(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
    )
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        taxon_column=taxon_column,
    )
    signal_k = compute_blombergs_k(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    signal_lambda = estimate_pagels_lambda(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    contrasts = compute_phylogenetic_independent_contrasts(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    brownian = fit_brownian_motion_model(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    ou = summarize_ou_trait_evolution(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    model_comparison = compare_brownian_and_ou_models(
        tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column
    )
    pgls_model = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    sensitivity = run_comparative_sensitivity_analysis(
        tree_path,
        traits_path,
        trait=pgls_inputs.response,
        model=model_comparison.better_model,
        taxon_column=taxon_column,
    )
    maturity = assess_comparative_method_maturity(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    limitations = build_comparative_limitations(
        readiness, pgls_inputs, pgls_model, brownian, ou
    )
    return ComparativeModelSnapshot(
        tree_path=tree_path,
        traits_path=traits_path,
        response=pgls_inputs.response,
        formula=pgls_inputs.formula,
        readiness=readiness,
        summary=summary,
        signal_k=signal_k,
        signal_lambda=signal_lambda,
        contrasts=contrasts,
        brownian=brownian,
        ou=ou,
        model_comparison=model_comparison,
        pgls_inputs=pgls_inputs,
        pgls_model=pgls_model,
        sensitivity=sensitivity,
        maturity=maturity,
        audit_rows=build_comparative_audit_rows(
            readiness, pgls_inputs, pgls_model, brownian, ou
        ),
        limitations=limitations,
    )


def build_comparative_audit_rows(
    readiness: ComparativeReadinessReport,
    pgls_inputs: PGLSInputReport,
    pgls_model: PGLSResult,
    brownian: BrownianMotionFitReport,
    ou: OUTraitEvolutionSummaryReport,
) -> list[ComparativeAuditRow]:
    excluded_taxa = sorted(
        set(readiness.missing_from_traits)
        | set(readiness.pruned_missing_value_taxa)
        | set(readiness.pruned_non_numeric_taxa)
    )
    return [
        ComparativeAuditRow(
            analysis="brownian",
            taxa_used=brownian.taxa,
            traits_used=[brownian.trait],
            excluded_taxa=excluded_taxa,
            assumptions=brownian.assumptions,
            warnings=list(brownian.residual_diagnostics.warnings),
        ),
        ComparativeAuditRow(
            analysis="ou",
            taxa_used=ou.analyzed_taxa,
            traits_used=[ou.trait],
            excluded_taxa=excluded_taxa,
            assumptions=ou.assumptions,
            warnings=[
                *ou.residual_diagnostics.warnings,
                *[warning.message for warning in ou.identifiability_warnings],
            ],
        ),
        ComparativeAuditRow(
            analysis="pgls",
            taxa_used=pgls_model.taxa,
            traits_used=[pgls_model.response, *pgls_inputs.formula.predictors],
            excluded_taxa=excluded_taxa,
            assumptions=[
                "PGLS assumes the specified predictors explain trait variation on the supplied phylogeny",
                "PGLS relies on dummy encoding for categorical predictors and explicit interaction expansion",
            ],
            warnings=pgls_inputs.warnings
            + [f"residual_mean={pgls_model.diagnostics.residual_mean:.6f}"],
        ),
    ]


def build_comparative_limitations(
    readiness: ComparativeReadinessReport,
    pgls_inputs: PGLSInputReport,
    pgls_model: PGLSResult,
    brownian: BrownianMotionFitReport,
    ou: OUTraitEvolutionSummaryReport,
) -> list[str]:
    limitations = [
        "comparative conclusions are conditioned on a single observed trait table and one supplied phylogeny",
        "causal interpretation is not warranted from comparative association alone",
        "comparative associations should not be treated as intervention or mechanism evidence without external biological support",
        "single-tree comparative fits do not quantify the full effect of phylogenetic uncertainty unless explicit tree-sensitivity checks are reviewed",
    ]
    limitations.extend(readiness.warnings)
    limitations.extend(brownian.residual_diagnostics.warnings)
    limitations.extend(ou.residual_diagnostics.warnings)
    limitations.extend(warning.message for warning in ou.identifiability_warnings)
    if pgls_inputs.residual_degrees_of_freedom <= 1:
        limitations.append(
            "PGLS residual degrees of freedom are minimal, so coefficient uncertainty is fragile"
        )
        limitations.append(
            "do not treat coefficient signs or p-values as robust when the model is close to saturation"
        )
    if pgls_inputs.categorical_predictors:
        limitations.append(
            "categorical predictors are interpreted relative to explicit reference levels"
        )
        limitations.append(
            "do not compare categorical coefficients as if they were absolute trait differences outside their encoded reference-level context"
        )
    if pgls_model.lambda_value in {0.0, 1.0}:
        limitations.append(
            "estimated or supplied lambda lies on the boundary of the supported search interval"
        )
        limitations.append(
            "do not over-interpret the phylogenetic dependence strength when lambda is supported only at a boundary value"
        )
    if pgls_inputs.formula_audit.transformed_terms:
        limitations.append(
            "transformed predictor terms are interpreted on transformed scales and should not be read as raw-unit effect sizes"
        )
    if pgls_inputs.formula_audit.excluded_taxa:
        limitations.append(
            "excluded taxa can materially change comparative conclusions and should be reviewed before publication use"
        )
    if ou.identifiability_warnings:
        limitations.append(
            "do not interpret apparent OU preference as strong evidence of stabilizing selection when OU identifiability warnings are present"
        )
    return sorted(set(limitations))
