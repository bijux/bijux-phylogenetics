from __future__ import annotations

from bijux_phylogenetics.comparative.reporting import ComparativeMethodReport

from .contracts import (
    ComparativeAnalysisSummaryRow,
    ComparativeAuditTableRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)


def summarize_comparative_analysis(
    report: ComparativeMethodReport,
) -> ComparativeAnalysisSummaryRow:
    model_rows = report.snapshot.model_comparison.rows
    selected_row = next(row for row in model_rows if row.selected)
    runner_up_aicc = min(
        row.aicc for row in model_rows if row.model != selected_row.model
    )
    return ComparativeAnalysisSummaryRow(
        response=report.snapshot.response,
        formula=report.snapshot.formula.formula,
        predictor_count=len(report.snapshot.formula.predictors),
        analysis_taxa=report.snapshot.pgls_model.taxon_count,
        excluded_taxa=len(report.snapshot.pgls_inputs.formula_audit.excluded_taxa),
        selected_model=report.snapshot.model_comparison.better_model,
        pgls_lambda=report.snapshot.pgls_model.lambda_value,
        pgls_log_likelihood=report.snapshot.pgls_model.log_likelihood,
        pgls_r_squared=report.snapshot.pgls_model.r_squared,
        phylogenetic_signal_k=report.snapshot.signal_k.k,
        phylogenetic_signal_lambda=report.snapshot.signal_lambda.lambda_value,
        independent_contrast_count=len(report.snapshot.contrasts.contrasts),
        better_model_aicc_delta=runner_up_aicc - selected_row.aicc,
    )


def summarize_comparative_coefficients(
    report: ComparativeMethodReport,
) -> list[ComparativeCoefficientTableRow]:
    rows: list[ComparativeCoefficientTableRow] = []
    for coefficient in report.snapshot.pgls_model.coefficients:
        rows.append(
            ComparativeCoefficientTableRow(
                term=coefficient.name,
                estimate=coefficient.estimate,
                standard_error=coefficient.standard_error,
                test_statistic=coefficient.test_statistic,
                p_value=coefficient.p_value,
                lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
                upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
                degrees_of_freedom=coefficient.degrees_of_freedom,
                inference_distribution=coefficient.inference_distribution,
                significant=coefficient.p_value <= 0.05,
            )
        )
    return rows


def summarize_comparative_residuals(
    report: ComparativeMethodReport,
) -> list[ComparativeResidualTableRow]:
    rows: list[ComparativeResidualTableRow] = []
    for surface in report.snapshot.maturity.residual_diagnostics:
        rows.append(
            ComparativeResidualTableRow(
                analysis=surface.analysis,
                residual_variance=surface.residual_variance,
                max_abs_standardized_residual=surface.max_abs_standardized_residual,
                phylogenetic_residual_lambda=surface.phylogenetic_residual_lambda,
                max_leverage=surface.max_leverage,
                outlier_taxa=tuple(surface.outlier_taxa),
                high_leverage_taxa=tuple(surface.high_leverage_taxa),
                warnings=tuple(surface.warnings),
            )
        )
    return rows


def summarize_comparative_signal(
    report: ComparativeMethodReport,
) -> ComparativeSignalTableRow:
    return ComparativeSignalTableRow(
        trait=report.snapshot.signal_k.trait,
        taxon_count=report.snapshot.signal_k.taxon_count,
        blombergs_k=report.snapshot.signal_k.k,
        pagels_lambda=report.snapshot.signal_lambda.lambda_value,
        lambda_log_likelihood=report.snapshot.signal_lambda.log_likelihood,
        lambda_null_log_likelihood=report.snapshot.signal_lambda.null_log_likelihood,
        lambda_brownian_log_likelihood=(
            report.snapshot.signal_lambda.brownian_log_likelihood
        ),
        independent_contrast_count=len(report.snapshot.contrasts.contrasts),
        independent_contrast_root_estimate=report.snapshot.contrasts.root_estimate,
    )


def summarize_comparative_interpretation(
    report: ComparativeMethodReport,
) -> list[ComparativeInterpretationRow]:
    rows: list[ComparativeInterpretationRow] = []
    summary = summarize_comparative_analysis(report)
    rows.append(
        ComparativeInterpretationRow(
            topic="formula",
            claim=f"comparative regression was fit as `{summary.formula}`",
            evidence=(
                f"analysis_taxa={summary.analysis_taxa}; predictors={summary.predictor_count}"
            ),
            caution=(
                "formula interpretation is conditional on the supplied tree, taxon overlap, and encoded predictor structure"
            ),
        )
    )
    signal_level = (
        "high"
        if report.snapshot.signal_lambda.lambda_value >= 0.8
        else "moderate"
        if report.snapshot.signal_lambda.lambda_value >= 0.3
        else "low"
    )
    rows.append(
        ComparativeInterpretationRow(
            topic="phylogenetic-signal",
            claim=f"the response trait retains {signal_level} phylogenetic structure",
            evidence=(
                f"blombergs_k={report.snapshot.signal_k.k:.6f}; "
                f"pagels_lambda={report.snapshot.signal_lambda.lambda_value:.6f}"
            ),
            caution=(
                "signal metrics describe covariance with the supplied phylogeny and are not by themselves a process explanation"
            ),
        )
    )
    rows.append(
        ComparativeInterpretationRow(
            topic="model-comparison",
            claim=(
                f"{report.snapshot.model_comparison.better_model} is the preferred continuous trait model on AICc"
            ),
            evidence=(
                f"aicc_delta={summary.better_model_aicc_delta:.6f}; "
                f"selected_from={len(report.snapshot.model_comparison.rows)} models"
            ),
            caution=(
                "small AICc differences should not be overstated as decisive biological process separation"
            ),
        )
    )
    for row in summarize_comparative_coefficients(report):
        if row.term == "intercept":
            continue
        direction = "positive" if row.estimate > 0.0 else "negative"
        strength = (
            "nominally supported" if row.significant else "not nominally supported"
        )
        rows.append(
            ComparativeInterpretationRow(
                topic="coefficient",
                claim=f"`{row.term}` shows a {direction} association and is {strength}",
                evidence=(
                    f"estimate={row.estimate:.6f}; p_value={row.p_value:.6f}; "
                    f"95%_ci=[{row.lower_95_confidence_interval:.6f}, {row.upper_95_confidence_interval:.6f}]"
                ),
                caution=(
                    "coefficient interpretation remains associational and depends on lambda choice, sample size, and residual diagnostics"
                ),
            )
        )
    for residual in summarize_comparative_residuals(report):
        if residual.warnings:
            rows.append(
                ComparativeInterpretationRow(
                    topic="diagnostics",
                    claim=f"{residual.analysis} residual diagnostics flagged review concerns",
                    evidence="; ".join(residual.warnings),
                    caution=(
                        "diagnostic warnings weaken confidence in coefficient interpretation and should be reviewed before publication use"
                    ),
                )
            )
    return rows


def summarize_comparative_audit(
    report: ComparativeMethodReport,
) -> list[ComparativeAuditTableRow]:
    return [
        ComparativeAuditTableRow(
            analysis=row.analysis,
            taxa_used=tuple(row.taxa_used),
            traits_used=tuple(row.traits_used),
            excluded_taxa=tuple(row.excluded_taxa),
            assumptions=tuple(row.assumptions),
            warnings=tuple(row.warnings),
        )
        for row in report.snapshot.audit_rows
    ]
