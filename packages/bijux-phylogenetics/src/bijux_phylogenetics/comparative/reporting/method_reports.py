from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import PGLSTaxonExclusion
from bijux_phylogenetics.comparative.reporting.contracts import (
    ComparativeMethodReport,
    ComparativeMethodsSummaryTextResult,
)
from bijux_phylogenetics.comparative.reporting.influence import (
    build_trait_influence_report,
)
from bijux_phylogenetics.comparative.reporting.snapshot import (
    build_comparative_model_snapshot,
)
from bijux_phylogenetics.render.html import write_html_report


def build_comparative_method_report(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeMethodReport:
    """Build an integrated comparative-method report for one response trait."""
    snapshot = build_comparative_model_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    influence = build_trait_influence_report(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    return ComparativeMethodReport(snapshot=snapshot, influence=influence)


def write_comparative_method_report(
    path: Path, report: ComparativeMethodReport
) -> Path:
    """Render the comparative-method report to a standalone HTML artifact."""
    methods_summary = build_comparative_methods_summary_text(report)
    sections = [
        ("methods-summary-text", methods_summary),
        ("readiness", str(asdict(report.snapshot.readiness))),
        ("summary", str(asdict(report.snapshot.summary))),
        ("signal-k", str(asdict(report.snapshot.signal_k))),
        ("signal-lambda", str(asdict(report.snapshot.signal_lambda))),
        ("contrasts", str(asdict(report.snapshot.contrasts))),
        ("brownian", str(asdict(report.snapshot.brownian))),
        ("ou", str(asdict(report.snapshot.ou))),
        ("model-comparison", str(asdict(report.snapshot.model_comparison))),
        ("pgls-inputs", str(asdict(report.snapshot.pgls_inputs))),
        ("formula-audit", str(asdict(report.snapshot.pgls_inputs.formula_audit))),
        ("pgls-model", str(asdict(report.snapshot.pgls_model))),
        ("maturity", str(asdict(report.snapshot.maturity))),
        ("audit-table", str([asdict(row) for row in report.snapshot.audit_rows])),
        ("limitations", str(report.snapshot.limitations)),
        ("influence", str(asdict(report.influence))),
    ]
    return write_html_report(
        title="Bijux Comparative Method Report",
        sections=sections,
        out_path=path,
        embedded_json=json.loads(json.dumps(asdict(report), default=str)),
    )


def write_comparative_methods_summary_text(
    path: Path, report: ComparativeMethodReport
) -> ComparativeMethodsSummaryTextResult:
    """Write reviewer-facing Markdown methods text for one comparative analysis."""
    text = build_comparative_methods_summary_text(report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return ComparativeMethodsSummaryTextResult(
        output_path=path,
        title="Comparative Analysis Methods Summary",
        selected_model=report.snapshot.model_comparison.better_model,
        predictor_count=len(report.snapshot.formula.predictors),
        analysis_taxa=report.snapshot.pgls_model.taxon_count,
        excluded_taxa=len(report.snapshot.pgls_inputs.formula_audit.excluded_taxa),
        warning_count=len(_comparative_methods_summary_warnings(report)),
        text=text,
        report=report,
    )


def build_comparative_methods_summary_text(report: ComparativeMethodReport) -> str:
    """Build reviewer-facing Markdown methods text for one comparative analysis."""
    snapshot = report.snapshot
    summary = snapshot.summary
    readiness = snapshot.readiness
    formula_audit = snapshot.pgls_inputs.formula_audit
    selected_row = next(row for row in snapshot.model_comparison.rows if row.selected)
    runner_up_aicc = min(
        row.aicc
        for row in snapshot.model_comparison.rows
        if row.model != selected_row.model
    )
    categorical_predictors = [
        predictor.name
        for predictor in snapshot.pgls_inputs.predictors
        if predictor.kind != "numeric"
    ]
    transformed_terms = list(formula_audit.transformed_terms)
    warnings = _comparative_methods_summary_warnings(report)
    return (
        "# Comparative Analysis Methods Summary\n\n"
        f"This comparative analysis evaluated response trait `{snapshot.response}` on tree "
        f"`{snapshot.tree_path.name}` and trait table `{snapshot.traits_path.name}` using comparative formula "
        f"`{snapshot.formula.formula}`. Trait-process comparison reviewed Brownian and Ornstein-Uhlenbeck models, "
        f"preferred `{snapshot.model_comparison.better_model}` on AICc, and then fit a PGLS regression with Pagel lambda "
        f"mode `{snapshot.pgls_model.lambda_fit.mode}` resolved to `{format(snapshot.pgls_model.lambda_value, '.15g')}`.\n\n"
        "## Tree And Trait Pruning\n\n"
        f"- tree taxon count: `{readiness.tree_taxa}`\n"
        f"- rooted tree required: `{'yes' if readiness.rooted else 'no'}`\n"
        f"- complete branch lengths required: `{'yes' if readiness.complete_branch_lengths else 'no'}`\n"
        f"- analysis taxa retained after overlap and numeric pruning: `{summary.taxon_count}`\n"
        f"- overlapping taxa missing from the trait table: `{len(readiness.missing_from_traits)}`\n"
        f"- overlapping taxa pruned for missing response or predictor values: `{len(readiness.pruned_missing_value_taxa)}`\n"
        f"- overlapping taxa pruned for non-numeric values: `{len(readiness.pruned_non_numeric_taxa)}`\n"
        f"- formula-level excluded taxa: `{len(formula_audit.excluded_taxa)}`\n"
        f"- excluded taxa detail: {_comparative_exclusion_text(formula_audit.excluded_taxa)}\n"
        f"- retained response range: `{format(summary.minimum, '.15g')}` to `{format(summary.maximum, '.15g')}`\n\n"
        "## Predictors And Encoding\n\n"
        f"- response term: `{snapshot.formula.response}`\n"
        f"- predictor terms: {_comparative_bullet_list(snapshot.formula.predictors)}\n"
        f"- encoded design columns: {_comparative_bullet_list(formula_audit.encoded_columns)}\n"
        f"- residual degrees of freedom: `{formula_audit.residual_degrees_of_freedom}`\n"
        f"- minimum required taxa under the encoded design: `{formula_audit.minimum_required_taxa}`\n"
        f"- categorical predictors: {_comparative_bullet_list(categorical_predictors)}\n"
        f"- transformed predictor terms: {_comparative_bullet_list(transformed_terms)}\n\n"
        "## Comparative Model And Signal\n\n"
        f"- selected comparative process model: `{snapshot.model_comparison.better_model}`\n"
        f"- selected model AICc: `{format(selected_row.aicc, '.15g')}`\n"
        f"- runner-up AICc delta: `{format(runner_up_aicc - selected_row.aicc, '.15g')}`\n"
        f"- Blomberg's K: `{format(snapshot.signal_k.k, '.15g')}`\n"
        f"- Pagel's lambda signal estimate: `{format(snapshot.signal_lambda.lambda_value, '.15g')}`\n"
        f"- PGLS log likelihood: `{format(snapshot.pgls_model.log_likelihood, '.15g')}`\n"
        f"- PGLS R-squared: `{format(snapshot.pgls_model.r_squared, '.15g')}`\n\n"
        "## Assumptions And Interpretation Boundaries\n\n"
        + "\n".join(
            f"- {item}" for item in _comparative_methods_summary_assumptions(report)
        )
        + "\n\n## Diagnostics And Sensitivity\n\n"
        f"- maturity warnings: `{len(snapshot.maturity.warnings)}`\n"
        f"- influential taxa from leave-one-taxon-out sensitivity: {_comparative_bullet_list(snapshot.sensitivity.most_influential_taxa)}\n"
        f"- top high-influence taxa from combined leverage and sensitivity review: {_comparative_bullet_list(report.influence.top_taxa)}\n"
        f"- predictor terms with the strongest fitted influence: {_comparative_bullet_list(report.influence.top_predictor_terms)}\n"
        f"- diagnostic and reviewer warnings: {_comparative_bullet_list(warnings)}\n"
    )


def _comparative_exclusion_text(exclusions: list[PGLSTaxonExclusion]) -> str:
    if not exclusions:
        return "none"
    return "; ".join(
        f"`{row.taxon}` ({row.reason}: {row.details})" for row in exclusions
    )


def _comparative_bullet_list(items: list[str]) -> str:
    if not items:
        return "none"
    return ", ".join(f"`{item}`" for item in items)


def _comparative_methods_summary_assumptions(
    report: ComparativeMethodReport,
) -> list[str]:
    assumptions: list[str] = []
    for audit_row in report.snapshot.audit_rows:
        for assumption in audit_row.assumptions:
            if assumption not in assumptions:
                assumptions.append(assumption)
    for limitation in report.snapshot.limitations:
        if limitation not in assumptions:
            assumptions.append(limitation)
    return assumptions


def _comparative_methods_summary_warnings(
    report: ComparativeMethodReport,
) -> list[str]:
    warnings: list[str] = []
    for warning in report.snapshot.maturity.warnings:
        if warning not in warnings:
            warnings.append(warning)
    for audit_row in report.snapshot.audit_rows:
        for warning in audit_row.warnings:
            if warning not in warnings:
                warnings.append(warning)
    for warning in report.influence.warnings:
        if warning not in warnings:
            warnings.append(warning)
    return warnings
