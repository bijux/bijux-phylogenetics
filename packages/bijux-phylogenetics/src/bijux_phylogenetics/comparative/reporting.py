from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.common import ComparativeReadinessReport, NumericTraitSummary, summarize_numeric_trait, summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.models import (
    BrownianMotionFitReport,
    ComparativeModelComparisonReport,
    ComparativeSensitivityReport,
    OUTraitModelReport,
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.pgls import ComparativeFormulaSpecification, PGLSInputReport, PGLSResult, inspect_pgls_inputs, run_pgls
from bijux_phylogenetics.comparative.signal import BlombergKReport, IndependentContrastReport, PagelLambdaReport, compute_blombergs_k, compute_phylogenetic_independent_contrasts, estimate_pagels_lambda
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.html import write_html_report


@dataclass(slots=True)
class ComparativeAuditRow:
    """Auditable record for one comparative analysis surface."""

    analysis: str
    taxa_used: list[str]
    traits_used: list[str]
    excluded_taxa: list[str]
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativePredictorInfluenceRow:
    """One predictor term ranked by effect and test strength."""

    term: str
    estimate: float
    z_score: float
    p_value: float


@dataclass(slots=True)
class ComparativeTaxonInfluenceRow:
    """Combined taxon influence summary across PGLS and leave-one-out sensitivity."""

    taxon: str
    leverage: float
    standardized_residual: float
    sensitivity_delta_log_likelihood: float
    sensitivity_delta_primary_parameter: float
    influence_score: float


@dataclass(slots=True)
class ComparativeInfluenceReport:
    """Predictor and taxon drivers for one comparative analysis."""

    tree_path: Path
    traits_path: Path
    response: str
    selected_model: str
    predictor_rows: list[ComparativePredictorInfluenceRow]
    taxon_rows: list[ComparativeTaxonInfluenceRow]


@dataclass(slots=True)
class ComparativeModelSnapshot:
    """Shared comparative snapshot used by reports and comparisons."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    readiness: ComparativeReadinessReport
    summary: NumericTraitSummary
    signal_k: BlombergKReport
    signal_lambda: PagelLambdaReport
    contrasts: IndependentContrastReport
    brownian: BrownianMotionFitReport
    ou: OUTraitModelReport
    model_comparison: ComparativeModelComparisonReport
    pgls_inputs: PGLSInputReport
    pgls_model: PGLSResult
    sensitivity: ComparativeSensitivityReport
    audit_rows: list[ComparativeAuditRow]
    limitations: list[str]


@dataclass(slots=True)
class ComparativeMethodReport:
    """Integrated comparative-method report for one response trait."""

    snapshot: ComparativeModelSnapshot
    influence: ComparativeInfluenceReport


@dataclass(slots=True)
class ComparativeCoefficientDeltaRow:
    """Difference in one encoded coefficient between two comparative analyses."""

    term: str
    left_estimate: float | None
    right_estimate: float | None
    delta: float | None


@dataclass(slots=True)
class ComparativeTreeComparisonReport:
    """Compare comparative results across two alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    response: str
    predictors: list[str]
    left_selected_model: str
    right_selected_model: str
    delta_blombergs_k: float
    delta_pagels_lambda: float
    delta_brownian_rate: float
    delta_ou_alpha: float
    coefficient_deltas: list[ComparativeCoefficientDeltaRow]


@dataclass(slots=True)
class ComparativePruningComparisonReport:
    """Compare comparative results before and after explicit taxon pruning."""

    tree_path: Path
    response: str
    predictors: list[str]
    baseline_taxa: list[str]
    pruned_taxa: list[str]
    dropped_taxa: list[str]
    delta_blombergs_k: float
    delta_pagels_lambda: float
    coefficient_deltas: list[ComparativeCoefficientDeltaRow]
    baseline_selected_model: str
    pruned_selected_model: str


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
    snapshot = _build_snapshot(
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


def write_comparative_method_report(path: Path, report: ComparativeMethodReport) -> Path:
    """Render the comparative-method report to a standalone HTML artifact."""
    sections = [
        ("readiness", str(asdict(report.snapshot.readiness))),
        ("summary", str(asdict(report.snapshot.summary))),
        ("signal-k", str(asdict(report.snapshot.signal_k))),
        ("signal-lambda", str(asdict(report.snapshot.signal_lambda))),
        ("contrasts", str(asdict(report.snapshot.contrasts))),
        ("brownian", str(asdict(report.snapshot.brownian))),
        ("ou", str(asdict(report.snapshot.ou))),
        ("model-comparison", str(asdict(report.snapshot.model_comparison))),
        ("pgls-inputs", str(asdict(report.snapshot.pgls_inputs))),
        ("pgls-model", str(asdict(report.snapshot.pgls_model))),
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


def build_trait_influence_report(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeInfluenceReport:
    """Identify predictor terms and taxa driving one comparative result."""
    snapshot = _build_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    predictor_rows = sorted(
        [
            ComparativePredictorInfluenceRow(
                term=coefficient.name,
                estimate=coefficient.estimate,
                z_score=coefficient.z_score,
                p_value=coefficient.p_value,
            )
            for coefficient in snapshot.pgls_model.coefficients
            if coefficient.name != "intercept"
        ],
        key=lambda row: (abs(row.z_score), abs(row.estimate), row.term),
        reverse=True,
    )
    leverage_by_taxon = {
        row.taxon: row for row in snapshot.pgls_model.diagnostics.leverage_rows
    }
    sensitivity_by_taxon = {
        row.dropped_taxon: row for row in snapshot.sensitivity.rows
    }
    taxon_rows = [
        ComparativeTaxonInfluenceRow(
            taxon=taxon,
            leverage=leverage_by_taxon[taxon].leverage,
            standardized_residual=leverage_by_taxon[taxon].standardized_residual,
            sensitivity_delta_log_likelihood=sensitivity_by_taxon[taxon].delta_log_likelihood,
            sensitivity_delta_primary_parameter=sensitivity_by_taxon[taxon].delta_primary_parameter,
            influence_score=(
                abs(leverage_by_taxon[taxon].standardized_residual)
                + leverage_by_taxon[taxon].leverage
                + abs(sensitivity_by_taxon[taxon].delta_log_likelihood)
            ),
        )
        for taxon in snapshot.pgls_model.taxa
    ]
    taxon_rows.sort(key=lambda row: (row.influence_score, row.taxon), reverse=True)
    return ComparativeInfluenceReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response=snapshot.response,
        selected_model=snapshot.model_comparison.better_model,
        predictor_rows=predictor_rows,
        taxon_rows=taxon_rows,
    )


def compare_comparative_results_across_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeTreeComparisonReport:
    """Compare comparative fits across two alternative trees."""
    left = _build_snapshot(
        left_tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    right = _build_snapshot(
        right_tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    return ComparativeTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        response=left.response,
        predictors=left.formula.predictors,
        left_selected_model=left.model_comparison.better_model,
        right_selected_model=right.model_comparison.better_model,
        delta_blombergs_k=right.signal_k.k - left.signal_k.k,
        delta_pagels_lambda=right.signal_lambda.lambda_value - left.signal_lambda.lambda_value,
        delta_brownian_rate=right.brownian.rate - left.brownian.rate,
        delta_ou_alpha=right.ou.alpha - left.ou.alpha,
        coefficient_deltas=_coefficient_delta_rows(left.pgls_model, right.pgls_model),
    )


def compare_comparative_results_across_pruning(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    drop_taxa: list[str] | None = None,
    keep_taxa: list[str] | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativePruningComparisonReport:
    """Compare comparative fits before and after explicit tree/trait pruning."""
    if bool(drop_taxa) == bool(keep_taxa):
        raise ValueError("provide exactly one of drop_taxa or keep_taxa for pruning comparison")
    baseline = _build_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    tree = load_tree(tree_path)
    kept_taxa = sorted(set(keep_taxa or tree.tip_names) - set(drop_taxa or []))
    pruned_tree, pruning_report = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    pruned_tree_path = Path(tempfile.mkstemp(prefix="bijux-comparative-pruned-", suffix=".nwk")[1])
    pruned_traits_path = Path(tempfile.mkstemp(prefix="bijux-comparative-pruned-", suffix=Path(traits_path).suffix)[1])
    try:
        write_newick(pruned_tree_path, pruned_tree)
        _write_pruned_trait_table(
            traits_path,
            pruned_traits_path,
            kept_taxa=pruning_report.kept_taxa,
            taxon_column=taxon_column,
        )
        pruned = _build_snapshot(
            pruned_tree_path,
            pruned_traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
    finally:
        pruned_tree_path.unlink(missing_ok=True)
        pruned_traits_path.unlink(missing_ok=True)
    return ComparativePruningComparisonReport(
        tree_path=tree_path,
        response=baseline.response,
        predictors=baseline.formula.predictors,
        baseline_taxa=baseline.pgls_model.taxa,
        pruned_taxa=pruned.pgls_model.taxa,
        dropped_taxa=sorted(set(baseline.pgls_model.taxa) - set(pruned.pgls_model.taxa)),
        delta_blombergs_k=pruned.signal_k.k - baseline.signal_k.k,
        delta_pagels_lambda=pruned.signal_lambda.lambda_value - baseline.signal_lambda.lambda_value,
        coefficient_deltas=_coefficient_delta_rows(baseline.pgls_model, pruned.pgls_model),
        baseline_selected_model=baseline.model_comparison.better_model,
        pruned_selected_model=pruned.model_comparison.better_model,
    )


def _build_snapshot(
    tree_path: Path,
    traits_path: Path,
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
    signal_k = compute_blombergs_k(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
    signal_lambda = estimate_pagels_lambda(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
    contrasts = compute_phylogenetic_independent_contrasts(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
    brownian = fit_brownian_motion_model(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
    ou = fit_ornstein_uhlenbeck_model(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
    model_comparison = compare_brownian_and_ou_models(tree_path, traits_path, trait=pgls_inputs.response, taxon_column=taxon_column)
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
    limitations = _build_limitations(readiness, pgls_inputs, pgls_model, brownian, ou)
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
        audit_rows=_build_audit_rows(readiness, pgls_inputs, pgls_model, brownian, ou),
        limitations=limitations,
    )


def _build_audit_rows(
    readiness: ComparativeReadinessReport,
    pgls_inputs: PGLSInputReport,
    pgls_model: PGLSResult,
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
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
            taxa_used=ou.taxa,
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
            warnings=pgls_inputs.warnings + [f"residual_mean={pgls_model.diagnostics.residual_mean:.6f}"],
        ),
    ]


def _build_limitations(
    readiness: ComparativeReadinessReport,
    pgls_inputs: PGLSInputReport,
    pgls_model: PGLSResult,
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
) -> list[str]:
    limitations = [
        "comparative conclusions are conditioned on a single observed trait table and one supplied phylogeny",
        "causal interpretation is not warranted from comparative association alone",
    ]
    limitations.extend(readiness.warnings)
    limitations.extend(brownian.residual_diagnostics.warnings)
    limitations.extend(ou.residual_diagnostics.warnings)
    limitations.extend(warning.message for warning in ou.identifiability_warnings)
    if pgls_inputs.residual_degrees_of_freedom <= 1:
        limitations.append("PGLS residual degrees of freedom are minimal, so coefficient uncertainty is fragile")
    if pgls_inputs.categorical_predictors:
        limitations.append("categorical predictors are interpreted relative to explicit reference levels")
    if pgls_model.lambda_value in {0.0, 1.0}:
        limitations.append("estimated or supplied lambda lies on the boundary of the supported search interval")
    return sorted(set(limitations))


def _coefficient_delta_rows(left: PGLSResult, right: PGLSResult) -> list[ComparativeCoefficientDeltaRow]:
    left_coefficients = {row.name: row.estimate for row in left.coefficients}
    right_coefficients = {row.name: row.estimate for row in right.coefficients}
    rows: list[ComparativeCoefficientDeltaRow] = []
    for term in sorted(set(left_coefficients) | set(right_coefficients)):
        if term == "intercept":
            continue
        left_estimate = left_coefficients.get(term)
        right_estimate = right_coefficients.get(term)
        delta = None
        if left_estimate is not None and right_estimate is not None:
            delta = right_estimate - left_estimate
        rows.append(
            ComparativeCoefficientDeltaRow(
                term=term,
                left_estimate=left_estimate,
                right_estimate=right_estimate,
                delta=delta,
            )
        )
    return rows


def _write_pruned_trait_table(
    source_path: Path,
    out_path: Path,
    *,
    kept_taxa: list[str],
    taxon_column: str | None,
) -> Path:
    table = load_taxon_table(source_path, taxon_column=taxon_column)
    kept = set(kept_taxa)
    rows = [row for row in table.rows if row[table.taxon_column] in kept]
    return write_taxon_rows(out_path, columns=table.columns, rows=rows)
