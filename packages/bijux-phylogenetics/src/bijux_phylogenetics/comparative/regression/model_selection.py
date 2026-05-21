from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.pgls.formula import (
    coerce_numeric_value,
    parse_term_descriptor,
    resolve_formula_specification,
)
from bijux_phylogenetics.comparative.pgls.models import (
    ComparativeFormulaSpecification,
)
from bijux_phylogenetics.datasets.study_inputs import (
    load_taxon_table,
    validate_traits_table,
    write_taxon_rows,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .logistic import summarize_phylogenetic_logistic


@dataclass(slots=True)
class ComparativeRegressionModelExclusion:
    """One taxon excluded from shared comparative model selection."""

    taxon: str
    reason: str
    missing_columns: list[str]


@dataclass(slots=True)
class ComparativeRegressionModelRow:
    """One ranked candidate model in a comparative model-selection report."""

    formula: str
    model_family: str
    parameter_count: int
    taxon_count: int
    phylogenetic_parameter_name: str
    phylogenetic_parameter_value: float
    phylogenetic_parameter_estimated: bool
    log_likelihood: float
    aic: float
    aicc: float
    bic: float
    delta_aicc: float
    delta_bic: float
    akaike_weight: float
    rank: int
    selected: bool
    encoded_columns: list[str]
    warning_count: int
    separation_detected: bool


@dataclass(slots=True)
class ComparativeRegressionPairwiseComparisonRow:
    """One pairwise comparison between two comparative regression candidates."""

    left_formula: str
    right_formula: str
    comparison_kind: str
    preferred_formula: str
    left_rank: int
    right_rank: int
    left_parameter_count: int
    right_parameter_count: int
    delta_parameter_count: int
    left_log_likelihood: float
    right_log_likelihood: float
    left_aicc: float
    right_aicc: float
    left_bic: float
    right_bic: float
    likelihood_ratio_statistic: float | None


@dataclass(slots=True)
class ComparativeRegressionModelSelectionReport:
    """Ranked comparative regression hypotheses on one shared taxon set."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    response: str
    model_family: str
    selected_criterion: str
    analysis_taxa: list[str]
    excluded_taxa: list[ComparativeRegressionModelExclusion]
    rows: list[ComparativeRegressionModelRow]
    pairwise_rows: list[ComparativeRegressionPairwiseComparisonRow]
    best_formula: str
    warnings: list[str]


def compare_comparative_regression_models(
    tree_path: Path,
    traits_path: Path,
    *,
    formulas: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeRegressionModelSelectionReport:
    """Compare multiple comparative regression hypotheses on one shared taxon set."""
    if len(formulas) < 2:
        raise ComparativeMethodError(
            "comparative model selection requires at least two candidate formulas"
        )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    trait_report = validate_traits_table(traits_path, taxon_column=taxon_column)
    column_kinds = {column.name: column.kind for column in trait_report.trait_columns}
    candidate_specs = [_candidate_specification(formula) for formula in formulas]
    response_columns = {candidate.response_column for candidate in candidate_specs}
    if len(response_columns) != 1:
        raise ComparativeMethodError(
            "comparative model selection requires every candidate formula to use the same response column"
        )
    response_column = next(iter(response_columns))
    missing_columns = sorted(
        {
            column
            for candidate in candidate_specs
            for column in candidate.required_columns
            if column not in table.columns
        }
    )
    if missing_columns:
        raise ComparativeMethodError(
            "trait table does not contain required model-selection columns: "
            + ", ".join(missing_columns)
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    table_taxa = set(table.taxa)
    excluded_taxa: list[ComparativeRegressionModelExclusion] = []
    for taxon in sorted(tree_taxa - table_taxa):
        excluded_taxa.append(
            ComparativeRegressionModelExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
                missing_columns=sorted(
                    {
                        column
                        for candidate in candidate_specs
                        for column in candidate.required_columns
                    }
                ),
            )
        )
    for taxon in sorted(table_taxa - tree_taxa):
        excluded_taxa.append(
            ComparativeRegressionModelExclusion(
                taxon=taxon,
                reason="missing_from_tree",
                missing_columns=[],
            )
        )

    overlap_taxa = sorted(tree_taxa & table_taxa)
    analysis_taxa: list[str] = []
    for taxon in overlap_taxa:
        row = rows_by_taxon[taxon]
        missing_columns_for_taxon: set[str] = set()
        invalid_numeric_columns: set[str] = set()
        for candidate in candidate_specs:
            for column in candidate.required_columns:
                if not row.get(column, ""):
                    missing_columns_for_taxon.add(column)
            if not row.get(candidate.response_column, ""):
                continue
            try:
                coerce_numeric_value(
                    row[candidate.response_column],
                    descriptor=candidate.response_descriptor,
                )
            except ValueError:
                invalid_numeric_columns.add(candidate.response_column)
            for predictor_descriptor in candidate.predictor_descriptors:
                if not row.get(predictor_descriptor.source_column, ""):
                    continue
                predictor_kind = column_kinds.get(predictor_descriptor.source_column)
                if predictor_kind is None:
                    continue
                if (
                    predictor_descriptor.transformation is not None
                    and predictor_kind != "numeric"
                ):
                    invalid_numeric_columns.add(predictor_descriptor.source_column)
                    continue
                if predictor_kind != "numeric":
                    continue
                try:
                    coerce_numeric_value(
                        row[predictor_descriptor.source_column],
                        descriptor=predictor_descriptor,
                    )
                except ValueError:
                    invalid_numeric_columns.add(predictor_descriptor.source_column)
        if missing_columns_for_taxon:
            excluded_taxa.append(
                ComparativeRegressionModelExclusion(
                    taxon=taxon,
                    reason="missing_required_values",
                    missing_columns=sorted(missing_columns_for_taxon),
                )
            )
            continue
        if invalid_numeric_columns:
            excluded_taxa.append(
                ComparativeRegressionModelExclusion(
                    taxon=taxon,
                    reason="invalid_numeric_values",
                    missing_columns=sorted(invalid_numeric_columns),
                )
            )
            continue
        analysis_taxa.append(taxon)

    if len(analysis_taxa) < 3:
        raise ComparativeMethodError(
            "comparative model selection does not retain enough shared complete-case taxa"
        )
    model_family = _shared_response_family(
        [rows_by_taxon[taxon][response_column] for taxon in analysis_taxa]
    )
    if model_family == "logistic" and lambda_value == "estimate":
        raise ComparativeMethodError(
            "phylogenetic logistic model selection requires a numeric lambda value"
        )

    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, analysis_taxa)
    reduced_rows = [rows_by_taxon[taxon] for taxon in analysis_taxa]
    reduced_warnings: list[str] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-model-selection-"
    ) as tmp_dir:
        reduced_tree_path = Path(tmp_dir) / "comparative-model-selection-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "comparative-model-selection-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n",
            encoding="utf-8",
        )
        write_taxon_rows(
            reduced_table_path,
            columns=table.columns,
            rows=reduced_rows,
        )
        fitted_rows = [
            _fit_candidate_row(
                reduced_tree_path=reduced_tree_path,
                reduced_table_path=reduced_table_path,
                taxon_column=table.taxon_column,
                analysis_taxa=analysis_taxa,
                candidate=candidate,
                model_family=model_family,
                lambda_value=lambda_value,
                warnings=reduced_warnings,
            )
            for candidate in candidate_specs
        ]
    selected_criterion = (
        "AICc" if any(math.isfinite(row.aicc) for row in fitted_rows) else "AIC"
    )
    _rank_rows(fitted_rows, selected_criterion=selected_criterion)
    pairwise_rows = _build_pairwise_rows(fitted_rows)
    selected_row = next(row for row in fitted_rows if row.selected)
    warnings = list(reduced_warnings)
    if selected_row.separation_detected:
        warnings.append(
            f"selected model '{selected_row.formula}' carries logistic separation-risk warnings"
        )
    return ComparativeRegressionModelSelectionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=response_column,
        model_family=model_family,
        selected_criterion=selected_criterion,
        analysis_taxa=analysis_taxa,
        excluded_taxa=excluded_taxa,
        rows=fitted_rows,
        pairwise_rows=pairwise_rows,
        best_formula=selected_row.formula,
        warnings=warnings,
    )


def write_comparative_regression_model_ranking_table(
    path: Path, report: ComparativeRegressionModelSelectionReport
) -> Path:
    """Write one ranked comparative model table."""
    return write_taxon_rows(
        path,
        columns=[
            "formula",
            "model_family",
            "parameter_count",
            "taxon_count",
            "phylogenetic_parameter_name",
            "phylogenetic_parameter_value",
            "phylogenetic_parameter_estimated",
            "log_likelihood",
            "aic",
            "aicc",
            "bic",
            "delta_aicc",
            "delta_bic",
            "akaike_weight",
            "rank",
            "selected",
            "encoded_columns",
            "warning_count",
            "separation_detected",
        ],
        rows=[
            {
                "formula": row.formula,
                "model_family": row.model_family,
                "parameter_count": str(row.parameter_count),
                "taxon_count": str(row.taxon_count),
                "phylogenetic_parameter_name": row.phylogenetic_parameter_name,
                "phylogenetic_parameter_value": format(
                    row.phylogenetic_parameter_value, ".15g"
                ),
                "phylogenetic_parameter_estimated": str(
                    row.phylogenetic_parameter_estimated
                ).lower(),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "aic": format(row.aic, ".15g"),
                "aicc": "inf" if math.isinf(row.aicc) else format(row.aicc, ".15g"),
                "bic": format(row.bic, ".15g"),
                "delta_aicc": "inf"
                if math.isinf(row.delta_aicc)
                else format(row.delta_aicc, ".15g"),
                "delta_bic": format(row.delta_bic, ".15g"),
                "akaike_weight": format(row.akaike_weight, ".15g"),
                "rank": str(row.rank),
                "selected": str(row.selected).lower(),
                "encoded_columns": ",".join(row.encoded_columns),
                "warning_count": str(row.warning_count),
                "separation_detected": str(row.separation_detected).lower(),
            }
            for row in report.rows
        ],
    )


def write_comparative_regression_pairwise_table(
    path: Path, report: ComparativeRegressionModelSelectionReport
) -> Path:
    """Write one pairwise nested-versus-non-nested comparison table."""
    return write_taxon_rows(
        path,
        columns=[
            "left_formula",
            "right_formula",
            "comparison_kind",
            "preferred_formula",
            "left_rank",
            "right_rank",
            "left_parameter_count",
            "right_parameter_count",
            "delta_parameter_count",
            "left_log_likelihood",
            "right_log_likelihood",
            "left_aicc",
            "right_aicc",
            "left_bic",
            "right_bic",
            "likelihood_ratio_statistic",
        ],
        rows=[
            {
                "left_formula": row.left_formula,
                "right_formula": row.right_formula,
                "comparison_kind": row.comparison_kind,
                "preferred_formula": row.preferred_formula,
                "left_rank": str(row.left_rank),
                "right_rank": str(row.right_rank),
                "left_parameter_count": str(row.left_parameter_count),
                "right_parameter_count": str(row.right_parameter_count),
                "delta_parameter_count": str(row.delta_parameter_count),
                "left_log_likelihood": format(row.left_log_likelihood, ".15g"),
                "right_log_likelihood": format(row.right_log_likelihood, ".15g"),
                "left_aicc": "inf"
                if math.isinf(row.left_aicc)
                else format(row.left_aicc, ".15g"),
                "right_aicc": "inf"
                if math.isinf(row.right_aicc)
                else format(row.right_aicc, ".15g"),
                "left_bic": format(row.left_bic, ".15g"),
                "right_bic": format(row.right_bic, ".15g"),
                "likelihood_ratio_statistic": ""
                if row.likelihood_ratio_statistic is None
                else format(row.likelihood_ratio_statistic, ".15g"),
            }
            for row in report.pairwise_rows
        ],
    )


def write_comparative_regression_excluded_taxa_table(
    path: Path, report: ComparativeRegressionModelSelectionReport
) -> Path:
    """Write one explicit excluded-taxa ledger for comparative model selection."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "missing_columns"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "missing_columns": ",".join(row.missing_columns),
            }
            for row in report.excluded_taxa
        ],
    )


@dataclass(slots=True)
class _CandidateSpecification:
    formula: str
    specification: ComparativeFormulaSpecification
    response_column: str
    response_descriptor: object
    required_columns: list[str]
    predictor_descriptors: list[object]


def _candidate_specification(formula: str) -> _CandidateSpecification:
    specification = resolve_formula_specification(
        response=None,
        predictors=None,
        formula=formula,
    )
    response_descriptor = parse_term_descriptor(specification.response)
    predictor_descriptors = [
        parse_term_descriptor(term) for term in specification.predictors
    ]
    required_columns = [response_descriptor.source_column]
    required_columns.extend(
        descriptor.source_column for descriptor in predictor_descriptors
    )
    return _CandidateSpecification(
        formula=formula,
        specification=specification,
        response_column=response_descriptor.source_column,
        response_descriptor=response_descriptor,
        required_columns=sorted(set(required_columns)),
        predictor_descriptors=predictor_descriptors,
    )


def _shared_response_family(response_values: list[str]) -> str:
    numeric_values = [float(value) for value in response_values]
    if all(
        math.isclose(value, round(value), abs_tol=1e-12) for value in numeric_values
    ) and {int(round(value)) for value in numeric_values} <= {0, 1}:
        return "logistic"
    return "pgls"


def _fit_candidate_row(
    *,
    reduced_tree_path: Path,
    reduced_table_path: Path,
    taxon_column: str,
    analysis_taxa: list[str],
    candidate: _CandidateSpecification,
    model_family: str,
    lambda_value: float | str,
    warnings: list[str],
) -> ComparativeRegressionModelRow:
    if model_family == "logistic":
        logistic_report = summarize_phylogenetic_logistic(
            reduced_tree_path,
            reduced_table_path,
            formula=candidate.formula,
            taxon_column=taxon_column,
            lambda_value=float(lambda_value),
        )
        for warning in logistic_report.warnings:
            warnings.append(f"{candidate.formula}: {warning.message}")
        parameter_count = len(logistic_report.coefficients)
        log_likelihood = logistic_report.binomial_log_likelihood
        aic, aicc, bic = _information_criteria(
            log_likelihood=log_likelihood,
            parameter_count=parameter_count,
            sample_size=len(analysis_taxa),
        )
        return ComparativeRegressionModelRow(
            formula=candidate.formula,
            model_family="logistic",
            parameter_count=parameter_count,
            taxon_count=len(analysis_taxa),
            phylogenetic_parameter_name="lambda",
            phylogenetic_parameter_value=logistic_report.lambda_value,
            phylogenetic_parameter_estimated=False,
            log_likelihood=log_likelihood,
            aic=aic,
            aicc=aicc,
            bic=bic,
            delta_aicc=0.0,
            delta_bic=0.0,
            akaike_weight=0.0,
            rank=0,
            selected=False,
            encoded_columns=[row.name for row in logistic_report.coefficients],
            warning_count=len(logistic_report.warnings),
            separation_detected=logistic_report.separation_detected,
        )
    pgls_report = run_pgls(
        reduced_tree_path,
        reduced_table_path,
        formula=candidate.formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    parameter_count = len(pgls_report.coefficients) + 1
    if pgls_report.lambda_fit.mode == "estimated":
        parameter_count += 1
    aic, aicc, bic = _information_criteria(
        log_likelihood=pgls_report.log_likelihood,
        parameter_count=parameter_count,
        sample_size=len(analysis_taxa),
    )
    return ComparativeRegressionModelRow(
        formula=candidate.formula,
        model_family="pgls",
        parameter_count=parameter_count,
        taxon_count=len(analysis_taxa),
        phylogenetic_parameter_name="lambda",
        phylogenetic_parameter_value=pgls_report.lambda_value,
        phylogenetic_parameter_estimated=pgls_report.lambda_fit.mode == "estimated",
        log_likelihood=pgls_report.log_likelihood,
        aic=aic,
        aicc=aicc,
        bic=bic,
        delta_aicc=0.0,
        delta_bic=0.0,
        akaike_weight=0.0,
        rank=0,
        selected=False,
        encoded_columns=list(pgls_report.encoded_columns),
        warning_count=0,
        separation_detected=False,
    )


def _information_criteria(
    *, log_likelihood: float, parameter_count: int, sample_size: int
) -> tuple[float, float, float]:
    aic = (2.0 * parameter_count) - (2.0 * log_likelihood)
    if sample_size <= parameter_count + 1:
        aicc = math.inf
    else:
        aicc = aic + (
            (2.0 * parameter_count * (parameter_count + 1))
            / (sample_size - parameter_count - 1)
        )
    bic = (math.log(sample_size) * parameter_count) - (2.0 * log_likelihood)
    return aic, aicc, bic


def _rank_rows(
    rows: list[ComparativeRegressionModelRow], *, selected_criterion: str
) -> None:
    if selected_criterion == "AICc":
        best_score = min(row.aicc for row in rows)
        for row in rows:
            row.delta_aicc = row.aicc - best_score
    else:
        best_score = min(row.aic for row in rows)
        for row in rows:
            row.delta_aicc = row.aic - best_score
    best_bic = min(row.bic for row in rows)
    for row in rows:
        row.delta_bic = row.bic - best_bic
        row.selected = False
    if selected_criterion == "AICc":
        raw_weights = [
            0.0 if math.isinf(row.delta_aicc) else math.exp(-0.5 * row.delta_aicc)
            for row in rows
        ]
    else:
        raw_weights = [math.exp(-0.5 * row.delta_aicc) for row in rows]
    weight_total = sum(raw_weights)
    for row, raw_weight in zip(rows, raw_weights, strict=True):
        row.akaike_weight = raw_weight / weight_total if weight_total else 0.0
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            row.delta_aicc,
            row.delta_bic,
            -row.log_likelihood,
            row.formula,
        ),
    )
    for rank, row in enumerate(sorted_rows, start=1):
        row.rank = rank
    sorted_rows[0].selected = True


def _build_pairwise_rows(
    rows: list[ComparativeRegressionModelRow],
) -> list[ComparativeRegressionPairwiseComparisonRow]:
    pairwise_rows: list[ComparativeRegressionPairwiseComparisonRow] = []
    for left_index, left in enumerate(rows):
        for right in rows[left_index + 1 :]:
            comparison_kind = _comparison_kind(
                left.encoded_columns, right.encoded_columns
            )
            preferred = left.formula if left.rank < right.rank else right.formula
            likelihood_ratio_statistic: float | None = None
            if comparison_kind in {"left_nested_in_right", "right_nested_in_left"}:
                likelihood_ratio_statistic = 2.0 * abs(
                    right.log_likelihood - left.log_likelihood
                )
            pairwise_rows.append(
                ComparativeRegressionPairwiseComparisonRow(
                    left_formula=left.formula,
                    right_formula=right.formula,
                    comparison_kind=comparison_kind,
                    preferred_formula=preferred,
                    left_rank=left.rank,
                    right_rank=right.rank,
                    left_parameter_count=left.parameter_count,
                    right_parameter_count=right.parameter_count,
                    delta_parameter_count=abs(
                        right.parameter_count - left.parameter_count
                    ),
                    left_log_likelihood=left.log_likelihood,
                    right_log_likelihood=right.log_likelihood,
                    left_aicc=left.aicc,
                    right_aicc=right.aicc,
                    left_bic=left.bic,
                    right_bic=right.bic,
                    likelihood_ratio_statistic=likelihood_ratio_statistic,
                )
            )
    return pairwise_rows


def _comparison_kind(left_columns: list[str], right_columns: list[str]) -> str:
    left_set = set(left_columns)
    right_set = set(right_columns)
    if left_set == right_set:
        return "identical"
    if left_set < right_set:
        return "left_nested_in_right"
    if right_set < left_set:
        return "right_nested_in_left"
    return "non_nested"
