from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.datasets.study_inputs import (
    load_taxon_table,
    validate_traits_table,
    write_taxon_rows,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .formula import (
    coerce_numeric_value,
    column_supports_numeric_values,
    interaction_column_names,
    interaction_values,
    parse_term_descriptor,
    resolve_formula_specification,
)
from .models import (
    ComparativeFormulaSpecification,
    PGLSFormulaAudit,
    PGLSInputReport,
    PGLSInteractionAudit,
    PGLSModelMatrixReport,
    PGLSModelMatrixRow,
    PGLSPredictorClassification,
    PGLSTaxonExclusion,
    _FormulaTermDescriptor,
)


def inspect_pgls_inputs(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
) -> PGLSInputReport:
    """Inspect whether a PGLS request is valid for the given tree and trait table."""
    specification = resolve_formula_specification(
        response=response,
        predictors=predictors,
        formula=formula,
    )
    response_descriptor = parse_term_descriptor(specification.response)
    if response_descriptor.transformation is not None:
        raise ComparativeMethodError(
            "transformed response terms are not supported for PGLS"
        )
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=response_descriptor.source_column,
        taxon_column=taxon_column,
    )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    trait_report = validate_traits_table(traits_path, taxon_column=taxon_column)
    column_kinds = {column.name: column.kind for column in trait_report.trait_columns}
    blockers: list[str] = []
    warnings = list(readiness.warnings)

    if not readiness.rooted:
        blockers.append("PGLS requires a rooted tree")
    if not readiness.complete_branch_lengths:
        blockers.append("PGLS requires complete tree branch lengths")
    if response_descriptor.source_column not in column_kinds:
        blockers.append(
            f"trait table does not contain response column '{response_descriptor.source_column}'"
        )
    elif not column_supports_numeric_values(
        table.rows,
        response_descriptor.source_column,
    ):
        blockers.append(
            f"response column '{response_descriptor.source_column}' must be numeric for PGLS"
        )

    predictor_reports: list[PGLSPredictorClassification] = []
    categorical_predictors: list[str] = []
    encoded_columns = ["intercept"] if specification.include_intercept else []
    for predictor in specification.predictors:
        term_descriptor = parse_term_descriptor(predictor)
        kind = column_kinds.get(term_descriptor.source_column)
        numeric_candidate = column_supports_numeric_values(
            table.rows,
            term_descriptor.source_column,
        )
        if kind is None:
            blockers.append(
                f"trait table does not contain predictor column '{term_descriptor.source_column}'"
            )
            continue
        if not numeric_candidate and term_descriptor.transformation is not None:
            blockers.append(
                f"transformation '{term_descriptor.transformation}' requires numeric predictor column '{term_descriptor.source_column}'"
            )
            continue
        if numeric_candidate:
            if term_descriptor.transformation is not None:
                warnings.append(
                    f"predictor term '{predictor}' applies {term_descriptor.transformation} transformation to column '{term_descriptor.source_column}'"
                )
            predictor_reports.append(
                PGLSPredictorClassification(
                    name=predictor,
                    kind="numeric",
                    raw_term=predictor,
                    source_column=term_descriptor.source_column,
                    transformation=term_descriptor.transformation,
                )
            )
            encoded_columns.append(predictor)
            continue

        categorical_predictors.append(predictor)
        levels = sorted(
            {
                row[term_descriptor.source_column]
                for row in table.rows
                if row[table.taxon_column] in tree.tip_names
                and row.get(term_descriptor.source_column, "")
            }
        )
        level_counts = {
            level: sum(
                1
                for row in table.rows
                if row[table.taxon_column] in tree.tip_names
                and row.get(term_descriptor.source_column, "") == level
            )
            for level in levels
        }
        if len(levels) < 2:
            blockers.append(
                f"categorical predictor '{term_descriptor.source_column}' requires at least two observed levels"
            )
            predictor_reports.append(
                PGLSPredictorClassification(
                    name=predictor,
                    kind=kind,
                    raw_term=predictor,
                    source_column=term_descriptor.source_column,
                    reference_level=levels[0] if levels else None,
                    encoded_columns=[],
                    observed_levels=levels,
                    level_counts=level_counts,
                )
            )
            continue
        reference_level = levels[0] if specification.include_intercept else None
        dummy_levels = levels[1:] if specification.include_intercept else levels
        dummy_columns = [f"{predictor}[{level}]" for level in dummy_levels]
        encoded_columns.extend(dummy_columns)
        if specification.include_intercept:
            warnings.append(
                f"categorical predictor '{predictor}' will be dummy-encoded with reference level '{reference_level}'"
            )
        else:
            warnings.append(
                f"categorical predictor '{predictor}' will be fully indicator-encoded because the formula excludes an intercept"
            )
        predictor_reports.append(
            PGLSPredictorClassification(
                name=predictor,
                kind=kind,
                raw_term=predictor,
                source_column=term_descriptor.source_column,
                reference_level=reference_level,
                encoded_columns=dummy_columns,
                observed_levels=levels,
                level_counts=level_counts,
            )
        )

    report_by_name = {report.name: report for report in predictor_reports}
    interaction_audits: list[PGLSInteractionAudit] = []
    for interaction in specification.interaction_terms:
        factor_names = interaction.split(":")
        missing_factors = [name for name in factor_names if name not in report_by_name]
        if missing_factors:
            blockers.append(
                f"interaction term '{interaction}' references unknown predictor(s): {', '.join(missing_factors)}"
            )
            continue
        interaction_columns = interaction_column_names(interaction, report_by_name)
        encoded_columns.extend(interaction_columns)
        interaction_audits.append(
            PGLSInteractionAudit(
                term=interaction,
                component_terms=factor_names,
                encoded_columns=interaction_columns,
            )
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    analysis_taxa: list[str] = []
    missing_tree_taxa: list[str] = []
    excluded_taxa: list[PGLSTaxonExclusion] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            missing_tree_taxa.append(taxon)
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="missing_from_trait_table",
                    details="taxon is present in the tree but absent from the trait table",
                )
            )
            continue
        required_columns = [response_descriptor.source_column] + [
            report.source_column or report.name for report in predictor_reports
        ]
        missing_columns = [
            column for column in required_columns if not row.get(column, "")
        ]
        if missing_columns:
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="missing_value",
                    details=f"taxon is missing required value(s): {', '.join(sorted(set(missing_columns)))}",
                )
            )
            continue
        try:
            coerce_numeric_value(
                row[response_descriptor.source_column], descriptor=response_descriptor
            )
            for predictor_report in predictor_reports:
                if predictor_report.kind == "numeric":
                    source_column = (
                        predictor_report.source_column or predictor_report.name
                    )
                    coerce_numeric_value(
                        row[source_column],
                        descriptor=parse_term_descriptor(
                            predictor_report.raw_term or predictor_report.name
                        ),
                    )
        except ValueError:
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="non_numeric_or_invalid_value",
                    details="taxon has non-numeric or transformation-invalid value(s) required by the model",
                )
            )
            continue
        analysis_taxa.append(taxon)
    if missing_tree_taxa:
        blockers.append(
            "PGLS requires all analyzed taxa to be resolved against the trait table"
        )
    residual_degrees_of_freedom = len(analysis_taxa) - len(encoded_columns)
    if residual_degrees_of_freedom <= 0:
        blockers.append(
            "PGLS overfit guard requires at least one residual degree of freedom after predictor encoding"
        )
    elif residual_degrees_of_freedom == 1:
        warnings.append(
            "PGLS has only one residual degree of freedom after predictor encoding"
        )
    transformed_terms = sorted(
        [
            predictor_report.name
            for predictor_report in predictor_reports
            if predictor_report.transformation is not None
        ]
    )
    formula_audit = PGLSFormulaAudit(
        response_term=specification.response,
        response_column=response_descriptor.source_column,
        predictor_terms=predictor_reports,
        interaction_terms=interaction_audits,
        transformed_terms=transformed_terms,
        excluded_taxa=excluded_taxa,
        includes_intercept=specification.include_intercept,
        encoded_columns=encoded_columns,
        analysis_taxa=analysis_taxa,
        parameter_count=len(encoded_columns),
        minimum_required_taxa=len(encoded_columns) + 1,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        overfit_guard_triggered=residual_degrees_of_freedom <= 0,
        warnings=warnings,
    )
    model_matrix = _build_model_matrix_report(
        rows_by_taxon=rows_by_taxon,
        taxa=analysis_taxa,
        specification=specification,
        predictor_reports=predictor_reports,
        response_descriptor=response_descriptor,
        response_column=response_descriptor.source_column,
    )

    return PGLSInputReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=specification.response,
        formula=specification,
        predictors=predictor_reports,
        formula_audit=formula_audit,
        categorical_predictors=categorical_predictors,
        encoded_columns=encoded_columns,
        analysis_taxa=analysis_taxa,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        model_matrix=model_matrix,
        ready=not blockers,
        blockers=blockers,
        warnings=warnings,
    )


def build_pgls_model_matrix(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
) -> PGLSModelMatrixReport:
    """Build the encoded design matrix implied by one PGLS request."""
    return inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    ).model_matrix


def write_pgls_model_matrix_table(path: Path, report: PGLSModelMatrixReport) -> Path:
    """Write a comparative model matrix as CSV or TSV."""
    return write_taxon_rows(
        path,
        columns=["taxon", "response_value", *report.encoded_columns],
        rows=[
            {
                "taxon": row.taxon,
                "response_value": format(row.response_value, ".15g"),
                **{
                    column: format(row.encoded_values[column], ".15g")
                    for column in report.encoded_columns
                },
            }
            for row in report.rows
        ],
    )


def _build_design_matrix(
    rows_by_taxon: dict[str, dict[str, str]],
    taxa: list[str],
    predictors: list[str],
    predictor_reports: list[PGLSPredictorClassification],
    interaction_terms: list[str],
    *,
    include_intercept: bool,
) -> tuple[list[list[float]], list[str]]:
    encoded_columns = ["intercept"] if include_intercept else []
    report_by_name = {report.name: report for report in predictor_reports}
    for predictor in predictors:
        report = report_by_name[predictor]
        if report.kind == "numeric":
            encoded_columns.append(predictor)
        else:
            encoded_columns.extend(report.encoded_columns or [])
    interaction_columns = {
        interaction: interaction_column_names(interaction, report_by_name)
        for interaction in interaction_terms
    }
    for columns in interaction_columns.values():
        encoded_columns.extend(columns)
    matrix: list[list[float]] = []
    for taxon in taxa:
        row = rows_by_taxon[taxon]
        encoded_row = [1.0] if include_intercept else []
        encoded_main_effects: dict[str, list[tuple[str, float]]] = {}
        for predictor in predictors:
            report = report_by_name[predictor]
            if report.kind == "numeric":
                source_column = report.source_column or predictor
                numeric_value = coerce_numeric_value(
                    row[source_column],
                    descriptor=parse_term_descriptor(report.raw_term or predictor),
                )
                encoded_row.append(numeric_value)
                encoded_main_effects[predictor] = [(predictor, numeric_value)]
                continue
            categorical_rows: list[tuple[str, float]] = []
            for encoded_name in report.encoded_columns or []:
                level = encoded_name.removeprefix(f"{predictor}[").removesuffix("]")
                source_column = report.source_column or predictor
                value = 1.0 if row[source_column] == level else 0.0
                categorical_rows.append((encoded_name, value))
                encoded_row.append(value)
            encoded_main_effects[predictor] = categorical_rows
        for interaction in interaction_terms:
            encoded_row.extend(
                value
                for _, value in interaction_values(interaction, encoded_main_effects)
            )
        matrix.append(encoded_row)
    return matrix, encoded_columns


def _build_model_matrix_report(
    *,
    rows_by_taxon: dict[str, dict[str, str]],
    taxa: list[str],
    specification: ComparativeFormulaSpecification,
    predictor_reports: list[PGLSPredictorClassification],
    response_descriptor: _FormulaTermDescriptor,
    response_column: str,
) -> PGLSModelMatrixReport:
    design_matrix, encoded_columns = _build_design_matrix(
        rows_by_taxon,
        taxa,
        specification.predictors,
        predictor_reports,
        specification.interaction_terms,
        include_intercept=specification.include_intercept,
    )
    rows = [
        PGLSModelMatrixRow(
            taxon=taxon,
            response_value=coerce_numeric_value(
                rows_by_taxon[taxon][response_column],
                descriptor=response_descriptor,
            ),
            encoded_values={
                column: design_matrix[row_index][column_index]
                for column_index, column in enumerate(encoded_columns)
            },
        )
        for row_index, taxon in enumerate(taxa)
    ]
    return PGLSModelMatrixReport(
        formula=specification,
        response_column=response_column,
        encoded_columns=encoded_columns,
        row_count=len(rows),
        rows=rows,
    )
