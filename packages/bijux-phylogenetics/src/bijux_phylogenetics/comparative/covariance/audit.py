from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    matrix_condition_number,
    matrix_rank,
    stable_covariance,
)
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    build_ou_covariance_matrix,
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import _alpha_grid
from bijux_phylogenetics.comparative.pgls import inspect_pgls_inputs
from bijux_phylogenetics.datasets.study_inputs import (
    inspect_taxon_table_index,
    write_taxon_rows,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

COVARIANCE_AUDIT_CONDITION_THRESHOLD = 1e12
COVARIANCE_AUDIT_REGULARIZATION_EPSILON = 1e-8


@dataclass(slots=True)
class CovarianceAuditCandidateRow:
    """One covariance candidate audited before comparative fitting."""

    candidate_label: str
    parameter_name: str | None
    parameter_value: float | None
    matrix_dimension: int
    matrix_rank: int
    condition_number: float
    fit_condition_number: float | None
    positive_definite_before_fit: bool
    singular: bool
    near_singular: bool
    fit_strategy: str
    fit_strategy_details: str


@dataclass(slots=True)
class CovarianceAuditExcludedTaxon:
    """One taxon excluded from one covariance-audit analysis."""

    taxon: str
    reason: str
    details: str


@dataclass(slots=True)
class ComparativeCovarianceAuditReport:
    """Pre-fit covariance audit for one comparative analysis surface."""

    tree_path: Path
    traits_path: Path
    analysis: str
    covariance_model: str
    analysis_label: str
    taxon_column: str | None
    tree_taxon_count: int
    trait_taxon_count: int
    matched_taxa: list[str]
    missing_from_traits: list[str]
    extra_trait_taxa: list[str]
    duplicate_tree_taxa: list[str]
    duplicate_trait_taxa: list[str]
    empty_trait_taxon_rows: list[int]
    analysis_taxa: list[str]
    excluded_taxa: list[CovarianceAuditExcludedTaxon]
    matrix_dimension: int
    matrix_rank: int
    condition_number: float
    fit_strategy: str
    fit_strategy_details: str
    singular: bool
    near_singular: bool
    zero_length_branch_count: int
    negative_branch_length_count: int
    minimum_branch_length: float | None
    maximum_branch_length: float | None
    blockers: list[str]
    warnings: list[str]
    candidate_rows: list[CovarianceAuditCandidateRow]


def summarize_comparative_covariance_audit(
    tree_path: Path,
    traits_path: Path,
    *,
    analysis: str,
    trait: str | None = None,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
    alpha: float | str = "estimate",
) -> ComparativeCovarianceAuditReport:
    """Audit comparative covariance behavior before model coefficients are trusted."""
    if analysis == "pgls":
        return _summarize_pgls_covariance_audit(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
    if analysis == "brownian-trait":
        if trait is None:
            raise ComparativeMethodError(
                "brownian-trait covariance audit requires --trait"
            )
        return _summarize_brownian_trait_covariance_audit(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
    if analysis == "ou-trait":
        if trait is None:
            raise ComparativeMethodError("ou-trait covariance audit requires --trait")
        return _summarize_ou_trait_covariance_audit(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            alpha=alpha,
        )
    raise ComparativeMethodError(
        "supported covariance audit analyses are 'pgls', 'brownian-trait', and 'ou-trait'"
    )


def write_comparative_covariance_audit_summary_table(
    path: Path, report: ComparativeCovarianceAuditReport
) -> Path:
    """Write one covariance-audit summary row as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "analysis",
            "covariance_model",
            "analysis_label",
            "matrix_dimension",
            "matrix_rank",
            "condition_number",
            "fit_strategy",
            "singular",
            "near_singular",
            "tree_taxon_count",
            "trait_taxon_count",
            "matched_taxon_count",
            "analysis_taxon_count",
            "missing_from_traits_count",
            "extra_trait_taxon_count",
            "duplicate_tree_taxon_count",
            "duplicate_trait_taxon_count",
            "empty_trait_taxon_row_count",
            "zero_length_branch_count",
            "negative_branch_length_count",
            "minimum_branch_length",
            "maximum_branch_length",
            "blocker_count",
            "warning_count",
        ],
        rows=[
            {
                "analysis": report.analysis,
                "covariance_model": report.covariance_model,
                "analysis_label": report.analysis_label,
                "matrix_dimension": report.matrix_dimension,
                "matrix_rank": report.matrix_rank,
                "condition_number": format(report.condition_number, ".15g"),
                "fit_strategy": report.fit_strategy,
                "singular": str(report.singular).lower(),
                "near_singular": str(report.near_singular).lower(),
                "tree_taxon_count": report.tree_taxon_count,
                "trait_taxon_count": report.trait_taxon_count,
                "matched_taxon_count": len(report.matched_taxa),
                "analysis_taxon_count": len(report.analysis_taxa),
                "missing_from_traits_count": len(report.missing_from_traits),
                "extra_trait_taxon_count": len(report.extra_trait_taxa),
                "duplicate_tree_taxon_count": len(report.duplicate_tree_taxa),
                "duplicate_trait_taxon_count": len(report.duplicate_trait_taxa),
                "empty_trait_taxon_row_count": len(report.empty_trait_taxon_rows),
                "zero_length_branch_count": report.zero_length_branch_count,
                "negative_branch_length_count": report.negative_branch_length_count,
                "minimum_branch_length": (
                    ""
                    if report.minimum_branch_length is None
                    else format(report.minimum_branch_length, ".15g")
                ),
                "maximum_branch_length": (
                    ""
                    if report.maximum_branch_length is None
                    else format(report.maximum_branch_length, ".15g")
                ),
                "blocker_count": len(report.blockers),
                "warning_count": len(report.warnings),
            }
        ],
    )


def write_comparative_covariance_audit_candidate_table(
    path: Path, report: ComparativeCovarianceAuditReport
) -> Path:
    """Write one candidate-level covariance audit ledger as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "candidate_label",
            "parameter_name",
            "parameter_value",
            "matrix_dimension",
            "matrix_rank",
            "condition_number",
            "fit_condition_number",
            "positive_definite_before_fit",
            "singular",
            "near_singular",
            "fit_strategy",
            "fit_strategy_details",
        ],
        rows=[
            {
                "candidate_label": row.candidate_label,
                "parameter_name": row.parameter_name or "",
                "parameter_value": (
                    ""
                    if row.parameter_value is None
                    else format(row.parameter_value, ".15g")
                ),
                "matrix_dimension": row.matrix_dimension,
                "matrix_rank": row.matrix_rank,
                "condition_number": format(row.condition_number, ".15g"),
                "fit_condition_number": (
                    ""
                    if row.fit_condition_number is None
                    else format(row.fit_condition_number, ".15g")
                ),
                "positive_definite_before_fit": str(
                    row.positive_definite_before_fit
                ).lower(),
                "singular": str(row.singular).lower(),
                "near_singular": str(row.near_singular).lower(),
                "fit_strategy": row.fit_strategy,
                "fit_strategy_details": row.fit_strategy_details,
            }
            for row in report.candidate_rows
        ],
    )


def write_comparative_covariance_audit_excluded_taxa_table(
    path: Path, report: ComparativeCovarianceAuditReport
) -> Path:
    """Write one explicit excluded-taxa ledger for covariance audit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "details"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "details": row.details,
            }
            for row in report.excluded_taxa
        ],
    )


def _summarize_pgls_covariance_audit(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
    taxon_column: str | None,
    lambda_value: float | str,
) -> ComparativeCovarianceAuditReport:
    tree_validation, tree = _validated_tree_inputs(tree_path)
    table_audit = inspect_taxon_table_index(traits_path, taxon_column=taxon_column)
    overlap = _taxon_overlap(tree, table_audit)
    blockers, warnings = _base_input_messages(tree_validation, table_audit)
    if table_audit.duplicate_taxa:
        blockers.append("trait table contains duplicate taxon keys")
    if table_audit.empty_taxon_rows:
        blockers.append("trait table contains empty taxon keys")

    analysis_taxa: list[str] = []
    excluded_taxa: list[CovarianceAuditExcludedTaxon] = []
    analysis_label = formula or (
        ""
        if response is None or not predictors
        else f"{response} ~ {' + '.join(predictors)}"
    )
    if not blockers:
        input_report = inspect_pgls_inputs(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
        )
        analysis_taxa = list(input_report.analysis_taxa)
        excluded_taxa = [
            CovarianceAuditExcludedTaxon(
                taxon=row.taxon,
                reason=row.reason,
                details=row.details,
            )
            for row in input_report.formula_audit.excluded_taxa
        ]
        blockers.extend(input_report.blockers)
        warnings.extend(input_report.warnings)
        analysis_label = input_report.formula.formula
    candidate_rows = (
        []
        if blockers
        else _pgls_candidate_rows(tree, analysis_taxa, lambda_value=lambda_value)
    )
    blockers.extend(_candidate_failure_messages(candidate_rows))
    return _finalize_covariance_audit(
        tree_path=tree_path,
        traits_path=traits_path,
        analysis="pgls",
        covariance_model="pagel-lambda",
        analysis_label=analysis_label,
        taxon_column=table_audit.taxon_column,
        tree=tree,
        tree_validation=tree_validation,
        table_audit=table_audit,
        overlap=overlap,
        analysis_taxa=analysis_taxa,
        excluded_taxa=excluded_taxa,
        blockers=blockers,
        warnings=warnings,
        candidate_rows=candidate_rows,
    )


def _summarize_brownian_trait_covariance_audit(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
) -> ComparativeCovarianceAuditReport:
    tree_validation, tree = _validated_tree_inputs(tree_path)
    table_audit = inspect_taxon_table_index(traits_path, taxon_column=taxon_column)
    overlap = _taxon_overlap(tree, table_audit)
    blockers, warnings = _base_input_messages(tree_validation, table_audit)
    if table_audit.duplicate_taxa:
        blockers.append("trait table contains duplicate taxon keys")
    if table_audit.empty_taxon_rows:
        blockers.append("trait table contains empty taxon keys")

    analysis_taxa: list[str] = []
    excluded_taxa: list[CovarianceAuditExcludedTaxon] = []
    if not blockers:
        readiness = summarize_numeric_trait_readiness(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        analysis_taxa = list(readiness.analysis_taxa)
        excluded_taxa = _trait_readiness_exclusions(readiness, trait=trait)
        blockers.extend(readiness.blockers)
        warnings.extend(readiness.warnings)
    candidate_rows = (
        []
        if blockers
        else [
            _assess_covariance_candidate(
                candidate_label="brownian",
                parameter_name=None,
                parameter_value=None,
                covariance_matrix=build_brownian_covariance_matrix(tree, analysis_taxa),
                fit_strategy="regularization",
                fit_strategy_details=(
                    "Brownian trait fitting applies diagonal stabilization epsilon=1e-8 before covariance inversion"
                ),
            )
        ]
    )
    blockers.extend(_candidate_failure_messages(candidate_rows))
    return _finalize_covariance_audit(
        tree_path=tree_path,
        traits_path=traits_path,
        analysis="brownian-trait",
        covariance_model="brownian",
        analysis_label=trait,
        taxon_column=table_audit.taxon_column,
        tree=tree,
        tree_validation=tree_validation,
        table_audit=table_audit,
        overlap=overlap,
        analysis_taxa=analysis_taxa,
        excluded_taxa=excluded_taxa,
        blockers=blockers,
        warnings=warnings,
        candidate_rows=candidate_rows,
    )


def _summarize_ou_trait_covariance_audit(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    alpha: float | str,
) -> ComparativeCovarianceAuditReport:
    tree_validation, tree = _validated_tree_inputs(tree_path)
    table_audit = inspect_taxon_table_index(traits_path, taxon_column=taxon_column)
    overlap = _taxon_overlap(tree, table_audit)
    blockers, warnings = _base_input_messages(tree_validation, table_audit)
    if table_audit.duplicate_taxa:
        blockers.append("trait table contains duplicate taxon keys")
    if table_audit.empty_taxon_rows:
        blockers.append("trait table contains empty taxon keys")

    analysis_taxa: list[str] = []
    excluded_taxa: list[CovarianceAuditExcludedTaxon] = []
    candidate_rows: list[CovarianceAuditCandidateRow] = []
    if not blockers:
        readiness = summarize_numeric_trait_readiness(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        analysis_taxa = list(readiness.analysis_taxa)
        excluded_taxa = _trait_readiness_exclusions(readiness, trait=trait)
        blockers.extend(readiness.blockers)
        warnings.extend(readiness.warnings)
        if not blockers:
            dataset = load_comparative_dataset(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                minimum_taxa=3,
                require_rooted=True,
                require_binary=False,
            )
            alpha_candidates = _resolved_alpha_candidates(dataset, alpha=alpha)
            candidate_rows = [
                _assess_covariance_candidate(
                    candidate_label=f"alpha={candidate:.6g}",
                    parameter_name="alpha",
                    parameter_value=candidate,
                    covariance_matrix=build_ou_covariance_matrix(
                        tree, analysis_taxa, alpha=candidate
                    ),
                    fit_strategy="regularization",
                    fit_strategy_details=(
                        "OU trait fitting applies diagonal stabilization epsilon=1e-8 before covariance inversion"
                    ),
                )
                for candidate in alpha_candidates
            ]
    blockers.extend(_candidate_failure_messages(candidate_rows))
    return _finalize_covariance_audit(
        tree_path=tree_path,
        traits_path=traits_path,
        analysis="ou-trait",
        covariance_model="ou",
        analysis_label=trait,
        taxon_column=table_audit.taxon_column,
        tree=tree,
        tree_validation=tree_validation,
        table_audit=table_audit,
        overlap=overlap,
        analysis_taxa=analysis_taxa,
        excluded_taxa=excluded_taxa,
        blockers=blockers,
        warnings=warnings,
        candidate_rows=candidate_rows,
    )


def _validated_tree_inputs(tree_path: Path):
    tree_validation = validate_tree_path(
        tree_path,
        allow_duplicates=True,
        allow_negative_branch_lengths=True,
    )
    tree = load_tree(tree_path)
    return tree_validation, tree


def _taxon_overlap(tree, table_audit) -> dict[str, list[str]]:
    tree_taxa = set(tree.tip_names)
    trait_taxa = set(table_audit.taxa)
    return {
        "matched_taxa": sorted(tree_taxa & trait_taxa),
        "missing_from_traits": sorted(tree_taxa - trait_taxa),
        "extra_trait_taxa": sorted(trait_taxa - tree_taxa),
    }


def _base_input_messages(tree_validation, table_audit) -> tuple[list[str], list[str]]:
    blockers: list[str] = []
    warnings = list(tree_validation.warnings)
    if tree_validation.duplicate_taxa:
        blockers.append("tree contains duplicate tip labels")
    if tree_validation.missing_taxa:
        blockers.append("tree contains unnamed tip labels")
    if not tree_validation.rooted:
        blockers.append("tree must be rooted for comparative covariance audit")
    if not tree_validation.has_complete_branch_lengths:
        blockers.append(
            "tree requires complete branch lengths for comparative covariance audit"
        )
    if tree_validation.negative_branch_lengths:
        blockers.append(
            "tree contains negative branch lengths that invalidate comparative covariance"
        )
    return blockers, warnings


def _trait_readiness_exclusions(
    readiness, *, trait: str
) -> list[CovarianceAuditExcludedTaxon]:
    rows: list[CovarianceAuditExcludedTaxon] = []
    rows.extend(
        CovarianceAuditExcludedTaxon(
            taxon=taxon,
            reason="missing_from_trait_table",
            details="taxon is present in the tree but absent from the trait table",
        )
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        CovarianceAuditExcludedTaxon(
            taxon=taxon,
            reason="missing_trait_value",
            details=f"trait '{trait}' is missing for this taxon",
        )
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        CovarianceAuditExcludedTaxon(
            taxon=taxon,
            reason="non_numeric_trait_value",
            details=f"trait '{trait}' is not numeric for this taxon",
        )
        for taxon in readiness.pruned_non_numeric_taxa
    )
    return rows


def _pgls_candidate_rows(
    tree, taxa: list[str], *, lambda_value: float | str
) -> list[CovarianceAuditCandidateRow]:
    base_covariance = build_brownian_covariance_matrix(tree, taxa)
    if isinstance(lambda_value, (float, int)):
        candidate_values = [float(lambda_value)]
    elif lambda_value == "estimate":
        candidate_values = [round(index / 100.0, 2) for index in range(101)]
    else:
        raise ComparativeMethodError(
            "PGLS lambda must be 'estimate' or a numeric value"
        )
    return [
        _assess_covariance_candidate(
            candidate_label=f"lambda={candidate:.2f}",
            parameter_name="lambda",
            parameter_value=candidate,
            covariance_matrix=_raw_lambda_transform_covariance(
                base_covariance, candidate
            ),
            fit_strategy="regularization",
            fit_strategy_details=(
                "PGLS applies diagonal stabilization epsilon=1e-8 before covariance inversion"
            ),
        )
        for candidate in candidate_values
    ]


def _resolved_alpha_candidates(dataset, *, alpha: float | str) -> list[float]:
    if isinstance(alpha, (float, int)):
        resolved = float(alpha)
        if resolved <= 0.0:
            raise ComparativeMethodError("OU alpha must be positive")
        return [resolved]
    if alpha != "estimate":
        raise ComparativeMethodError(
            "OU alpha must be 'estimate' or a positive numeric value"
        )
    return _alpha_grid(dataset)


def _assess_covariance_candidate(
    *,
    candidate_label: str,
    parameter_name: str | None,
    parameter_value: float | None,
    covariance_matrix: list[list[float]],
    fit_strategy: str,
    fit_strategy_details: str,
) -> CovarianceAuditCandidateRow:
    matrix_dimension = len(covariance_matrix)
    covariance_rank = matrix_rank(covariance_matrix, tolerance=1e-12)
    singular = covariance_rank < matrix_dimension
    positive_definite_before_fit = _is_positive_definite(covariance_matrix)
    condition_number = math.inf
    if positive_definite_before_fit and not singular:
        condition_number = matrix_condition_number(covariance_matrix)
    near_singular = singular or condition_number >= COVARIANCE_AUDIT_CONDITION_THRESHOLD
    fit_condition_number: float | None = None
    actual_fit_strategy = fit_strategy
    if positive_definite_before_fit and not near_singular:
        actual_fit_strategy = "exact"
        fit_strategy_details = "raw covariance is positive definite and well-conditioned enough for direct inversion without stabilization"
    elif fit_strategy == "regularization":
        stabilized = stable_covariance(
            covariance_matrix,
            epsilon=COVARIANCE_AUDIT_REGULARIZATION_EPSILON,
        )
        if _is_positive_definite(stabilized):
            fit_condition_number = matrix_condition_number(stabilized)
            fit_strategy_details = (
                fit_strategy_details
                + "; raw covariance is singular or ill-conditioned enough to require stabilization before inversion"
            )
        else:
            actual_fit_strategy = "failure"
            fit_strategy_details = (
                fit_strategy_details
                + "; stabilized covariance still cannot be inverted with a positive determinant"
            )
    return CovarianceAuditCandidateRow(
        candidate_label=candidate_label,
        parameter_name=parameter_name,
        parameter_value=parameter_value,
        matrix_dimension=matrix_dimension,
        matrix_rank=covariance_rank,
        condition_number=condition_number,
        fit_condition_number=fit_condition_number,
        positive_definite_before_fit=positive_definite_before_fit,
        singular=singular,
        near_singular=near_singular,
        fit_strategy=actual_fit_strategy,
        fit_strategy_details=fit_strategy_details,
    )


def _finalize_covariance_audit(
    *,
    tree_path: Path,
    traits_path: Path,
    analysis: str,
    covariance_model: str,
    analysis_label: str,
    taxon_column: str | None,
    tree,
    tree_validation,
    table_audit,
    overlap: dict[str, list[str]],
    analysis_taxa: list[str],
    excluded_taxa: list[CovarianceAuditExcludedTaxon],
    blockers: list[str],
    warnings: list[str],
    candidate_rows: list[CovarianceAuditCandidateRow],
) -> ComparativeCovarianceAuditReport:
    matrix_rank = min((row.matrix_rank for row in candidate_rows), default=0)
    condition_number = max(
        (row.condition_number for row in candidate_rows),
        default=math.inf if blockers else 0.0,
    )
    singular = any(row.singular for row in candidate_rows)
    near_singular = any(row.near_singular for row in candidate_rows)
    fit_strategy = _overall_fit_strategy(candidate_rows, blockers)
    fit_strategy_details = _overall_fit_strategy_details(candidate_rows, blockers)
    minimum_branch_length = None
    maximum_branch_length = None
    if tree_validation.has_complete_branch_lengths:
        branch_lengths = [
            node.branch_length
            for node in tree.iter_nodes()
            if node is not tree.root and node.branch_length is not None
        ]
        if branch_lengths:
            minimum_branch_length = min(branch_lengths)
            maximum_branch_length = max(branch_lengths)
    return ComparativeCovarianceAuditReport(
        tree_path=tree_path,
        traits_path=traits_path,
        analysis=analysis,
        covariance_model=covariance_model,
        analysis_label=analysis_label,
        taxon_column=taxon_column,
        tree_taxon_count=tree_validation.tip_count,
        trait_taxon_count=len(table_audit.taxa),
        matched_taxa=overlap["matched_taxa"],
        missing_from_traits=overlap["missing_from_traits"],
        extra_trait_taxa=overlap["extra_trait_taxa"],
        duplicate_tree_taxa=tree_validation.duplicate_taxa,
        duplicate_trait_taxa=table_audit.duplicate_taxa,
        empty_trait_taxon_rows=table_audit.empty_taxon_rows,
        analysis_taxa=analysis_taxa,
        excluded_taxa=sorted(excluded_taxa, key=lambda row: row.taxon),
        matrix_dimension=len(analysis_taxa),
        matrix_rank=matrix_rank,
        condition_number=condition_number,
        fit_strategy=fit_strategy,
        fit_strategy_details=fit_strategy_details,
        singular=singular,
        near_singular=near_singular,
        zero_length_branch_count=tree_validation.zero_length_branches,
        negative_branch_length_count=tree_validation.negative_branch_lengths,
        minimum_branch_length=minimum_branch_length,
        maximum_branch_length=maximum_branch_length,
        blockers=_deduplicate(blockers),
        warnings=_deduplicate(warnings),
        candidate_rows=candidate_rows,
    )


def _candidate_failure_messages(
    candidate_rows: list[CovarianceAuditCandidateRow],
) -> list[str]:
    if not candidate_rows:
        return []
    failed = [
        row.candidate_label for row in candidate_rows if row.fit_strategy == "failure"
    ]
    if not failed:
        return []
    return [
        "one or more covariance candidates cannot be solved after governed regularization: "
        + ", ".join(failed)
    ]


def _overall_fit_strategy(
    candidate_rows: list[CovarianceAuditCandidateRow], blockers: list[str]
) -> str:
    if blockers:
        return "failure"
    strategies = {row.fit_strategy for row in candidate_rows}
    if "failure" in strategies:
        return "failure"
    if "regularization" in strategies:
        return "regularization"
    if "pseudoinverse" in strategies:
        return "pseudoinverse"
    return "exact"


def _overall_fit_strategy_details(
    candidate_rows: list[CovarianceAuditCandidateRow], blockers: list[str]
) -> str:
    if blockers:
        return "comparative covariance audit found blockers before model fitting could proceed"
    if not candidate_rows:
        return "no covariance candidates were available for audit"
    if all(
        row.fit_strategy == candidate_rows[0].fit_strategy
        and row.fit_strategy_details == candidate_rows[0].fit_strategy_details
        for row in candidate_rows
    ):
        return candidate_rows[0].fit_strategy_details
    return "candidate covariance matrices required more than one solve path across the audited parameter surface"


def _raw_lambda_transform_covariance(
    covariance_matrix: list[list[float]], lambda_value: float
) -> list[list[float]]:
    return [
        [
            value if row_index == column_index else value * lambda_value
            for column_index, value in enumerate(row)
        ]
        for row_index, row in enumerate(covariance_matrix)
    ]


def _is_positive_definite(matrix: list[list[float]]) -> bool:
    try:
        invert_matrix(matrix)
        log_determinant(matrix)
    except ValueError:
        return False
    return True


def _deduplicate(values: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result
