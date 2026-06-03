from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import (
    ComparativeFormulaSpecification,
    inspect_pgls_inputs,
    run_pgls,
)
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_logistic,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_MAJOR_CLADE_FRACTION = 0.25


@dataclass(slots=True)
class ComparativeCladeStabilityRow:
    """One leave-one-clade-out model-stability summary row."""

    clade_id: str
    node_label: str | None
    dropped_taxon_count: int
    dropped_taxa: list[str]
    retained_taxon_count: int
    fit_status: str
    blocker: str | None
    baseline_term_count: int
    coefficient_comparison_count: int
    missing_baseline_term_count: int
    missing_baseline_terms: list[str]
    sign_changed_term_count: int
    significance_changed_term_count: int
    max_abs_delta_estimate: float | None
    max_abs_delta_p_value: float | None
    delta_log_likelihood: float | None
    influence_score: float | None
    rank: int


@dataclass(slots=True)
class ComparativeCladeCoefficientChangeRow:
    """One coefficient-level comparison after removing one clade."""

    clade_id: str
    node_label: str | None
    term: str
    baseline_estimate: float
    dropped_estimate: float
    delta_estimate: float
    baseline_p_value: float
    dropped_p_value: float
    delta_p_value: float
    baseline_significant: bool
    dropped_significant: bool
    sign_changed: bool
    significance_changed: bool


@dataclass(slots=True)
class ComparativeCladeStabilityReport:
    """Leave-one-clade-out stability review for one comparative model."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    model_family: str
    baseline_taxa: list[str]
    baseline_log_likelihood: float
    baseline_term_count: int
    major_clade_fraction: float
    minimum_major_clade_size: int
    candidate_clade_count: int
    blocked_clade_count: int
    influential_clades: list[str]
    clade_rows: list[ComparativeCladeStabilityRow]
    coefficient_rows: list[ComparativeCladeCoefficientChangeRow]
    warnings: list[str]


@dataclass(slots=True)
class _CoefficientSnapshot:
    estimate: float
    p_value: float


@dataclass(slots=True)
class _ComparativeFitSnapshot:
    formula: ComparativeFormulaSpecification
    model_family: str
    taxa: list[str]
    log_likelihood: float
    coefficients: dict[str, _CoefficientSnapshot]


@dataclass(slots=True)
class _CandidateClade:
    clade_id: str
    node_label: str | None
    taxa: list[str]


@dataclass(slots=True)
class _CoefficientComparison:
    rows: list[ComparativeCladeCoefficientChangeRow]
    missing_baseline_terms: list[str]


def analyze_comparative_clade_stability(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeCladeStabilityReport:
    """Refit one comparative model after removing each major clade."""
    baseline_input = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    if not baseline_input.ready:
        raise ComparativeMethodError("; ".join(baseline_input.blockers))

    baseline_fit = _fit_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, baseline_fit.taxa)
    minimum_major_clade_size = _minimum_major_clade_size(len(baseline_fit.taxa))
    candidate_clades = _collect_major_clades(
        pruned_tree.root,
        baseline_fit.taxa,
        minimum_major_clade_size=minimum_major_clade_size,
    )

    rows: list[ComparativeCladeStabilityRow] = []
    coefficient_rows: list[ComparativeCladeCoefficientChangeRow] = []
    blocked_clade_count = 0
    baseline_term_count = len(baseline_fit.coefficients)
    for candidate in candidate_clades:
        retained_taxa = [
            taxon for taxon in baseline_fit.taxa if taxon not in candidate.taxa
        ]
        try:
            dropped_fit = _refit_on_taxa(
                tree_path,
                traits_path,
                kept_taxa=retained_taxa,
                taxon_column=taxon_column,
                response=response,
                predictors=predictors,
                formula=formula,
                lambda_value=lambda_value,
            )
        except ComparativeMethodError as error:
            blocked_clade_count += 1
            rows.append(
                ComparativeCladeStabilityRow(
                    clade_id=candidate.clade_id,
                    node_label=candidate.node_label,
                    dropped_taxon_count=len(candidate.taxa),
                    dropped_taxa=list(candidate.taxa),
                    retained_taxon_count=len(retained_taxa),
                    fit_status="blocked",
                    blocker=str(error),
                    baseline_term_count=baseline_term_count,
                    coefficient_comparison_count=0,
                    missing_baseline_term_count=baseline_term_count,
                    missing_baseline_terms=sorted(baseline_fit.coefficients),
                    sign_changed_term_count=0,
                    significance_changed_term_count=0,
                    max_abs_delta_estimate=None,
                    max_abs_delta_p_value=None,
                    delta_log_likelihood=None,
                    influence_score=None,
                    rank=0,
                )
            )
            continue

        comparison = _coefficient_change_rows(
            clade_id=candidate.clade_id,
            node_label=candidate.node_label,
            baseline=baseline_fit,
            dropped=dropped_fit,
        )
        coefficient_rows.extend(comparison.rows)
        max_abs_delta_estimate = (
            max(abs(row.delta_estimate) for row in comparison.rows)
            if comparison.rows
            else 0.0
        )
        max_abs_delta_p_value = (
            max(abs(row.delta_p_value) for row in comparison.rows)
            if comparison.rows
            else 0.0
        )
        sign_changed_term_count = sum(1 for row in comparison.rows if row.sign_changed)
        significance_changed_term_count = sum(
            1 for row in comparison.rows if row.significance_changed
        )
        delta_log_likelihood = dropped_fit.log_likelihood - baseline_fit.log_likelihood
        influence_score = (
            (10.0 * sign_changed_term_count)
            + (5.0 * significance_changed_term_count)
            + (3.0 * len(comparison.missing_baseline_terms))
            + max_abs_delta_estimate
            + max_abs_delta_p_value
            + abs(delta_log_likelihood)
        )
        rows.append(
            ComparativeCladeStabilityRow(
                clade_id=candidate.clade_id,
                node_label=candidate.node_label,
                dropped_taxon_count=len(candidate.taxa),
                dropped_taxa=list(candidate.taxa),
                retained_taxon_count=len(retained_taxa),
                fit_status="fit",
                blocker=None,
                baseline_term_count=baseline_term_count,
                coefficient_comparison_count=len(comparison.rows),
                missing_baseline_term_count=len(comparison.missing_baseline_terms),
                missing_baseline_terms=comparison.missing_baseline_terms,
                sign_changed_term_count=sign_changed_term_count,
                significance_changed_term_count=significance_changed_term_count,
                max_abs_delta_estimate=max_abs_delta_estimate,
                max_abs_delta_p_value=max_abs_delta_p_value,
                delta_log_likelihood=delta_log_likelihood,
                influence_score=influence_score,
                rank=0,
            )
        )

    _rank_clade_rows(rows)
    influential_clades = [
        row.clade_id for row in rows if row.fit_status == "fit" and row.rank > 0
    ][:3]
    warnings: list[str] = []
    if blocked_clade_count:
        warnings.append(
            "one or more candidate clade removals could not be refit on the remaining taxa"
        )
    if any(row.significance_changed_term_count > 0 for row in rows):
        warnings.append(
            "one or more clade removals changed coefficient significance calls"
        )
    if any(row.sign_changed_term_count > 0 for row in rows):
        warnings.append("one or more clade removals changed coefficient direction")
    if any(row.missing_baseline_term_count > 0 for row in rows):
        warnings.append(
            "one or more clade removals changed the estimable baseline coefficient set"
        )
    return ComparativeCladeStabilityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response=baseline_fit.formula.response,
        formula=baseline_fit.formula,
        model_family=baseline_fit.model_family,
        baseline_taxa=list(baseline_fit.taxa),
        baseline_log_likelihood=baseline_fit.log_likelihood,
        baseline_term_count=baseline_term_count,
        major_clade_fraction=_MAJOR_CLADE_FRACTION,
        minimum_major_clade_size=minimum_major_clade_size,
        candidate_clade_count=len(candidate_clades),
        blocked_clade_count=blocked_clade_count,
        influential_clades=influential_clades,
        clade_rows=rows,
        coefficient_rows=coefficient_rows,
        warnings=warnings,
    )


def write_comparative_clade_stability_table(
    path: Path,
    report: ComparativeCladeStabilityReport,
) -> Path:
    """Write one leave-one-clade-out stability summary ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "node_label",
            "dropped_taxon_count",
            "dropped_taxa",
            "retained_taxon_count",
            "fit_status",
            "blocker",
            "baseline_term_count",
            "coefficient_comparison_count",
            "missing_baseline_term_count",
            "missing_baseline_terms",
            "sign_changed_term_count",
            "significance_changed_term_count",
            "max_abs_delta_estimate",
            "max_abs_delta_p_value",
            "delta_log_likelihood",
            "influence_score",
            "rank",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "node_label": row.node_label or "",
                "dropped_taxon_count": str(row.dropped_taxon_count),
                "dropped_taxa": ",".join(row.dropped_taxa),
                "retained_taxon_count": str(row.retained_taxon_count),
                "fit_status": row.fit_status,
                "blocker": row.blocker or "",
                "baseline_term_count": str(row.baseline_term_count),
                "coefficient_comparison_count": str(row.coefficient_comparison_count),
                "missing_baseline_term_count": str(row.missing_baseline_term_count),
                "missing_baseline_terms": ",".join(row.missing_baseline_terms),
                "sign_changed_term_count": str(row.sign_changed_term_count),
                "significance_changed_term_count": str(
                    row.significance_changed_term_count
                ),
                "max_abs_delta_estimate": ""
                if row.max_abs_delta_estimate is None
                else format(row.max_abs_delta_estimate, ".15g"),
                "max_abs_delta_p_value": ""
                if row.max_abs_delta_p_value is None
                else format(row.max_abs_delta_p_value, ".15g"),
                "delta_log_likelihood": ""
                if row.delta_log_likelihood is None
                else format(row.delta_log_likelihood, ".15g"),
                "influence_score": ""
                if row.influence_score is None
                else format(row.influence_score, ".15g"),
                "rank": str(row.rank),
            }
            for row in report.clade_rows
        ],
    )


def write_comparative_clade_coefficient_change_table(
    path: Path,
    report: ComparativeCladeStabilityReport,
) -> Path:
    """Write one coefficient-change ledger across leave-one-clade-out refits."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "node_label",
            "term",
            "baseline_estimate",
            "dropped_estimate",
            "delta_estimate",
            "baseline_p_value",
            "dropped_p_value",
            "delta_p_value",
            "baseline_significant",
            "dropped_significant",
            "sign_changed",
            "significance_changed",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "node_label": row.node_label or "",
                "term": row.term,
                "baseline_estimate": format(row.baseline_estimate, ".15g"),
                "dropped_estimate": format(row.dropped_estimate, ".15g"),
                "delta_estimate": format(row.delta_estimate, ".15g"),
                "baseline_p_value": format(row.baseline_p_value, ".15g"),
                "dropped_p_value": format(row.dropped_p_value, ".15g"),
                "delta_p_value": format(row.delta_p_value, ".15g"),
                "baseline_significant": str(row.baseline_significant).lower(),
                "dropped_significant": str(row.dropped_significant).lower(),
                "sign_changed": str(row.sign_changed).lower(),
                "significance_changed": str(row.significance_changed).lower(),
            }
            for row in report.coefficient_rows
        ],
    )


def _fit_snapshot(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
    taxon_column: str | None,
    lambda_value: float | str,
) -> _ComparativeFitSnapshot:
    input_report = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    if not input_report.ready:
        raise ComparativeMethodError("; ".join(input_report.blockers))
    response_values = [row.response_value for row in input_report.model_matrix.rows]
    model_family = _shared_response_family(response_values)
    if model_family == "logistic":
        if lambda_value == "estimate":
            raise ComparativeMethodError(
                "leave-one-clade-out stability requires a numeric lambda value for binary-response comparative models"
            )
        report = summarize_phylogenetic_logistic(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
            lambda_value=float(lambda_value),
        )
        return _ComparativeFitSnapshot(
            formula=report.formula,
            model_family="logistic",
            taxa=[row.taxon for row in report.fitted_rows],
            log_likelihood=report.binomial_log_likelihood,
            coefficients={
                row.name: _CoefficientSnapshot(
                    estimate=row.estimate,
                    p_value=row.p_value,
                )
                for row in report.coefficients
            },
        )
    report = run_pgls(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    return _ComparativeFitSnapshot(
        formula=report.formula,
        model_family="pgls",
        taxa=list(report.taxa),
        log_likelihood=report.log_likelihood,
        coefficients={
            row.name: _CoefficientSnapshot(
                estimate=row.estimate,
                p_value=row.p_value,
            )
            for row in report.coefficients
        },
    )


def _shared_response_family(response_values: list[float]) -> str:
    if all(
        math.isclose(value, round(value), abs_tol=1e-12) for value in response_values
    ) and {int(round(value)) for value in response_values} <= {0, 1}:
        return "logistic"
    return "pgls"


def _minimum_major_clade_size(total_taxa: int) -> int:
    return max(2, math.ceil(total_taxa * _MAJOR_CLADE_FRACTION))


def _collect_major_clades(
    root: TreeNode,
    baseline_taxa: list[str],
    *,
    minimum_major_clade_size: int,
) -> list[_CandidateClade]:
    total_taxa = len(baseline_taxa)
    rows: list[_CandidateClade] = []

    def visit(node: TreeNode, *, is_root: bool) -> list[str]:
        if node.is_leaf():
            return [node.name] if node.name is not None else []
        taxa: list[str] = []
        for child in node.children:
            taxa.extend(visit(child, is_root=False))
        ordered_taxa = sorted(taxa)
        if (
            not is_root
            and len(ordered_taxa) >= minimum_major_clade_size
            and len(ordered_taxa) <= total_taxa - 2
        ):
            rows.append(
                _CandidateClade(
                    clade_id="|".join(ordered_taxa),
                    node_label=node.name,
                    taxa=ordered_taxa,
                )
            )
        return ordered_taxa

    visit(root, is_root=True)
    return rows


def _refit_on_taxa(
    tree_path: Path,
    traits_path: Path,
    *,
    kept_taxa: list[str],
    taxon_column: str | None,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
    lambda_value: float | str,
) -> _ComparativeFitSnapshot:
    if len(kept_taxa) < 3:
        raise ComparativeMethodError(
            "remaining taxa are insufficient for comparative refitting after clade removal"
        )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    missing_taxa = [taxon for taxon in kept_taxa if taxon not in rows_by_taxon]
    if missing_taxa:
        raise ComparativeMethodError(
            "trait table is missing one or more retained taxa after clade removal"
        )
    rows = [rows_by_taxon[taxon] for taxon in kept_taxa]
    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-clade-stability-"
    ) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        reduced_tree_path = tmp_dir_path / "comparative-clade-stability-tree.nwk"
        reduced_table_path = tmp_dir_path / "comparative-clade-stability-traits.tsv"
        reduced_tree_path.write_text(dumps_newick(pruned_tree) + "\n", encoding="utf-8")
        write_taxon_rows(
            reduced_table_path,
            columns=table.columns,
            rows=rows,
        )
        return _fit_snapshot(
            reduced_tree_path,
            reduced_table_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=table.taxon_column,
            lambda_value=lambda_value,
        )


def _coefficient_change_rows(
    *,
    clade_id: str,
    node_label: str | None,
    baseline: _ComparativeFitSnapshot,
    dropped: _ComparativeFitSnapshot,
) -> _CoefficientComparison:
    shared_terms = sorted(set(baseline.coefficients) & set(dropped.coefficients))
    missing_baseline_terms = sorted(
        set(baseline.coefficients) - set(dropped.coefficients)
    )
    rows: list[ComparativeCladeCoefficientChangeRow] = []
    for term in shared_terms:
        baseline_term = baseline.coefficients[term]
        dropped_term = dropped.coefficients[term]
        rows.append(
            ComparativeCladeCoefficientChangeRow(
                clade_id=clade_id,
                node_label=node_label,
                term=term,
                baseline_estimate=baseline_term.estimate,
                dropped_estimate=dropped_term.estimate,
                delta_estimate=dropped_term.estimate - baseline_term.estimate,
                baseline_p_value=baseline_term.p_value,
                dropped_p_value=dropped_term.p_value,
                delta_p_value=dropped_term.p_value - baseline_term.p_value,
                baseline_significant=baseline_term.p_value <= 0.05,
                dropped_significant=dropped_term.p_value <= 0.05,
                sign_changed=(
                    _estimate_sign(baseline_term.estimate)
                    != _estimate_sign(dropped_term.estimate)
                ),
                significance_changed=(
                    (baseline_term.p_value <= 0.05) != (dropped_term.p_value <= 0.05)
                ),
            )
        )
    return _CoefficientComparison(
        rows=rows,
        missing_baseline_terms=missing_baseline_terms,
    )


def _estimate_sign(value: float) -> int:
    if math.isclose(value, 0.0, abs_tol=1e-12):
        return 0
    return 1 if value > 0.0 else -1


def _rank_clade_rows(rows: list[ComparativeCladeStabilityRow]) -> None:
    ranked_rows = sorted(
        [
            row
            for row in rows
            if row.fit_status == "fit" and row.influence_score is not None
        ],
        key=lambda row: (
            -(row.influence_score or 0.0),
            -(abs(row.delta_log_likelihood or 0.0)),
            -(row.significance_changed_term_count),
            -(row.sign_changed_term_count),
            row.clade_id,
        ),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row.rank = rank
