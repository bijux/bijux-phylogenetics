from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.clades import robinson_foulds_metrics

from .balanced_minimum_evolution import score_balanced_minimum_evolution
from .imported import (
    build_tree_from_imported_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import (
    DistanceMethodComparisonReport,
    DistanceMethodComparisonRow,
    DistanceMethodRfRow,
    DistanceMethodWarningRow,
    MissingDistancePolicy,
)
from .ordinary_least_squares import fit_ordinary_least_squares_tree
from .patristic_residuals import compute_patristic_residual_diagnostics
from .taxon_leave_one_out import (
    ROUND_DIGITS,
    defined_distance_lookup_from_entries,
    resolve_distance_lookup,
)

_COMPARED_METHODS = [
    "neighbor-joining",
    "bionj",
    "upgma",
    "wpgma",
]


def _method_warning_rows(
    *,
    ultrametric_compatible: bool,
    rows: list[DistanceMethodComparisonRow],
    matrix_warnings: list[str],
) -> list[DistanceMethodWarningRow]:
    warnings: list[DistanceMethodWarningRow] = []
    for warning in matrix_warnings:
        warnings.append(
            DistanceMethodWarningRow(
                warning_rank=0,
                scope="matrix",
                method=None,
                warning=warning,
            )
        )
    for row in rows:
        for warning in row.assumption_warnings:
            warnings.append(
                DistanceMethodWarningRow(
                    warning_rank=0,
                    scope="method",
                    method=row.method,
                    warning=warning,
                )
            )
    if not ultrametric_compatible:
        for method in ("upgma", "wpgma"):
            warnings.append(
                DistanceMethodWarningRow(
                    warning_rank=0,
                    scope="method",
                    method=method,
                    warning=(
                        "source matrix is not ultrametric, so this clock-like clustering method is assumption-violating on this input"
                    ),
                )
            )
    return [
        DistanceMethodWarningRow(
            warning_rank=index,
            scope=row.scope,
            method=row.method,
            warning=row.warning,
        )
        for index, row in enumerate(
            sorted(
                warnings,
                key=lambda row: (
                    row.scope,
                    "" if row.method is None else row.method,
                    row.warning,
                ),
            ),
            start=1,
        )
    ]


def _rf_rows(
    method_to_tree: dict[str, object],
    taxa: list[str],
) -> list[DistanceMethodRfRow]:
    shared_taxa = set(taxa)
    rows: list[DistanceMethodRfRow] = []
    for left_index, left_method in enumerate(_COMPARED_METHODS):
        left_tree = method_to_tree[left_method]
        for right_method in _COMPARED_METHODS[left_index + 1 :]:
            right_tree = method_to_tree[right_method]
            metrics = robinson_foulds_metrics(
                left_tree,
                right_tree,
                shared_taxa,
                rf_mode="rooted",
            )
            rows.append(
                DistanceMethodRfRow(
                    left_method=left_method,
                    right_method=right_method,
                    rooted_robinson_foulds_distance=metrics.distance,
                    rooted_normalized_robinson_foulds=round(
                        metrics.normalized_distance,
                        ROUND_DIGITS,
                    ),
                )
            )
    return rows


def compare_distance_tree_methods_from_imported_distance_matrix(
    matrix_path: Path,
    *,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> DistanceMethodComparisonReport:
    """Build owned distance trees from one imported matrix and compare them under shared scoring surfaces."""
    validation = validate_imported_distance_matrix(matrix_path)
    entries = load_imported_distance_matrix(matrix_path)
    defined_lookup = defined_distance_lookup_from_entries(entries)
    resolved_lookup = resolve_distance_lookup(
        validation.identifiers,
        defined_lookup,
        missing_distance_policy=missing_distance_policy,
    )

    method_to_tree = {}
    rows: list[DistanceMethodComparisonRow] = []
    matrix_warnings: list[str] = []
    ultrametric_compatible = True
    for method in _COMPARED_METHODS:
        tree, build_report = build_tree_from_imported_distance_matrix(
            matrix_path,
            method=method,
            missing_distance_policy=missing_distance_policy,
        )
        method_to_tree[method] = tree
        if not matrix_warnings:
            matrix_warnings = list(build_report.assumptions.warnings)
            ultrametric_compatible = build_report.assumptions.ultrametric_compatible
        patristic_residuals = compute_patristic_residual_diagnostics(
            tree,
            validation.identifiers,
            resolved_lookup,
        )
        _ols_tree, ols_report = fit_ordinary_least_squares_tree(
            tree,
            validation.identifiers,
            resolved_lookup,
        )
        rows.append(
            DistanceMethodComparisonRow(
                method=method,
                tree_newick=dumps_newick(tree),
                patristic_residual_sum_squares=patristic_residuals.residual_sum_squares,
                balanced_minimum_evolution_score=score_balanced_minimum_evolution(
                    tree,
                    validation.identifiers,
                    resolved_lookup,
                ),
                ordinary_least_squares_residual_sum_squares=(
                    ols_report.residual_sum_squares
                ),
                ordinary_least_squares_negative_branch_count=(
                    ols_report.negative_branch_count
                ),
                assumption_warnings=list(build_report.method_policy.limitations),
            )
        )

    return DistanceMethodComparisonReport(
        source_path=matrix_path,
        source_kind="imported-distance-matrix",
        missing_distance_policy=missing_distance_policy,
        taxa=list(validation.identifiers),
        compared_methods=list(_COMPARED_METHODS),
        rows=rows,
        rf_rows=_rf_rows(method_to_tree, validation.identifiers),
        warning_rows=_method_warning_rows(
            ultrametric_compatible=ultrametric_compatible,
            rows=rows,
            matrix_warnings=matrix_warnings,
        ),
    )


def write_distance_method_comparison_scores(
    path: Path,
    report: DistanceMethodComparisonReport,
) -> Path:
    """Write the per-method score and residual ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "method",
        "patristic_residual_sum_squares",
        "balanced_minimum_evolution_score",
        "ordinary_least_squares_residual_sum_squares",
        "ordinary_least_squares_negative_branch_count",
        "tree_newick",
    ]
    lines = ["\t".join(columns)]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    row.method,
                    format(row.patristic_residual_sum_squares, ".12g"),
                    format(row.balanced_minimum_evolution_score, ".12g"),
                    format(row.ordinary_least_squares_residual_sum_squares, ".12g"),
                    str(row.ordinary_least_squares_negative_branch_count),
                    row.tree_newick,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_method_comparison_rf_matrix(
    path: Path,
    report: DistanceMethodComparisonReport,
) -> Path:
    """Write the rooted RF distance matrix across all compared distance-tree methods."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rf_lookup = {
        (row.left_method, row.right_method): row.rooted_robinson_foulds_distance
        for row in report.rf_rows
    }
    methods = report.compared_methods
    lines = ["method\t" + "\t".join(methods)]
    for left_method in methods:
        values = [left_method]
        for right_method in methods:
            if left_method == right_method:
                values.append("0")
                continue
            pair = (
                (left_method, right_method)
                if (left_method, right_method) in rf_lookup
                else (right_method, left_method)
            )
            values.append(str(rf_lookup[pair]))
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_method_comparison_warnings(
    path: Path,
    report: DistanceMethodComparisonReport,
) -> Path:
    """Write the explicit matrix-level and method-level assumption warnings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["warning_rank", "scope", "method", "warning"]
    lines = ["\t".join(columns)]
    for row in report.warning_rows:
        lines.append(
            "\t".join(
                [
                    str(row.warning_rank),
                    row.scope,
                    "" if row.method is None else row.method,
                    row.warning,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_method_comparison_run_json(
    path: Path,
    report: DistanceMethodComparisonReport,
) -> Path:
    """Write the full method-comparison payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source_path": str(report.source_path),
        "source_kind": report.source_kind,
        "missing_distance_policy": report.missing_distance_policy,
        "taxa": report.taxa,
        "compared_methods": report.compared_methods,
        "rows": [asdict(row) for row in report.rows],
        "rf_rows": [asdict(row) for row in report.rf_rows],
        "warning_rows": [asdict(row) for row in report.warning_rows],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_distance_method_comparison_artifacts(
    out_dir: Path,
    report: DistanceMethodComparisonReport,
) -> dict[str, Path]:
    """Write governed distance-method comparison artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_paths: dict[str, Path] = {}
    for row in report.rows:
        tree_path = out_dir / f"{row.method}.nwk"
        tree_path.write_text(row.tree_newick + "\n", encoding="utf-8")
        tree_paths[row.method] = tree_path
    score_table_path = write_distance_method_comparison_scores(
        out_dir / "method_scores.tsv",
        report,
    )
    rf_matrix_path = write_distance_method_comparison_rf_matrix(
        out_dir / "rf_matrix.tsv",
        report,
    )
    warnings_path = write_distance_method_comparison_warnings(
        out_dir / "assumption_warnings.tsv",
        report,
    )
    run_json_path = write_distance_method_comparison_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        **{f"{method}_tree": path for method, path in tree_paths.items()},
        "method_scores": score_table_path,
        "rf_matrix": rf_matrix_path,
        "assumption_warnings": warnings_path,
        "run_json": run_json_path,
    }
