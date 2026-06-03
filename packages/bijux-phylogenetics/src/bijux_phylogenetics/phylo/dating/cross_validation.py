from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .calibrations import load_fixed_dating_calibrations
from .inputs import load_tip_dates_for_tree
from .models import (
    DatingCalibrationAnchor,
    PenalizedLikelihoodCrossValidationCandidateRow,
    PenalizedLikelihoodCrossValidationPredictionRow,
    PenalizedLikelihoodCrossValidationReport,
)
from .penalized_likelihood import (
    fit_penalized_likelihood_dating,
    write_penalized_likelihood_branch_rate_tsv,
    write_penalized_likelihood_node_dates_tsv,
)


def cross_validate_penalized_likelihood_smoothing(
    tree,
    tip_dates: dict[str, float],
    calibration_rows: Sequence[DatingCalibrationAnchor],
    *,
    smoothing_parameters: Sequence[float],
    max_coordinate_passes: int = 8,
    tree_path: Path | None = None,
    metadata_path: Path | None = None,
    calibration_path: Path | None = None,
    taxon_column: str = "taxon",
    date_column: str = "date",
) -> PenalizedLikelihoodCrossValidationReport:
    """Select one penalized-dating smoothing value by held-out calibration prediction error."""
    candidate_parameters = _normalize_smoothing_parameters(smoothing_parameters)
    if len(calibration_rows) < 2:
        raise PhylogeneticsError(
            "smoothing cross-validation requires at least two fixed internal-node calibrations",
            code="penalized_likelihood_cross_validation_error",
        )

    prediction_rows: list[PenalizedLikelihoodCrossValidationPredictionRow] = []
    candidate_rows: list[PenalizedLikelihoodCrossValidationCandidateRow] = []
    for smoothing_parameter in candidate_parameters:
        candidate_predictions: list[
            PenalizedLikelihoodCrossValidationPredictionRow
        ] = []
        for held_out in calibration_rows:
            training_node_dates = {
                row.node_id: row.fixed_date
                for row in calibration_rows
                if row.calibration_id != held_out.calibration_id
            }
            held_out_fit = fit_penalized_likelihood_dating(
                tree,
                tip_dates,
                fixed_node_dates=training_node_dates,
                smoothing_parameter=smoothing_parameter,
                max_coordinate_passes=max_coordinate_passes,
                tree_path=tree_path,
                metadata_path=metadata_path,
                taxon_column=taxon_column,
                date_column=date_column,
            )
            predicted_date = _node_date_from_selected_fit(
                held_out_fit,
                node_id=held_out.node_id,
            )
            absolute_error = abs(predicted_date - held_out.fixed_date)
            candidate_predictions.append(
                PenalizedLikelihoodCrossValidationPredictionRow(
                    smoothing_parameter=smoothing_parameter,
                    held_out_calibration_id=held_out.calibration_id,
                    held_out_target_label=held_out.target_label,
                    held_out_descendant_taxa=held_out.descendant_taxa,
                    held_out_node_id=held_out.node_id,
                    training_calibration_count=len(calibration_rows) - 1,
                    observed_date=held_out.fixed_date,
                    predicted_date=predicted_date,
                    absolute_error=absolute_error,
                    squared_error=absolute_error * absolute_error,
                    optimization_pass_count=held_out_fit.optimization_pass_count,
                    function_evaluation_count=held_out_fit.function_evaluation_count,
                    converged=held_out_fit.converged,
                )
            )
        prediction_rows.extend(candidate_predictions)
        candidate_rows.append(
            _build_candidate_row(
                smoothing_parameter=smoothing_parameter,
                prediction_rows=candidate_predictions,
            )
        )

    selected_smoothing_parameter = _select_smoothing_parameter(candidate_rows)
    selected_candidate_row = next(
        row
        for row in candidate_rows
        if math.isclose(
            row.smoothing_parameter,
            selected_smoothing_parameter,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    candidate_rows = [
        PenalizedLikelihoodCrossValidationCandidateRow(
            smoothing_parameter=row.smoothing_parameter,
            fold_count=row.fold_count,
            mean_absolute_error=row.mean_absolute_error,
            mean_squared_error=row.mean_squared_error,
            root_mean_squared_error=row.root_mean_squared_error,
            max_absolute_error=row.max_absolute_error,
            selected=math.isclose(
                row.smoothing_parameter,
                selected_smoothing_parameter,
                rel_tol=0.0,
                abs_tol=1e-12,
            ),
        )
        for row in candidate_rows
    ]
    selected_fit = fit_penalized_likelihood_dating(
        tree,
        tip_dates,
        fixed_node_dates={row.node_id: row.fixed_date for row in calibration_rows},
        smoothing_parameter=selected_smoothing_parameter,
        max_coordinate_passes=max_coordinate_passes,
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    return PenalizedLikelihoodCrossValidationReport(
        tree_newick=dumps_newick(tree),
        taxa=sorted(tree.tip_names),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        branch_count=sum(1 for _parent, _child in tree.iter_edges()),
        tree_path=None if tree_path is None else str(tree_path),
        metadata_path=None if metadata_path is None else str(metadata_path),
        calibration_path=None if calibration_path is None else str(calibration_path),
        taxon_column=taxon_column,
        date_column=date_column,
        usable_calibration_count=len(calibration_rows),
        candidate_count=len(candidate_rows),
        selected_smoothing_parameter=selected_smoothing_parameter,
        selected_mean_absolute_error=selected_candidate_row.mean_absolute_error,
        selected_mean_squared_error=selected_candidate_row.mean_squared_error,
        selected_root_mean_squared_error=selected_candidate_row.root_mean_squared_error,
        calibration_rows=list(calibration_rows),
        candidate_rows=candidate_rows,
        prediction_rows=prediction_rows,
        selected_fit=selected_fit,
    )


def cross_validate_penalized_likelihood_smoothing_from_metadata(
    tree_path: Path,
    metadata_path: Path,
    calibration_path: Path,
    *,
    smoothing_parameters: Sequence[float],
    max_coordinate_passes: int = 8,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> PenalizedLikelihoodCrossValidationReport:
    """Run held-out calibration cross-validation from one tree, tip-date table, and calibration table."""
    validate_tree_path(tree_path, require_rooted=True)
    tree = load_tree(tree_path)
    tree.rooted = True
    tip_dates, resolved_taxon_column = load_tip_dates_for_tree(
        metadata_path,
        tree_taxa=tree.tip_names,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    calibration_rows = load_fixed_dating_calibrations(tree_path, calibration_path)
    return cross_validate_penalized_likelihood_smoothing(
        tree,
        tip_dates,
        calibration_rows,
        smoothing_parameters=smoothing_parameters,
        max_coordinate_passes=max_coordinate_passes,
        tree_path=tree_path,
        metadata_path=metadata_path,
        calibration_path=calibration_path,
        taxon_column=resolved_taxon_column,
        date_column=date_column,
    )


def _normalize_smoothing_parameters(
    smoothing_parameters: Sequence[float],
) -> list[float]:
    normalized: list[float] = []
    for raw_value in smoothing_parameters:
        try:
            smoothing_parameter = float(raw_value)
        except (TypeError, ValueError) as error:
            raise PhylogeneticsError(
                f"invalid smoothing parameter '{raw_value}'",
                code="penalized_likelihood_cross_validation_error",
            ) from error
        if smoothing_parameter <= 0.0:
            raise PhylogeneticsError(
                "smoothing cross-validation requires strictly positive smoothing parameters",
                code="penalized_likelihood_cross_validation_error",
            )
        if not any(
            math.isclose(
                smoothing_parameter,
                existing,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            for existing in normalized
        ):
            normalized.append(smoothing_parameter)
    if not normalized:
        raise PhylogeneticsError(
            "smoothing cross-validation requires at least one smoothing parameter candidate",
            code="penalized_likelihood_cross_validation_error",
        )
    return sorted(normalized)


def _build_candidate_row(
    *,
    smoothing_parameter: float,
    prediction_rows: Sequence[PenalizedLikelihoodCrossValidationPredictionRow],
) -> PenalizedLikelihoodCrossValidationCandidateRow:
    absolute_errors = [row.absolute_error for row in prediction_rows]
    squared_errors = [row.squared_error for row in prediction_rows]
    mean_absolute_error = sum(absolute_errors) / len(absolute_errors)
    mean_squared_error = sum(squared_errors) / len(squared_errors)
    return PenalizedLikelihoodCrossValidationCandidateRow(
        smoothing_parameter=smoothing_parameter,
        fold_count=len(prediction_rows),
        mean_absolute_error=mean_absolute_error,
        mean_squared_error=mean_squared_error,
        root_mean_squared_error=math.sqrt(mean_squared_error),
        max_absolute_error=max(absolute_errors),
        selected=False,
    )


def _select_smoothing_parameter(
    candidate_rows: Sequence[PenalizedLikelihoodCrossValidationCandidateRow],
) -> float:
    return min(
        candidate_rows,
        key=lambda row: (
            row.root_mean_squared_error,
            row.mean_absolute_error,
            row.smoothing_parameter,
        ),
    ).smoothing_parameter


def _node_date_from_selected_fit(selected_fit, *, node_id: str) -> float:
    for row in selected_fit.node_rows:
        if row.node_id == node_id:
            return row.estimated_date
    raise PhylogeneticsError(
        f"selected penalized dating fit is missing node '{node_id}'",
        code="penalized_likelihood_cross_validation_error",
    )


def write_penalized_likelihood_cross_validation_summary_tsv(
    path: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> Path:
    """Write one summary row for one smoothing cross-validation run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tree_path",
        "metadata_path",
        "calibration_path",
        "taxon_column",
        "date_column",
        "tip_count",
        "internal_node_count",
        "branch_count",
        "usable_calibration_count",
        "candidate_count",
        "selected_smoothing_parameter",
        "selected_mean_absolute_error",
        "selected_mean_squared_error",
        "selected_root_mean_squared_error",
        "final_root_date",
        "final_data_score",
        "final_penalty_score",
        "final_total_score",
        "final_converged",
    ]
    values = [
        report.tree_path or "",
        report.metadata_path or "",
        report.calibration_path or "",
        report.taxon_column,
        report.date_column,
        str(report.tip_count),
        str(report.internal_node_count),
        str(report.branch_count),
        str(report.usable_calibration_count),
        str(report.candidate_count),
        format(report.selected_smoothing_parameter, ".15g"),
        format(report.selected_mean_absolute_error, ".15g"),
        format(report.selected_mean_squared_error, ".15g"),
        format(report.selected_root_mean_squared_error, ".15g"),
        format(report.selected_fit.root_date, ".15g"),
        format(report.selected_fit.data_score, ".15g"),
        format(report.selected_fit.penalty_score, ".15g"),
        format(report.selected_fit.total_score, ".15g"),
        str(report.selected_fit.converged).lower(),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_penalized_likelihood_cross_validation_candidates_tsv(
    path: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> Path:
    """Write one candidate error row per evaluated smoothing parameter."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "smoothing_parameter",
        "fold_count",
        "mean_absolute_error",
        "mean_squared_error",
        "root_mean_squared_error",
        "max_absolute_error",
        "selected",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                format(row.smoothing_parameter, ".15g"),
                str(row.fold_count),
                format(row.mean_absolute_error, ".15g"),
                format(row.mean_squared_error, ".15g"),
                format(row.root_mean_squared_error, ".15g"),
                format(row.max_absolute_error, ".15g"),
                str(row.selected).lower(),
            ]
        )
        for row in report.candidate_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_penalized_likelihood_cross_validation_predictions_tsv(
    path: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> Path:
    """Write one held-out prediction row per smoothing parameter and calibration."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "smoothing_parameter",
        "held_out_calibration_id",
        "held_out_target_label",
        "held_out_descendant_taxa",
        "held_out_node_id",
        "training_calibration_count",
        "observed_date",
        "predicted_date",
        "absolute_error",
        "squared_error",
        "optimization_pass_count",
        "function_evaluation_count",
        "converged",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                format(row.smoothing_parameter, ".15g"),
                row.held_out_calibration_id,
                row.held_out_target_label,
                "|".join(row.held_out_descendant_taxa),
                row.held_out_node_id,
                str(row.training_calibration_count),
                format(row.observed_date, ".15g"),
                format(row.predicted_date, ".15g"),
                format(row.absolute_error, ".15g"),
                format(row.squared_error, ".15g"),
                str(row.optimization_pass_count),
                str(row.function_evaluation_count),
                str(row.converged).lower(),
            ]
        )
        for row in report.prediction_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_penalized_likelihood_cross_validation_calibrations_tsv(
    path: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> Path:
    """Write the resolved fixed calibrations used for cross-validation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "calibration_id",
        "target_kind",
        "target_label",
        "descendant_taxa",
        "node_id",
        "node_kind",
        "fixed_date",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.calibration_id,
                row.target_kind,
                row.target_label,
                "|".join(row.descendant_taxa),
                row.node_id,
                row.node_kind,
                format(row.fixed_date, ".15g"),
            ]
        )
        for row in report.calibration_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_penalized_likelihood_cross_validation_run_json(
    path: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> Path:
    """Write the full smoothing cross-validation report as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_penalized_likelihood_cross_validation_artifacts(
    out_dir: Path,
    report: PenalizedLikelihoodCrossValidationReport,
) -> dict[str, Path]:
    """Write governed artifact outputs for one smoothing cross-validation run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    dated_tree_path = write_newick(
        out_dir / "dated_tree.nwk",
        loads_newick(report.selected_fit.dated_tree_newick),
    )
    summary_path = write_penalized_likelihood_cross_validation_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    candidate_scores_path = write_penalized_likelihood_cross_validation_candidates_tsv(
        out_dir / "candidate_scores.tsv",
        report,
    )
    prediction_errors_path = (
        write_penalized_likelihood_cross_validation_predictions_tsv(
            out_dir / "prediction_errors.tsv",
            report,
        )
    )
    calibrations_path = write_penalized_likelihood_cross_validation_calibrations_tsv(
        out_dir / "calibrations.tsv",
        report,
    )
    node_dates_path = write_penalized_likelihood_node_dates_tsv(
        out_dir / "node_dates.tsv",
        report.selected_fit,
    )
    branch_rates_path = write_penalized_likelihood_branch_rate_tsv(
        out_dir / "branch_rates.tsv",
        report.selected_fit,
    )
    run_json_path = write_penalized_likelihood_cross_validation_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "dated_tree_path": dated_tree_path,
        "summary_path": summary_path,
        "candidate_scores_path": candidate_scores_path,
        "prediction_errors_path": prediction_errors_path,
        "calibrations_path": calibrations_path,
        "node_dates_path": node_dates_path,
        "branch_rates_path": branch_rates_path,
        "run_json_path": run_json_path,
    }
