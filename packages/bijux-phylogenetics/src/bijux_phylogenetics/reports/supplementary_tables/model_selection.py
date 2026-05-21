from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.artifacts.iqtree import (
    IqtreeModelCandidate,
    IqtreeModelSelectionSummary,
    parse_iqtree_model_selection_summary,
    resolve_iqtree_model_sidecar,
)

from .columns import model_selection_table_columns
from .models import (
    SupplementaryModelSelectionRow,
    SupplementaryModelSelectionTableResult,
)
from .shared import write_dict_rows


def _serialize_model_selection_row(
    row: SupplementaryModelSelectionRow,
) -> dict[str, object]:
    return {
        "iqtree_report_source": row.iqtree_report_source,
        "model_sidecar_source": ""
        if row.model_sidecar_source is None
        else row.model_sidecar_source,
        "rank": row.rank,
        "model": row.model,
        "log_likelihood": row.log_likelihood,
        "parameter_count": "" if row.parameter_count is None else row.parameter_count,
        "aic": row.aic,
        "aicc": row.aicc,
        "bic": row.bic,
        "best_aic": row.best_aic,
        "best_aicc": row.best_aicc,
        "best_bic": row.best_bic,
        "selected_model": row.selected_model,
        "selected_model_name": row.selected_model_name,
        "selected_criterion": ""
        if row.selected_criterion is None
        else row.selected_criterion,
    }


def _write_model_selection_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryModelSelectionRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_model_selection_row(row) for row in rows],
    )


def _resolve_model_sidecar_path(
    iqtree_report_path: Path,
    model_sidecar_path: Path | None,
) -> Path | None:
    if model_sidecar_path is not None:
        return model_sidecar_path
    return resolve_iqtree_model_sidecar(iqtree_report_path.with_suffix(""))


def _build_model_selection_row(
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None,
    candidate: IqtreeModelCandidate,
    summary: IqtreeModelSelectionSummary,
) -> SupplementaryModelSelectionRow:
    return SupplementaryModelSelectionRow(
        iqtree_report_source=str(iqtree_report_path),
        model_sidecar_source=(
            None if model_sidecar_path is None else str(model_sidecar_path)
        ),
        rank=candidate.rank,
        model=candidate.model,
        log_likelihood=candidate.log_likelihood,
        parameter_count=candidate.parameter_count,
        aic=candidate.aic,
        aicc=candidate.aicc,
        bic=candidate.bic,
        best_aic=candidate.model == summary.best_model_aic,
        best_aicc=candidate.model == summary.best_model_aicc,
        best_bic=candidate.model == summary.best_model_bic,
        selected_model=candidate.model == summary.selected_model,
        selected_model_name=summary.selected_model or candidate.model,
        selected_criterion=summary.selected_criterion,
    )


def write_supplementary_model_selection_table(
    path: Path,
    *,
    iqtree_report_path: Path,
    model_sidecar_path: Path | None = None,
) -> SupplementaryModelSelectionTableResult:
    """Write one supplementary model-selection table from parsed IQ-TREE artifacts."""
    resolved_sidecar_path = _resolve_model_sidecar_path(
        iqtree_report_path,
        model_sidecar_path,
    )
    summary = parse_iqtree_model_selection_summary(
        iqtree_report_path=iqtree_report_path,
        model_sidecar_path=resolved_sidecar_path,
    )
    if summary is None or summary.selected_model is None:
        raise ValueError(
            "iqtree model-selection artifacts do not expose a selected model"
        )
    if not summary.candidates:
        raise ValueError(
            "iqtree model-selection artifacts do not expose candidate model rows"
        )
    rows = [
        _build_model_selection_row(
            iqtree_report_path=iqtree_report_path,
            model_sidecar_path=resolved_sidecar_path,
            candidate=candidate,
            summary=summary,
        )
        for candidate in summary.candidates
    ]
    columns = model_selection_table_columns()
    _write_model_selection_rows(path, columns=columns, rows=rows)
    return SupplementaryModelSelectionTableResult(
        output_path=path,
        row_count=len(rows),
        selected_model=summary.selected_model,
        selected_criterion=summary.selected_criterion,
        candidate_count=summary.candidate_count,
        columns=columns,
        rows=rows,
    )
