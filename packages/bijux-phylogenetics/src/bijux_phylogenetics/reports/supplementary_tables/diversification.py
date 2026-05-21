from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.diversification import (
    CladeDiversificationObservation,
    DiversificationRateReport,
    SamplingFractionIssue,
    SamplingFractionReport,
    compare_diversification_models,
    detect_diversification_outlier_clades,
    detect_incomplete_taxon_sampling_metadata,
    estimate_diversification_rate,
)

from .columns import diversification_table_columns
from .models import (
    SupplementaryDiversificationRow,
    SupplementaryDiversificationTableResult,
)
from .shared import stringify_list, write_dict_rows


def _serialize_sampling_issue(row: SamplingFractionIssue) -> str:
    raw_value = row.raw_value if row.raw_value else "<missing>"
    return f"{row.taxon}:{row.code}:{raw_value}"


def _serialize_diversification_row(
    row: SupplementaryDiversificationRow,
) -> dict[str, str]:
    return {
        "tree_source": row.tree_source,
        "metadata_source": "" if row.metadata_source is None else row.metadata_source,
        "clade_model": row.clade_model,
        "better_model": row.better_model,
        "node": row.node,
        "node_name": "" if row.node_name is None else row.node_name,
        "descendant_taxa": stringify_list(row.descendant_taxa),
        "tip_count": str(row.tip_count),
        "crown_age": str(row.crown_age),
        "clade_diversification_rate": str(row.clade_diversification_rate),
        "clade_rate_z_score": str(row.clade_rate_z_score),
        "clade_classification": row.clade_classification,
        "global_diversification_rate": str(row.global_diversification_rate),
        "yule_log_likelihood": str(row.yule_log_likelihood),
        "yule_aic": str(row.yule_aic),
        "yule_corrected_tip_count": str(row.yule_corrected_tip_count),
        "yule_sampling_fraction": str(row.yule_sampling_fraction),
        "yule_net_diversification_rate": str(row.yule_net_diversification_rate),
        "yule_relative_extinction": str(row.yule_relative_extinction),
        "birth_death_log_likelihood": str(row.birth_death_log_likelihood),
        "birth_death_aic": str(row.birth_death_aic),
        "birth_death_corrected_tip_count": str(row.birth_death_corrected_tip_count),
        "birth_death_sampling_fraction": str(row.birth_death_sampling_fraction),
        "birth_death_net_diversification_rate": str(
            row.birth_death_net_diversification_rate
        ),
        "birth_death_relative_extinction": str(row.birth_death_relative_extinction),
        "sampling_metadata_complete": (
            ""
            if row.sampling_metadata_complete is None
            else str(row.sampling_metadata_complete).lower()
        ),
        "sampling_column": "" if row.sampling_column is None else row.sampling_column,
        "sampling_fraction": (
            "" if row.sampling_fraction is None else str(row.sampling_fraction)
        ),
        "sampling_heterogeneous": (
            ""
            if row.sampling_heterogeneous is None
            else str(row.sampling_heterogeneous).lower()
        ),
        "sampling_missing_taxa": stringify_list(row.sampling_missing_taxa),
        "sampling_invalid_rows": stringify_list(row.sampling_invalid_rows),
        "warning_count": str(row.warning_count),
        "warnings": stringify_list(row.warnings),
    }


def _write_diversification_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryDiversificationRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_diversification_row(row) for row in rows],
    )


def _diversification_model_by_name(
    reports: list[DiversificationRateReport],
) -> dict[str, DiversificationRateReport]:
    return {report.model: report for report in reports}


def _build_diversification_row(
    *,
    observation: CladeDiversificationObservation,
    tree_path: Path,
    metadata_path: Path | None,
    clade_model: str,
    better_model: str,
    global_diversification_rate: float,
    yule_report: DiversificationRateReport,
    birth_death_report: DiversificationRateReport,
    sampling_report: SamplingFractionReport | None,
    warnings: list[str],
) -> SupplementaryDiversificationRow:
    return SupplementaryDiversificationRow(
        tree_source=str(tree_path),
        metadata_source=None if metadata_path is None else str(metadata_path),
        clade_model=clade_model,
        better_model=better_model,
        node=observation.node,
        node_name=observation.node_name,
        descendant_taxa=list(observation.descendant_taxa),
        tip_count=observation.tip_count,
        crown_age=observation.crown_age,
        clade_diversification_rate=observation.diversification_rate,
        clade_rate_z_score=observation.z_score,
        clade_classification=observation.classification,
        global_diversification_rate=global_diversification_rate,
        yule_log_likelihood=yule_report.log_likelihood,
        yule_aic=yule_report.aic,
        yule_corrected_tip_count=yule_report.corrected_tip_count,
        yule_sampling_fraction=yule_report.sampling_fraction,
        yule_net_diversification_rate=yule_report.net_diversification_rate,
        yule_relative_extinction=yule_report.relative_extinction,
        birth_death_log_likelihood=birth_death_report.log_likelihood,
        birth_death_aic=birth_death_report.aic,
        birth_death_corrected_tip_count=birth_death_report.corrected_tip_count,
        birth_death_sampling_fraction=birth_death_report.sampling_fraction,
        birth_death_net_diversification_rate=(
            birth_death_report.net_diversification_rate
        ),
        birth_death_relative_extinction=birth_death_report.relative_extinction,
        sampling_metadata_complete=(
            None if sampling_report is None else sampling_report.complete
        ),
        sampling_column=None
        if sampling_report is None
        else sampling_report.sampling_column,
        sampling_fraction=(
            None if sampling_report is None else sampling_report.sampling_fraction
        ),
        sampling_heterogeneous=(
            None if sampling_report is None else sampling_report.heterogeneous_values
        ),
        sampling_missing_taxa=(
            [] if sampling_report is None else list(sampling_report.missing_taxa)
        ),
        sampling_invalid_rows=(
            []
            if sampling_report is None
            else [
                _serialize_sampling_issue(issue)
                for issue in sampling_report.invalid_rows
            ]
        ),
        warning_count=len(warnings),
        warnings=list(warnings),
    )


def write_supplementary_diversification_table(
    path: Path,
    *,
    tree_path: Path,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    clade_model: str = "birth-death",
) -> SupplementaryDiversificationTableResult:
    """Write one supplementary diversification table over clade and model evidence."""
    if clade_model not in {"yule", "birth-death"}:
        raise ValueError("clade_model must be 'yule' or 'birth-death'")
    sampling_report = (
        None
        if metadata_path is None
        else detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
    )
    comparison = compare_diversification_models(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    diversification_reports = [
        estimate_diversification_rate(
            tree_path,
            metadata_path=metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
            model="yule",
        ),
        estimate_diversification_rate(
            tree_path,
            metadata_path=metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
            model="birth-death",
        ),
    ]
    reports_by_name = _diversification_model_by_name(diversification_reports)
    clade_report = detect_diversification_outlier_clades(
        tree_path,
        model=clade_model,
    )
    warnings = sorted(
        set(
            clade_report.warnings
            + reports_by_name["yule"].warnings
            + reports_by_name["birth-death"].warnings
            + ([] if sampling_report is None else sampling_report.warnings)
        )
    )
    rows = [
        _build_diversification_row(
            observation=observation,
            tree_path=tree_path,
            metadata_path=metadata_path,
            clade_model=clade_model,
            better_model=comparison.better_model,
            global_diversification_rate=clade_report.global_rate,
            yule_report=reports_by_name["yule"],
            birth_death_report=reports_by_name["birth-death"],
            sampling_report=sampling_report,
            warnings=warnings,
        )
        for observation in clade_report.observations
    ]
    columns = diversification_table_columns()
    _write_diversification_rows(path, columns=columns, rows=rows)
    return SupplementaryDiversificationTableResult(
        output_path=path,
        row_count=len(rows),
        better_model=comparison.better_model,
        clade_model=clade_model,
        high_clade_count=len(clade_report.high_diversification_clades),
        low_clade_count=len(clade_report.low_diversification_clades),
        warning_count=len(warnings),
        sampling_metadata_complete=(
            None if sampling_report is None else sampling_report.complete
        ),
        columns=columns,
        rows=rows,
    )
