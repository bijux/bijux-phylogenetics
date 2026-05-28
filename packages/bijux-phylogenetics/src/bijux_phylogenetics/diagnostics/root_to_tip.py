from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.diagnostics.validation import _load_tree, validate_tree_path
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    MetadataJoinError,
    PhylogeneticsError,
)


@dataclass(slots=True)
class RootToTipDistance:
    tip: str
    distance: float | None


@dataclass(slots=True)
class RootToTipDistanceReport:
    path: Path
    source_format: str
    distances: list[RootToTipDistance]


@dataclass(slots=True)
class UltrametricReport:
    path: Path
    source_format: str
    tolerance: float
    ultrametric: bool | None
    max_deviation: float | None
    tip_count: int


@dataclass(slots=True)
class RootToTipRegressionRow:
    tip: str
    sampling_time: float
    root_to_tip_distance: float
    fitted_distance: float
    residual: float
    leverage: float
    studentized_residual: float
    outlier: bool


@dataclass(slots=True)
class RootToTipOutlier:
    tip: str
    sampling_time: float
    root_to_tip_distance: float
    residual: float
    studentized_residual: float
    leverage: float
    rank: int


@dataclass(slots=True)
class RootToTipRegressionReport:
    tree_path: Path
    metadata_path: Path
    source_format: str
    taxon_column: str
    date_column: str
    tip_count: int
    slope: float
    intercept: float
    r_squared: float
    residual_mean_square: float
    outlier_threshold: float
    sampling_time_min: float
    sampling_time_max: float
    root_to_tip_min: float
    root_to_tip_max: float
    rows: list[RootToTipRegressionRow]
    outliers: list[RootToTipOutlier]


def compute_root_to_tip_distances(
    path: Path, *, source_format: str | None = None
) -> RootToTipDistanceReport:
    """Compute one root-to-tip distance per leaf."""
    tree = _load_tree(path, source_format=source_format)
    return RootToTipDistanceReport(
        path=path,
        source_format=tree.source_format,
        distances=[
            RootToTipDistance(tip=tip_name, distance=distance)
            for tip_name, distance in tree.root_to_tip_pairs()
            if tip_name is not None
        ],
    )


def write_root_to_tip_tsv(path: Path, report: RootToTipDistanceReport) -> Path:
    """Write root-to-tip distances as a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["tip\tdistance"]
    lines.extend(
        f"{row.tip}\t{'' if row.distance is None else format(row.distance, '.15g')}"
        for row in report.distances
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def diagnose_root_to_tip_regression(
    tree_path: Path,
    metadata_path: Path,
    *,
    source_format: str | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
    outlier_threshold: float = 2.0,
) -> RootToTipRegressionReport:
    """Regress root-to-tip distance against one numeric sampling-time column."""
    validate_tree_path(tree_path, source_format=source_format, require_rooted=True)
    distance_report = compute_root_to_tip_distances(
        tree_path,
        source_format=source_format,
    )
    if any(row.distance is None for row in distance_report.distances):
        raise InvalidBranchLengthError(
            "root-to-tip regression requires complete branch lengths"
        )
    sampling_times, resolved_taxon_column = _load_sampling_times(
        metadata_path,
        tree_taxa=[row.tip for row in distance_report.distances],
        taxon_column=taxon_column,
        date_column=date_column,
    )
    if len(distance_report.distances) < 3:
        raise PhylogeneticsError(
            "root-to-tip regression requires at least three dated tips",
            code="root_to_tip_regression_error",
        )
    observations = [
        (
            row.tip,
            sampling_times[row.tip],
            float(row.distance),
        )
        for row in distance_report.distances
    ]
    sampling_values = [sampling_time for _tip, sampling_time, _distance in observations]
    distance_values = [distance for _tip, _sampling_time, distance in observations]
    mean_sampling_time = sum(sampling_values) / len(sampling_values)
    mean_distance = sum(distance_values) / len(distance_values)
    sum_squares_sampling = sum(
        (sampling_time - mean_sampling_time) ** 2
        for sampling_time in sampling_values
    )
    if sum_squares_sampling == 0.0:
        raise PhylogeneticsError(
            "root-to-tip regression requires variation in sampling time",
            code="root_to_tip_regression_error",
        )
    sum_products = sum(
        (sampling_time - mean_sampling_time) * (distance - mean_distance)
        for _tip, sampling_time, distance in observations
    )
    slope = sum_products / sum_squares_sampling
    intercept = mean_distance - (slope * mean_sampling_time)
    fitted_distances = [
        intercept + (slope * sampling_time) for sampling_time in sampling_values
    ]
    residuals = [
        distance - fitted_distance
        for distance, fitted_distance in zip(
            distance_values,
            fitted_distances,
            strict=True,
        )
    ]
    residual_sum_squares = sum(residual * residual for residual in residuals)
    total_sum_squares = sum(
        (distance - mean_distance) ** 2 for distance in distance_values
    )
    r_squared = (
        1.0 - (residual_sum_squares / total_sum_squares)
        if total_sum_squares
        else 1.0
    )
    residual_mean_square = residual_sum_squares / (len(observations) - 2)
    rows: list[RootToTipRegressionRow] = []
    for tip, sampling_time, distance, fitted_distance, residual in zip(
        [tip for tip, _sampling_time, _distance in observations],
        sampling_values,
        distance_values,
        fitted_distances,
        residuals,
        strict=True,
    ):
        leverage = min(
            max(
                (1.0 / len(observations))
                + (((sampling_time - mean_sampling_time) ** 2) / sum_squares_sampling),
                0.0,
            ),
            0.999999,
        )
        denominator = math.sqrt(
            max(residual_mean_square * (1.0 - leverage), 1e-12)
        )
        studentized_residual = residual / denominator
        rows.append(
            RootToTipRegressionRow(
                tip=tip,
                sampling_time=sampling_time,
                root_to_tip_distance=distance,
                fitted_distance=fitted_distance,
                residual=residual,
                leverage=leverage,
                studentized_residual=studentized_residual,
                outlier=abs(studentized_residual) >= outlier_threshold,
            )
        )
    ranked_outliers = sorted(
        (row for row in rows if row.outlier),
        key=lambda row: (
            -abs(row.studentized_residual),
            row.tip,
        ),
    )
    outliers = [
        RootToTipOutlier(
            tip=row.tip,
            sampling_time=row.sampling_time,
            root_to_tip_distance=row.root_to_tip_distance,
            residual=row.residual,
            studentized_residual=row.studentized_residual,
            leverage=row.leverage,
            rank=index,
        )
        for index, row in enumerate(ranked_outliers, start=1)
    ]
    return RootToTipRegressionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        source_format=distance_report.source_format,
        taxon_column=resolved_taxon_column,
        date_column=date_column,
        tip_count=len(observations),
        slope=slope,
        intercept=intercept,
        r_squared=r_squared,
        residual_mean_square=residual_mean_square,
        outlier_threshold=outlier_threshold,
        sampling_time_min=min(sampling_values),
        sampling_time_max=max(sampling_values),
        root_to_tip_min=min(distance_values),
        root_to_tip_max=max(distance_values),
        rows=rows,
        outliers=outliers,
    )


def write_root_to_tip_regression_summary_tsv(
    path: Path,
    report: RootToTipRegressionReport,
) -> Path:
    """Write one summary row for a root-to-tip regression run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tree_path",
        "metadata_path",
        "source_format",
        "taxon_column",
        "date_column",
        "tip_count",
        "slope",
        "intercept",
        "r_squared",
        "residual_mean_square",
        "outlier_threshold",
        "outlier_count",
        "outlier_tips",
        "sampling_time_min",
        "sampling_time_max",
        "root_to_tip_min",
        "root_to_tip_max",
    ]
    values = [
        str(report.tree_path),
        str(report.metadata_path),
        report.source_format,
        report.taxon_column,
        report.date_column,
        str(report.tip_count),
        format(report.slope, ".15g"),
        format(report.intercept, ".15g"),
        format(report.r_squared, ".15g"),
        format(report.residual_mean_square, ".15g"),
        format(report.outlier_threshold, ".15g"),
        str(len(report.outliers)),
        "|".join(outlier.tip for outlier in report.outliers),
        format(report.sampling_time_min, ".15g"),
        format(report.sampling_time_max, ".15g"),
        format(report.root_to_tip_min, ".15g"),
        format(report.root_to_tip_max, ".15g"),
    ]
    path.write_text(
        "\n".join(["\t".join(columns), "\t".join(values)]) + "\n",
        encoding="utf-8",
    )
    return path


def write_root_to_tip_regression_residuals_tsv(
    path: Path,
    report: RootToTipRegressionReport,
) -> Path:
    """Write one residual row per tip in tree order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "tip",
        "sampling_time",
        "root_to_tip_distance",
        "fitted_distance",
        "residual",
        "leverage",
        "studentized_residual",
        "outlier",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                row.tip,
                format(row.sampling_time, ".15g"),
                format(row.root_to_tip_distance, ".15g"),
                format(row.fitted_distance, ".15g"),
                format(row.residual, ".15g"),
                format(row.leverage, ".15g"),
                format(row.studentized_residual, ".15g"),
                str(row.outlier).lower(),
            ]
        )
        for row in report.rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_root_to_tip_regression_outliers_tsv(
    path: Path,
    report: RootToTipRegressionReport,
) -> Path:
    """Write only outlier tips ranked by absolute studentized residual."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "rank",
        "tip",
        "sampling_time",
        "root_to_tip_distance",
        "residual",
        "studentized_residual",
        "leverage",
    ]
    lines = ["\t".join(columns)]
    lines.extend(
        "\t".join(
            [
                str(row.rank),
                row.tip,
                format(row.sampling_time, ".15g"),
                format(row.root_to_tip_distance, ".15g"),
                format(row.residual, ".15g"),
                format(row.studentized_residual, ".15g"),
                format(row.leverage, ".15g"),
            ]
        )
        for row in report.outliers
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_root_to_tip_regression_run_json(
    path: Path,
    report: RootToTipRegressionReport,
) -> Path:
    """Write the full root-to-tip regression payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "tree_path": str(report.tree_path),
        "metadata_path": str(report.metadata_path),
        "source_format": report.source_format,
        "taxon_column": report.taxon_column,
        "date_column": report.date_column,
        "tip_count": report.tip_count,
        "slope": report.slope,
        "intercept": report.intercept,
        "r_squared": report.r_squared,
        "residual_mean_square": report.residual_mean_square,
        "outlier_threshold": report.outlier_threshold,
        "sampling_time_min": report.sampling_time_min,
        "sampling_time_max": report.sampling_time_max,
        "root_to_tip_min": report.root_to_tip_min,
        "root_to_tip_max": report.root_to_tip_max,
        "rows": [asdict(row) for row in report.rows],
        "outliers": [asdict(row) for row in report.outliers],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_root_to_tip_regression_artifacts(
    out_dir: Path,
    report: RootToTipRegressionReport,
) -> dict[str, Path]:
    """Write governed artifact outputs for one root-to-tip regression run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = write_root_to_tip_regression_summary_tsv(
        out_dir / "summary.tsv",
        report,
    )
    residuals_path = write_root_to_tip_regression_residuals_tsv(
        out_dir / "residuals.tsv",
        report,
    )
    outliers_path = write_root_to_tip_regression_outliers_tsv(
        out_dir / "outliers.tsv",
        report,
    )
    run_json_path = write_root_to_tip_regression_run_json(out_dir / "run.json", report)
    return {
        "summary": summary_path,
        "residuals": residuals_path,
        "outliers": outliers_path,
        "run_json": run_json_path,
    }


def diagnose_ultrametricity(
    path: Path,
    *,
    source_format: str | None = None,
    tolerance: float = 1e-6,
) -> UltrametricReport:
    """Assess whether a tree is ultrametric within the given tolerance."""
    report = compute_root_to_tip_distances(path, source_format=source_format)
    if any(row.distance is None for row in report.distances):
        return UltrametricReport(
            path=path,
            source_format=report.source_format,
            tolerance=tolerance,
            ultrametric=None,
            max_deviation=None,
            tip_count=len(report.distances),
        )
    summary = summarize_ultrametric_tip_depths(
        {
            row.tip: float(row.distance)
            for row in report.distances
            if row.distance is not None
        },
        tolerance=tolerance,
    )
    return UltrametricReport(
        path=path,
        source_format=report.source_format,
        tolerance=tolerance,
        ultrametric=summary.ultrametric,
        max_deviation=round(summary.max_tip_depth_deviation, 15),
        tip_count=len(report.distances),
    )


def _load_sampling_times(
    metadata_path: Path,
    *,
    tree_taxa: list[str],
    taxon_column: str | None,
    date_column: str,
) -> tuple[dict[str, float], str]:
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    if date_column not in table.columns:
        raise MetadataJoinError(
            f"missing date column '{date_column}' in {metadata_path}"
        )
    tree_taxa_set = set(tree_taxa)
    table_taxa_set = set(table.taxa)
    missing_tree_taxa = sorted(tree_taxa_set - table_taxa_set)
    if missing_tree_taxa:
        raise MetadataJoinError(
            "sampling-time table is missing tree taxa: "
            + ", ".join(missing_tree_taxa)
        )
    extra_table_taxa = sorted(table_taxa_set - tree_taxa_set)
    if extra_table_taxa:
        raise MetadataJoinError(
            "sampling-time table contains taxa absent from the tree: "
            + ", ".join(extra_table_taxa)
        )
    sampling_times: dict[str, float] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_value = row.get(date_column, "")
        if not raw_value:
            raise MetadataJoinError(
                f"taxon '{taxon}' is missing a numeric sampling time in column '{date_column}'"
            )
        try:
            sampling_times[taxon] = float(raw_value)
        except ValueError as error:
            raise MetadataJoinError(
                f"taxon '{taxon}' has a non-numeric sampling time '{raw_value}'"
            ) from error
    return sampling_times, table.taxon_column
