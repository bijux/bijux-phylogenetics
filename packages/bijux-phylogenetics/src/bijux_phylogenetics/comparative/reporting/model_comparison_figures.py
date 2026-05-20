from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from html import escape
import json
import math
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)

from ..models import (
    BrownianMotionFitReport,
    ComparativeModelComparisonReport,
    OUTraitModelReport,
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureLegendEntry:
    """One legend row for the comparative model-comparison package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureCaptionDraft:
    """Structured caption draft for one comparative model-comparison package."""

    title: str
    lead_sentence: str
    criteria_sentence: str
    likelihood_sentence: str
    parameter_sentence: str
    fit_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelCriteriaRow:
    """Reviewer-facing information-criterion summary for one fitted model."""

    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    delta_aic: float
    delta_aicc: float
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelLikelihoodRow:
    """Reviewer-facing likelihood summary for one fitted model."""

    model: str
    log_likelihood: float
    delta_log_likelihood: float
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelParameterRow:
    """Reviewer-facing parameter summary for one fitted model parameter."""

    model: str
    parameter: str
    estimate: float
    lower_95: float
    upper_95: float
    interval_method: str


@dataclass(frozen=True, slots=True)
class ComparativeModelFitRow:
    """Reviewer-facing fit diagnostics for one fitted comparative model."""

    model: str
    taxon_count: int
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float
    outlier_taxon_count: int
    warning_count: int
    convergence_status: str
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureAudit:
    """Publication-oriented audit for a comparative model-comparison package."""

    publication_ready: bool
    criteria_surface_visible: bool
    likelihood_surface_visible: bool
    parameter_surface_visible: bool
    fit_surface_visible: bool
    legend_complete: bool
    caption_ready: bool
    finite_aicc_model_count: int
    support_distinct: bool
    selected_model: str
    aicc_delta: float | None
    plotted_model_count: int
    rendered_parameter_count: int
    rendered_fit_row_count: int
    warning_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class ComparativeModelFigurePackageResult:
    output_dir: Path
    criteria_figure_path: Path
    likelihood_figure_path: Path
    parameter_figure_path: Path
    fit_figure_path: Path
    criteria_table_path: Path
    likelihood_table_path: Path
    parameter_table_path: Path
    fit_table_path: Path
    legend_path: Path
    caption_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    comparison_report: ComparativeModelComparisonReport
    brownian_report: BrownianMotionFitReport
    ou_report: OUTraitModelReport
    criteria_rows: list[ComparativeModelCriteriaRow]
    likelihood_rows: list[ComparativeModelLikelihoodRow]
    parameter_rows: list[ComparativeModelParameterRow]
    fit_rows: list[ComparativeModelFitRow]
    legend_entries: list[ComparativeModelFigureLegendEntry]
    caption_draft: ComparativeModelFigureCaptionDraft
    audit: ComparativeModelFigureAudit
    machine_manifest: dict[str, object]


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _format_number(value: float | None) -> str:
    if value is None:
        return "n/a"
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0.0 else "-inf"
    return format(value, ".15g")


def _delta_from_best(value: float, best: float) -> float:
    if math.isfinite(value) and math.isfinite(best):
        return value - best
    if math.isinf(value) and math.isinf(best):
        return 0.0
    if math.isinf(value):
        return math.inf
    if math.isinf(best):
        return -math.inf
    return math.nan


def _model_color(model: str, *, selected_model: str) -> str:
    if model == selected_model:
        return "#0f766e"
    if model == "brownian":
        return "#0369a1"
    return "#9333ea"


def _parameter_label(parameter: str) -> str:
    return {
        "root_state": "root state",
        "rate": "brownian rate",
        "alpha": "ou alpha",
        "theta": "ou theta",
        "sigma_squared": "ou sigma^2",
    }.get(parameter, parameter.replace("_", " "))


def _criteria_rows(
    report: ComparativeModelComparisonReport,
) -> list[ComparativeModelCriteriaRow]:
    best_aic = min(row.aic for row in report.rows)
    best_aicc = min(row.aicc for row in report.rows)
    return [
        ComparativeModelCriteriaRow(
            model=row.model,
            parameter_count=row.parameter_count,
            log_likelihood=row.log_likelihood,
            aic=row.aic,
            aicc=row.aicc,
            delta_aic=_delta_from_best(row.aic, best_aic),
            delta_aicc=_delta_from_best(row.aicc, best_aicc),
            selected=row.selected,
        )
        for row in report.rows
    ]


def _likelihood_rows(
    report: ComparativeModelComparisonReport,
) -> list[ComparativeModelLikelihoodRow]:
    best_log_likelihood = max(row.log_likelihood for row in report.rows)
    return [
        ComparativeModelLikelihoodRow(
            model=row.model,
            log_likelihood=row.log_likelihood,
            delta_log_likelihood=best_log_likelihood - row.log_likelihood,
            selected=row.selected,
        )
        for row in report.rows
    ]


def _parameter_rows(
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
) -> list[ComparativeModelParameterRow]:
    rows: list[ComparativeModelParameterRow] = []
    for interval in brownian.confidence_intervals:
        rows.append(
            ComparativeModelParameterRow(
                model="brownian",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                interval_method=interval.method,
            )
        )
    for interval in ou.confidence_intervals:
        rows.append(
            ComparativeModelParameterRow(
                model="ou",
                parameter=interval.name,
                estimate=interval.estimate,
                lower_95=interval.lower_95,
                upper_95=interval.upper_95,
                interval_method=interval.method,
            )
        )
    return rows


def _fit_rows(
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
    *,
    selected_model: str,
) -> list[ComparativeModelFitRow]:
    return [
        ComparativeModelFitRow(
            model="brownian",
            taxon_count=brownian.taxon_count,
            residual_variance=brownian.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=(
                brownian.residual_diagnostics.max_abs_standardized_residual
            ),
            phylogenetic_residual_lambda=(
                brownian.residual_diagnostics.phylogenetic_residual_lambda
            ),
            outlier_taxon_count=len(brownian.residual_diagnostics.outlier_taxa),
            warning_count=len(brownian.residual_diagnostics.warnings),
            convergence_status="analytic",
            selected=selected_model == "brownian",
        ),
        ComparativeModelFitRow(
            model="ou",
            taxon_count=ou.taxon_count,
            residual_variance=ou.residual_diagnostics.residual_variance,
            max_abs_standardized_residual=(
                ou.residual_diagnostics.max_abs_standardized_residual
            ),
            phylogenetic_residual_lambda=(
                ou.residual_diagnostics.phylogenetic_residual_lambda
            ),
            outlier_taxon_count=len(ou.residual_diagnostics.outlier_taxa),
            warning_count=len(ou.residual_diagnostics.warnings)
            + len(ou.identifiability_warnings),
            convergence_status=ou.convergence_status,
            selected=selected_model == "ou",
        ),
    ]


def _legend_entries() -> list[ComparativeModelFigureLegendEntry]:
    return [
        ComparativeModelFigureLegendEntry(
            surface="information-criteria",
            label="selected model",
            swatch="#0f766e",
            detail="the highlighted row marks the model selected by AICc when the candidate models retain finite information-criterion values",
        ),
        ComparativeModelFigureLegendEntry(
            surface="likelihood",
            label="log-likelihood bar",
            swatch="#0369a1",
            detail="bar length tracks relative log-likelihood so the gap between candidate models remains visible instead of being hidden inside a summary sentence",
        ),
        ComparativeModelFigureLegendEntry(
            surface="parameters",
            label="95% parameter interval",
            swatch="#9333ea",
            detail="horizontal whiskers retain the estimated parameter value with its fitted 95% interval for each continuous-trait model parameter",
        ),
        ComparativeModelFigureLegendEntry(
            surface="fit-summary",
            label="diagnostic row",
            swatch="#475569",
            detail="the fit-summary surface keeps residual scale, standardized outlier burden, phylogenetic residual lambda, and warning counts visible for both models",
        ),
    ]


def _write_table(path: Path, rows: list[object]) -> Path:
    return write_taxon_rows(
        path,
        columns=list(asdict(rows[0]).keys()),
        rows=[asdict(row) for row in rows],
    )


def _write_legend_table(
    path: Path, entries: list[ComparativeModelFigureLegendEntry]
) -> Path:
    return write_taxon_rows(
        path,
        columns=["surface", "label", "swatch", "detail"],
        rows=[
            {
                "surface": entry.surface,
                "label": entry.label,
                "swatch": entry.swatch,
                "detail": entry.detail,
            }
            for entry in entries
        ],
    )


def _write_information_criteria_svg(
    path: Path,
    rows: list[ComparativeModelCriteriaRow],
    *,
    selected_model: str,
) -> int:
    width = 1080
    row_height = 52
    top = 84
    height = top + len(rows) * row_height + 48
    columns = [
        ("model", 56),
        ("logL", 264),
        ("AIC", 430),
        ("AICc", 566),
        ("delta AIC", 710),
        ("delta AICc", 866),
        ("selected", 992),
    ]
    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model information criteria">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        '<text x="56" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Comparative model information criteria</text>',
        '<text x="56" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">AIC and AICc remain explicit so the selected model is reviewable from the package itself.</text>',
        '<rect x="40" y="64" width="1000" height="36" rx="14" fill="#e2e8f0" />',
    ]
    for label, x in columns:
        anchor = "middle" if label != "model" else "start"
        segments.append(
            f'<text x="{x}" y="87" text-anchor="{anchor}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#334155">{escape(label)}</text>'
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        fill = "#ecfeff" if row.model == selected_model else "#ffffff"
        stroke = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="40" y="{y - 24}" width="1000" height="40" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="56" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<text x="264" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.log_likelihood))}</text>',
                f'<text x="430" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.aic))}</text>',
                f'<text x="566" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.aicc))}</text>',
                f'<text x="710" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.delta_aic))}</text>',
                f'<text x="866" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.delta_aicc))}</text>',
                f'<text x="992" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{"yes" if row.selected else "no"}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_likelihood_svg(
    path: Path,
    rows: list[ComparativeModelLikelihoodRow],
    *,
    selected_model: str,
) -> int:
    width = 980
    height = max(250, 130 + len(rows) * 84)
    left = 240
    right = 52
    top = 52
    bottom = 48
    plot_width = width - left - right
    max_log_likelihood = max(row.log_likelihood for row in rows)
    min_log_likelihood = min(row.log_likelihood for row in rows)
    span = max(max_log_likelihood - min_log_likelihood, 1e-9)

    def x_position(value: float) -> float:
        return left + ((value - min_log_likelihood) / span) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model likelihood comparison">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        '<text x="48" y="26" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Comparative model likelihood comparison</text>',
        '<text x="48" y="50" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">Higher log-likelihood bars indicate better likelihood fit before the information-criterion penalty is applied.</text>',
        f'<line x1="{left}" y1="{height - bottom}" x2="{left + plot_width}" y2="{height - bottom}" stroke="#cbd5e1" stroke-width="1.5" />',
    ]
    for tick_value in (min_log_likelihood, (min_log_likelihood + max_log_likelihood) / 2.0, max_log_likelihood):
        x = x_position(tick_value)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{height - bottom}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{x}" y="{height - 12}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{escape(_format_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(sorted(rows, key=lambda item: item.log_likelihood)):
        y = top + index * 84 + 24
        bar_end = x_position(row.log_likelihood)
        color = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<text x="{left - 16}" y="{y + 4}" text-anchor="end" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<rect x="{left}" y="{y - 10}" width="{max(bar_end - left, 4):.2f}" height="20" rx="10" fill="{color}" opacity="0.9" />',
                f'<text x="{bar_end + 12}" y="{y + 4}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">logL={escape(_format_number(row.log_likelihood))}; gap={escape(_format_number(row.delta_log_likelihood))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_parameter_svg(
    path: Path,
    rows: list[ComparativeModelParameterRow],
    *,
    selected_model: str,
) -> int:
    width = 1120
    row_height = 74
    left = 310
    right = 64
    top = 72
    bottom = 48
    height = top + len(rows) * row_height + bottom
    plot_width = width - left - right
    minimum = min(min(row.lower_95, row.estimate) for row in rows)
    maximum = max(max(row.upper_95, row.estimate) for row in rows)
    if math.isclose(minimum, maximum, abs_tol=1e-12):
        minimum -= 1.0
        maximum += 1.0
    padding = max((maximum - minimum) * 0.08, 1e-6)
    minimum -= padding
    maximum += padding
    span = maximum - minimum

    def x_position(value: float) -> float:
        return left + ((value - minimum) / span) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model parameter intervals">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#faf5ff" />',
        '<text x="56" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#3b0764">Comparative model parameter summaries</text>',
        '<text x="56" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#6b21a8">Each parameter keeps its fitted estimate and interval visible so the model comparison stays interpretable beyond the winner call.</text>',
    ]
    for tick_value in (minimum, (minimum + maximum) / 2.0, maximum):
        x = x_position(tick_value)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top - 8}" x2="{x}" y2="{height - bottom}" stroke="#e9d5ff" stroke-width="1" />',
                f'<text x="{x}" y="{height - 14}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#6b21a8">{escape(_format_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        color = _model_color(row.model, selected_model=selected_model)
        lower = x_position(row.lower_95)
        estimate = x_position(row.estimate)
        upper = x_position(row.upper_95)
        segments.extend(
            [
                f'<text x="56" y="{y + 6}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#3b0764">{escape(row.model)}</text>',
                f'<text x="56" y="{y + 28}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#581c87">{escape(_parameter_label(row.parameter))}</text>',
                f'<text x="56" y="{y + 48}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#7e22ce">{escape(row.interval_method)}</text>',
                f'<line x1="{lower:.2f}" y1="{y + 22}" x2="{upper:.2f}" y2="{y + 22}" stroke="{color}" stroke-width="4" stroke-linecap="round" />',
                f'<line x1="{lower:.2f}" y1="{y + 12}" x2="{lower:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<line x1="{upper:.2f}" y1="{y + 12}" x2="{upper:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<circle cx="{estimate:.2f}" cy="{y + 22}" r="7" fill="{color}" stroke="#ffffff" stroke-width="1.5" />',
                f'<text x="{estimate:.2f}" y="{y + 54}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">{escape(_format_number(row.estimate))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_fit_summary_svg(
    path: Path,
    rows: list[ComparativeModelFitRow],
    *,
    selected_model: str,
) -> int:
    width = 1140
    row_height = 56
    top = 84
    height = top + len(rows) * row_height + 48
    columns = [
        ("model", 60),
        ("residual variance", 286),
        ("max |z|", 474),
        ("residual lambda", 644),
        ("outliers", 816),
        ("warnings", 934),
        ("status", 1042),
    ]
    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model fit summary">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fffaf0" />',
        '<text x="60" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">Comparative model fit summaries</text>',
        '<text x="60" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#9a3412">Residual scale, standardized outlier burden, and warning counts remain reviewable alongside the information criteria.</text>',
        '<rect x="44" y="64" width="1050" height="36" rx="14" fill="#fed7aa" />',
    ]
    for label, x in columns:
        anchor = "middle" if label != "model" else "start"
        segments.append(
            f'<text x="{x}" y="87" text-anchor="{anchor}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">{escape(label)}</text>'
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        fill = "#fffbeb" if row.model == selected_model else "#ffffff"
        stroke = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="44" y="{y - 26}" width="1050" height="42" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="60" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#451a03">{escape(row.model)}</text>',
                f'<text x="286" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.residual_variance))}</text>',
                f'<text x="474" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.max_abs_standardized_residual))}</text>',
                f'<text x="644" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.phylogenetic_residual_lambda))}</text>',
                f'<text x="816" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.outlier_taxon_count}</text>',
                f'<text x="934" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.warning_count}</text>',
                f'<text x="1042" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{escape(row.convergence_status)}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _build_audit(
    *,
    comparison_report: ComparativeModelComparisonReport,
    brownian: BrownianMotionFitReport,
    ou: OUTraitModelReport,
    criteria_rows: list[ComparativeModelCriteriaRow],
    likelihood_rows: list[ComparativeModelLikelihoodRow],
    parameter_rows: list[ComparativeModelParameterRow],
    fit_rows: list[ComparativeModelFitRow],
    legend_entries: list[ComparativeModelFigureLegendEntry],
) -> ComparativeModelFigureAudit:
    criteria_surface_visible = len(criteria_rows) == len(comparison_report.rows)
    likelihood_surface_visible = len(likelihood_rows) == len(comparison_report.rows)
    parameter_surface_visible = len(parameter_rows) >= 5
    fit_surface_visible = len(fit_rows) == len(comparison_report.rows)
    legend_complete = {
        entry.surface for entry in legend_entries
    } == {"information-criteria", "likelihood", "parameters", "fit-summary"}
    caption_ready = (
        criteria_surface_visible
        and likelihood_surface_visible
        and parameter_surface_visible
        and fit_surface_visible
    )
    finite_rows = [row for row in criteria_rows if math.isfinite(row.aicc)]
    selected_row = next(row for row in criteria_rows if row.model == comparison_report.better_model)
    runner_up_candidates = [
        row.aicc for row in criteria_rows if row.model != comparison_report.better_model
    ]
    runner_up_aicc = min(runner_up_candidates) if runner_up_candidates else math.nan
    aicc_delta = (
        runner_up_aicc - selected_row.aicc
        if math.isfinite(runner_up_aicc) and math.isfinite(selected_row.aicc)
        else None
    )
    support_distinct = aicc_delta is not None and aicc_delta >= 2.0
    warning_count = len(brownian.residual_diagnostics.warnings) + len(
        ou.residual_diagnostics.warnings
    ) + len(ou.identifiability_warnings)
    publication_ready = (
        criteria_surface_visible
        and likelihood_surface_visible
        and parameter_surface_visible
        and fit_surface_visible
        and legend_complete
        and caption_ready
        and len(finite_rows) == len(criteria_rows)
        and support_distinct
    )
    reviewer_summary = [
        f"criteria rows rendered: {len(criteria_rows)}/{len(comparison_report.rows)}",
        f"likelihood rows rendered: {len(likelihood_rows)}/{len(comparison_report.rows)}",
        f"parameter rows rendered: {len(parameter_rows)}",
        f"fit rows rendered: {len(fit_rows)}/{len(comparison_report.rows)}",
        f"selected model: {comparison_report.better_model}",
        f"aicc delta to runner-up: {_format_number(aicc_delta)}",
    ]
    limitations: list[str] = []
    if len(finite_rows) != len(criteria_rows):
        limitations.append(
            "one or more candidate models do not retain finite AICc values, so the model-comparison package stays reviewable but not publication-ready"
        )
    if not support_distinct:
        limitations.append(
            "the AICc separation between the selected model and the runner-up remains below the publication threshold of 2.0, so the package does not claim a decisively supported winner"
        )
    if not criteria_surface_visible:
        limitations.append(
            "the information-criterion surface does not retain a complete model row for every candidate model"
        )
    if not likelihood_surface_visible:
        limitations.append(
            "the likelihood surface does not retain a complete log-likelihood row for every candidate model"
        )
    if not parameter_surface_visible:
        limitations.append(
            "the parameter surface does not retain the fitted BM and OU parameter intervals"
        )
    if not fit_surface_visible:
        limitations.append(
            "the fit-summary surface does not retain both fitted model diagnostic rows"
        )
    if not legend_complete:
        limitations.append(
            "the figure legend does not cover all four rendered model-comparison surfaces"
        )
    if not limitations:
        limitations.append(
            "the current package keeps information criteria, likelihood, parameters, and fit diagnostics explicit enough for publication-oriented comparative review"
        )
    return ComparativeModelFigureAudit(
        publication_ready=publication_ready,
        criteria_surface_visible=criteria_surface_visible,
        likelihood_surface_visible=likelihood_surface_visible,
        parameter_surface_visible=parameter_surface_visible,
        fit_surface_visible=fit_surface_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        finite_aicc_model_count=len(finite_rows),
        support_distinct=support_distinct,
        selected_model=comparison_report.better_model,
        aicc_delta=aicc_delta,
        plotted_model_count=len(comparison_report.rows),
        rendered_parameter_count=len(parameter_rows),
        rendered_fit_row_count=len(fit_rows),
        warning_count=warning_count,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def _build_caption_draft(
    *,
    comparison_report: ComparativeModelComparisonReport,
    audit: ComparativeModelFigureAudit,
    fit_rows: list[ComparativeModelFitRow],
) -> ComparativeModelFigureCaptionDraft:
    selected_fit = next(row for row in fit_rows if row.model == audit.selected_model)
    return ComparativeModelFigureCaptionDraft(
        title="Comparative continuous-trait model review across information criteria, likelihood, parameters, and fit diagnostics",
        lead_sentence=(
            f"This package compares Brownian-motion and Ornstein-Uhlenbeck fits for `{comparison_report.trait}` across {comparison_report.taxon_count} overlapping taxa with reviewer-facing figure surfaces instead of leaving model choice buried in one JSON payload."
        ),
        criteria_sentence=(
            f"The information-criterion surface retains AIC and AICc for both candidate models and currently selects `{audit.selected_model}` with an AICc delta of {_format_number(audit.aicc_delta)} relative to the runner-up."
        ),
        likelihood_sentence=(
            "The likelihood surface keeps the relative log-likelihood gap explicit so likelihood fit can be reviewed separately from the information-criterion penalty."
        ),
        parameter_sentence=(
            f"The parameter surface renders {audit.rendered_parameter_count} fitted BM and OU parameter intervals, including Brownian root-state and rate estimates plus OU alpha, theta, and sigma-squared."
        ),
        fit_sentence=(
            f"The fit-summary surface keeps residual variance, standardized outlier burden, phylogenetic residual lambda, and warning counts visible for the selected `{audit.selected_model}` model, which currently carries {selected_fit.warning_count} warning rows."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )


def _write_caption(path: Path, draft: ComparativeModelFigureCaptionDraft) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.criteria_sentence,
                draft.likelihood_sentence,
                draft.parameter_sentence,
                draft.fit_sentence,
                draft.limitation_sentence,
                "",
                f"caption_ready: {'true' if draft.caption_ready else 'false'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def _build_review_html(
    *,
    criteria_figure_path: Path,
    likelihood_figure_path: Path,
    parameter_figure_path: Path,
    fit_figure_path: Path,
    criteria_table_path: Path,
    likelihood_table_path: Path,
    parameter_table_path: Path,
    fit_table_path: Path,
    legend_path: Path,
    caption_path: Path,
    audit: ComparativeModelFigureAudit,
) -> str:
    figures = {
        "criteria": criteria_figure_path.read_text(encoding="utf-8"),
        "likelihood": likelihood_figure_path.read_text(encoding="utf-8"),
        "parameters": parameter_figure_path.read_text(encoding="utf-8"),
        "fit": fit_figure_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>"
        + escape(label)
        + "</th><td>"
        + escape(value)
        + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("selected_model", audit.selected_model),
            ("support_distinct", str(audit.support_distinct).lower()),
            ("aicc_delta", _format_number(audit.aicc_delta)),
            ("finite_aicc_model_count", str(audit.finite_aicc_model_count)),
            ("warning_count", str(audit.warning_count)),
        ]
    )
    limitation_items = "".join(f"<li>{escape(item)}</li>" for item in audit.limitations)
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Comparative Model Comparison Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #f8fafc 0%, #fffaf0 100%); color: #0f172a; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1280px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { margin-top: 0; color: #0f172a; }",
            "    .grid { display: grid; grid-template-columns: 1fr; gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.92); border: 1px solid rgba(15,23,42,0.08); border-radius: 18px; padding: 18px; box-shadow: 0 18px 40px rgba(15,23,42,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(148,163,184,0.35); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #0f766e; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Comparative Model Comparison Review</h1>",
            "  <p>Reviewer-facing figure package for continuous-trait model support across information criteria, likelihood, parameters, and fit diagnostics.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Information Criteria</h2><div class="figure-shell">' + figures["criteria"] + "</div></section>",
            '    <section class="panel"><h2>Likelihood</h2><div class="figure-shell">' + figures["likelihood"] + "</div></section>",
            '    <section class="panel"><h2>Parameters</h2><div class="figure-shell">' + figures["parameters"] + "</div></section>",
            '    <section class="panel"><h2>Fit Summary</h2><div class="figure-shell">' + figures["fit"] + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(criteria_table_path.name)}">{escape(criteria_table_path.name)}</a></li>',
            f'      <li><a href="{escape(likelihood_table_path.name)}">{escape(likelihood_table_path.name)}</a></li>',
            f'      <li><a href="{escape(parameter_table_path.name)}">{escape(parameter_table_path.name)}</a></li>',
            f'      <li><a href="{escape(fit_table_path.name)}">{escape(fit_table_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def build_comparative_model_figure_package(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    out_dir: Path,
    taxon_column: str | None = None,
) -> ComparativeModelFigurePackageResult:
    """Build a reviewer-facing figure package for BM versus OU continuous-trait model comparison."""
    out_dir.mkdir(parents=True, exist_ok=True)
    criteria_figure_path = out_dir / "model-comparison-criteria.svg"
    likelihood_figure_path = out_dir / "model-comparison-likelihood.svg"
    parameter_figure_path = out_dir / "model-comparison-parameters.svg"
    fit_figure_path = out_dir / "model-comparison-fit-summary.svg"
    criteria_table_path = out_dir / "model-comparison-criteria.tsv"
    likelihood_table_path = out_dir / "model-comparison-likelihood.tsv"
    parameter_table_path = out_dir / "model-comparison-parameters.tsv"
    fit_table_path = out_dir / "model-comparison-fit-summary.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "model-comparison-review.html"
    manifest_path = out_dir / "model-comparison-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    comparison_report = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    brownian_report = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    ou_report = fit_ornstein_uhlenbeck_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    criteria_rows = _criteria_rows(comparison_report)
    likelihood_rows = _likelihood_rows(comparison_report)
    parameter_rows = _parameter_rows(brownian_report, ou_report)
    fit_rows = _fit_rows(
        brownian_report,
        ou_report,
        selected_model=comparison_report.better_model,
    )
    _write_table(criteria_table_path, criteria_rows)
    _write_table(likelihood_table_path, likelihood_rows)
    _write_table(parameter_table_path, parameter_rows)
    _write_table(fit_table_path, fit_rows)
    _write_information_criteria_svg(
        criteria_figure_path,
        criteria_rows,
        selected_model=comparison_report.better_model,
    )
    _write_likelihood_svg(
        likelihood_figure_path,
        likelihood_rows,
        selected_model=comparison_report.better_model,
    )
    _write_parameter_svg(
        parameter_figure_path,
        parameter_rows,
        selected_model=comparison_report.better_model,
    )
    _write_fit_summary_svg(
        fit_figure_path,
        fit_rows,
        selected_model=comparison_report.better_model,
    )
    legend_entries = _legend_entries()
    _write_legend_table(legend_path, legend_entries)
    audit = _build_audit(
        comparison_report=comparison_report,
        brownian=brownian_report,
        ou=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
    )
    caption_draft = _build_caption_draft(
        comparison_report=comparison_report,
        audit=audit,
        fit_rows=fit_rows,
    )
    _write_caption(caption_path, caption_draft)
    review_path.write_text(
        _build_review_html(
            criteria_figure_path=criteria_figure_path,
            likelihood_figure_path=likelihood_figure_path,
            parameter_figure_path=parameter_figure_path,
            fit_figure_path=fit_figure_path,
            criteria_table_path=criteria_table_path,
            likelihood_table_path=likelihood_table_path,
            parameter_table_path=parameter_table_path,
            fit_table_path=fit_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            audit=audit,
        ),
        encoding="utf-8",
    )
    artifact_paths = [
        criteria_figure_path,
        likelihood_figure_path,
        parameter_figure_path,
        fit_figure_path,
        criteria_table_path,
        likelihood_table_path,
        parameter_table_path,
        fit_table_path,
        legend_path,
        caption_path,
        review_path,
    ]
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="comparative_model_figure_package",
        input_files=[
            ("tree", tree_path),
            ("traits", traits_path),
        ],
        generated_figures=[
            ("information_criteria", criteria_figure_path),
            ("likelihood", likelihood_figure_path),
            ("parameter_intervals", parameter_figure_path),
            ("fit_diagnostics", fit_figure_path),
        ],
        generated_tables=[
            ("information_criteria", criteria_table_path),
            ("likelihood", likelihood_table_path),
            ("parameter_intervals", parameter_table_path),
            ("fit_diagnostics", fit_table_path),
        ],
        filters=None,
        model={
            "kind": "comparative_model_comparison",
            "name": comparison_report.better_model,
            "candidate_models": [row.model for row in comparison_report.rows],
            "aicc_delta": audit.aicc_delta,
        },
        settings={
            "trait": trait,
            "taxon_column": taxon_column,
            "taxon_count": comparison_report.taxon_count,
        },
        linked_artifacts=[
            ("legend", legend_path),
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "comparative_model_figure_package",
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(
            reproducibility_manifest_path
        ),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "trait": trait,
            "taxon_column": taxon_column,
        },
        "metrics": {
            "taxon_count": comparison_report.taxon_count,
            "selected_model": audit.selected_model,
            "publication_ready": audit.publication_ready,
            "support_distinct": audit.support_distinct,
            "aicc_delta": audit.aicc_delta,
            "finite_aicc_model_count": audit.finite_aicc_model_count,
            "plotted_model_count": audit.plotted_model_count,
            "rendered_parameter_count": audit.rendered_parameter_count,
            "rendered_fit_row_count": audit.rendered_fit_row_count,
            "warning_count": audit.warning_count,
        },
        "comparison_report": _json_ready(asdict(comparison_report)),
        "brownian_report": _json_ready(asdict(brownian_report)),
        "ou_report": _json_ready(asdict(ou_report)),
        "criteria_rows": _json_ready([asdict(row) for row in criteria_rows]),
        "likelihood_rows": _json_ready([asdict(row) for row in likelihood_rows]),
        "parameter_rows": _json_ready([asdict(row) for row in parameter_rows]),
        "fit_rows": _json_ready([asdict(row) for row in fit_rows]),
        "audit": _json_ready(asdict(audit)),
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return ComparativeModelFigurePackageResult(
        output_dir=out_dir,
        criteria_figure_path=criteria_figure_path,
        likelihood_figure_path=likelihood_figure_path,
        parameter_figure_path=parameter_figure_path,
        fit_figure_path=fit_figure_path,
        criteria_table_path=criteria_table_path,
        likelihood_table_path=likelihood_table_path,
        parameter_table_path=parameter_table_path,
        fit_table_path=fit_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        comparison_report=comparison_report,
        brownian_report=brownian_report,
        ou_report=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        machine_manifest=machine_manifest,
    )
