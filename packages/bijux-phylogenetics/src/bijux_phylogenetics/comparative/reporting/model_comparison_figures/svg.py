from __future__ import annotations

from html import escape
import math
from pathlib import Path

from .contracts import (
    ComparativeModelCriteriaRow,
    ComparativeModelFitRow,
    ComparativeModelLikelihoodRow,
    ComparativeModelParameterRow,
)
from .summaries import (
    format_model_figure_number,
    model_figure_color,
    parameter_label,
)


def write_information_criteria_svg(
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
            f'<text x="{x}" y="87" text-anchor="{anchor}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(label)}</text>'
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        fill = "#ecfeff" if row.model == selected_model else "#ffffff"
        stroke = model_figure_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="40" y="{y - 26}" width="1000" height="42" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="56" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<text x="264" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(format_model_figure_number(row.log_likelihood))}</text>',
                f'<text x="430" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(format_model_figure_number(row.aic))}</text>',
                f'<text x="566" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(format_model_figure_number(row.aicc))}</text>',
                f'<text x="710" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(format_model_figure_number(row.delta_aic))}</text>',
                f'<text x="866" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(format_model_figure_number(row.delta_aicc))}</text>',
                f'<text x="992" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{"yes" if row.selected else "no"}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def write_likelihood_svg(
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
    for tick_value in (
        min_log_likelihood,
        (min_log_likelihood + max_log_likelihood) / 2.0,
        max_log_likelihood,
    ):
        x = x_position(tick_value)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{height - bottom}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{x}" y="{height - 12}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{escape(format_model_figure_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(sorted(rows, key=lambda item: item.log_likelihood)):
        y = top + index * 84 + 24
        bar_end = x_position(row.log_likelihood)
        color = model_figure_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<text x="{left - 16}" y="{y + 4}" text-anchor="end" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<rect x="{left}" y="{y - 10}" width="{max(bar_end - left, 4):.2f}" height="20" rx="10" fill="{color}" opacity="0.9" />',
                f'<text x="{bar_end + 12}" y="{y + 4}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">logL={escape(format_model_figure_number(row.log_likelihood))}; gap={escape(format_model_figure_number(row.delta_log_likelihood))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def write_parameter_svg(
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
                f'<text x="{x}" y="{height - 14}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#6b21a8">{escape(format_model_figure_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        color = model_figure_color(row.model, selected_model=selected_model)
        lower = x_position(row.lower_95)
        estimate = x_position(row.estimate)
        upper = x_position(row.upper_95)
        segments.extend(
            [
                f'<text x="56" y="{y + 6}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#3b0764">{escape(row.model)}</text>',
                f'<text x="56" y="{y + 28}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#581c87">{escape(parameter_label(row.parameter))}</text>',
                f'<text x="56" y="{y + 48}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#7e22ce">{escape(row.interval_method)}</text>',
                f'<line x1="{lower:.2f}" y1="{y + 22}" x2="{upper:.2f}" y2="{y + 22}" stroke="{color}" stroke-width="4" stroke-linecap="round" />',
                f'<line x1="{lower:.2f}" y1="{y + 12}" x2="{lower:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<line x1="{upper:.2f}" y1="{y + 12}" x2="{upper:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<circle cx="{estimate:.2f}" cy="{y + 22}" r="7" fill="{color}" stroke="#ffffff" stroke-width="1.5" />',
                f'<text x="{estimate:.2f}" y="{y + 54}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">{escape(format_model_figure_number(row.estimate))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def write_fit_summary_svg(
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
        stroke = model_figure_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="44" y="{y - 26}" width="1050" height="42" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="60" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#451a03">{escape(row.model)}</text>',
                f'<text x="286" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(format_model_figure_number(row.residual_variance))}</text>',
                f'<text x="474" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(format_model_figure_number(row.max_abs_standardized_residual))}</text>',
                f'<text x="644" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(format_model_figure_number(row.phylogenetic_residual_lambda))}</text>',
                f'<text x="816" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.outlier_taxon_count}</text>',
                f'<text x="934" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.warning_count}</text>',
                f'<text x="1042" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{escape(row.convergence_status)}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)
