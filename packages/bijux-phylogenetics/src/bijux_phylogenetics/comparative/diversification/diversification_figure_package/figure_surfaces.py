from __future__ import annotations

from html import escape
from pathlib import Path

from ..models import (
    CladeDiversificationObservation,
    CladeDiversificationScanReport,
    DiversificationModelComparisonReport,
    DiversificationModelComparisonRow,
    LineageThroughTimeReport,
)


def format_figure_value(value: float) -> str:
    return format(value, ".15g")


def _lineage_color(point_count: int, point_index: int) -> str:
    if point_count <= 1:
        return "#0f766e"
    fraction = point_index / (point_count - 1)
    if fraction <= 0.25:
        return "#0f766e"
    if fraction <= 0.5:
        return "#0d9488"
    if fraction <= 0.75:
        return "#14b8a6"
    return "#2dd4bf"


def _classification_color(classification: str) -> str:
    if classification == "high":
        return "#b91c1c"
    if classification == "low":
        return "#1d4ed8"
    return "#94a3b8"


def _model_color(row: DiversificationModelComparisonRow, *, better_model: str) -> str:
    return "#ca8a04" if row.model == better_model else "#475569"


def write_ltt_svg(path: Path, report: LineageThroughTimeReport) -> int:
    width = 960
    height = 360
    left = 92
    right = 36
    top = 30
    bottom = 58
    plot_width = width - left - right
    plot_height = height - top - bottom
    maximum_time = max(report.root_age, 1e-9)
    maximum_lineages = max(point.lineage_count for point in report.points)
    time_ticks = [report.root_age, report.root_age / 2.0, 0.0]
    lineage_ticks = list(range(1, maximum_lineages + 1))

    def x_position(time_before_present: float) -> float:
        return left + ((maximum_time - time_before_present) / maximum_time) * plot_width

    def y_position(lineage_count: int) -> float:
        denominator = max(maximum_lineages - 1, 1)
        return top + ((maximum_lineages - lineage_count) / denominator) * plot_height

    segments: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="lineage through time curve">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        f'<text x="{left}" y="18" font-size="18" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Lineage-through-time curve</text>',
        f'<text x="{left}" y="{height - 18}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#475569">time before present</text>',
        f'<text x="18" y="{top - 8}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#475569">lineages</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#cbd5e1" stroke-width="1.5" />',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#cbd5e1" stroke-width="1.5" />',
    ]
    for tick in time_ticks:
        x = x_position(tick)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_height}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{x}" y="{top + plot_height + 22}" text-anchor="middle" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{escape(format_figure_value(tick))}</text>',
            ]
        )
    for tick in lineage_ticks:
        y = y_position(tick)
        segments.extend(
            [
                f'<line x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{left - 14}" y="{y + 4}" text-anchor="end" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{tick}</text>',
            ]
        )

    polyline = " ".join(
        f"{x_position(point.time_before_present):.2f},{y_position(point.lineage_count):.2f}"
        for point in report.points
    )
    segments.append(
        f'<polyline points="{polyline}" fill="none" stroke="#0f766e" stroke-width="4" stroke-linejoin="round" stroke-linecap="round" />'
    )
    for index, point in enumerate(report.points):
        x = x_position(point.time_before_present)
        y = y_position(point.lineage_count)
        color = _lineage_color(len(report.points), index)
        segments.extend(
            [
                f'<circle cx="{x}" cy="{y}" r="5.5" fill="{color}" stroke="#ffffff" stroke-width="1.5" />',
                f'<text x="{x}" y="{y - 10}" text-anchor="middle" font-size="11" font-family="Avenir Next, Segoe UI, sans-serif" fill="#134e4a">{escape(point.event)}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(report.points)


def _observation_label(observation: CladeDiversificationObservation) -> str:
    if observation.node_name:
        return observation.node_name
    if len(observation.descendant_taxa) <= 3:
        return ",".join(observation.descendant_taxa)
    return observation.node


def write_clade_outlier_svg(path: Path, report: CladeDiversificationScanReport) -> int:
    observations = sorted(
        report.observations,
        key=lambda row: (-row.z_score, row.tip_count, row.node),
    )
    row_height = 34
    width = 960
    height = max(240, 124 + len(observations) * row_height)
    left = 230
    right = 48
    top = 42
    plot_width = width - left - right
    plot_height = max(60, len(observations) * row_height)
    max_abs_z = max(max(abs(row.z_score) for row in observations), 1.0)

    def x_position(z_score: float) -> float:
        return left + ((z_score + max_abs_z) / (2.0 * max_abs_z)) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="clade diversification outliers">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fffaf5" />',
        f'<text x="{left}" y="20" font-size="18" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">Clade diversification outliers</text>',
        f'<text x="{left}" y="{height - 14}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">z-score relative to tree-wide diversification rate</text>',
    ]
    for tick in (-max_abs_z, -1.0, 0.0, 1.0, max_abs_z):
        x = x_position(tick)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_height}" stroke="#fde68a" stroke-width="1" />',
                f'<text x="{x}" y="{top + plot_height + 20}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#92400e">{escape(format_figure_value(tick))}</text>',
            ]
        )
    zero = x_position(0.0)
    segments.append(
        f'<line x1="{zero}" y1="{top}" x2="{zero}" y2="{top + plot_height}" stroke="#7c2d12" stroke-width="1.75" />'
    )
    for index, row in enumerate(observations):
        y = top + index * row_height + 18
        x = x_position(row.z_score)
        bar_left = min(zero, x)
        bar_width = max(abs(x - zero), 6.0)
        color = _classification_color(row.classification)
        segments.extend(
            [
                f'<text x="{left - 10}" y="{y + 4}" text-anchor="end" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#334155">{escape(_observation_label(row))}</text>',
                f'<rect x="{bar_left}" y="{y - 9}" width="{bar_width}" height="18" rx="6" fill="{color}" opacity="0.9" />',
                f'<text x="{left + plot_width + 8}" y="{y + 4}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#7c2d12">rate={escape(format_figure_value(row.diversification_rate))}, tips={row.tip_count}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(observations)


def write_model_comparison_svg(
    path: Path, report: DiversificationModelComparisonReport
) -> int:
    rows = sorted(report.rows, key=lambda row: (row.aic, row.model))
    width = 960
    height = max(250, 130 + len(rows) * 72)
    left = 220
    right = 54
    top = 44
    bottom = 46
    plot_width = width - left - right
    best_aic = min(row.aic for row in rows)
    max_delta = max(max(row.aic - best_aic for row in rows), 2.0)

    def x_position(delta_aic: float) -> float:
        return left + (delta_aic / max_delta) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="diversification model comparison">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fefce8" />',
        f'<text x="{left}" y="20" font-size="18" font-family="Avenir Next, Segoe UI, sans-serif" fill="#713f12">Diversification model comparison</text>',
        f'<text x="{left}" y="{height - 14}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#713f12">delta AIC relative to the better-supported model</text>',
    ]
    for tick in (0.0, max_delta / 2.0, max_delta):
        x = x_position(tick)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{height - bottom}" stroke="#fde68a" stroke-width="1" />',
                f'<text x="{x}" y="{height - bottom + 22}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#92400e">{escape(format_figure_value(tick))}</text>',
            ]
        )
    for index, row in enumerate(rows):
        y = top + index * 72 + 24
        delta_aic = row.aic - best_aic
        x = x_position(delta_aic)
        color = _model_color(row, better_model=report.better_model)
        segments.extend(
            [
                f'<text x="{left - 12}" y="{y + 4}" text-anchor="end" font-size="13" font-family="Avenir Next, Segoe UI, sans-serif" fill="#44403c">{escape(row.model)}</text>',
                f'<line x1="{left}" y1="{y}" x2="{x}" y2="{y}" stroke="{color}" stroke-width="4" stroke-linecap="round" />',
                f'<circle cx="{x}" cy="{y}" r="8" fill="{color}" stroke="#ffffff" stroke-width="1.5" />',
                f'<text x="{left}" y="{y - 14}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#713f12">AIC={escape(format_figure_value(row.aic))}</text>',
                f'<text x="{left}" y="{y + 22}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#713f12">net rate={escape(format_figure_value(row.net_diversification_rate))}, relative extinction={escape(format_figure_value(row.relative_extinction))}</text>',
            ]
        )
        if row.model == report.better_model:
            segments.append(
                f'<text x="{x + 14}" y="{y + 4}" font-size="11" font-family="Avenir Next, Segoe UI, sans-serif" fill="#854d0e">better supported</text>'
            )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)
