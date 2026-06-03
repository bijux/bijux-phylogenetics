from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.alignment import (
    AlignmentWindowSummary,
    SequenceQualityRankingRow,
)

from .contracts import AlignmentFigureLegendEntry, AlignmentHeatmapCell


def heatmap_color(fraction: float) -> str:
    """Return the display color for one heatmap burden fraction."""
    clipped = max(0.0, min(fraction, 1.0))
    if clipped <= 0.05:
        return "#ecfeff"
    if clipped <= 0.15:
        return "#cffafe"
    if clipped <= 0.30:
        return "#67e8f9"
    if clipped <= 0.50:
        return "#22c55e"
    if clipped <= 0.70:
        return "#f59e0b"
    if clipped <= 0.90:
        return "#f97316"
    return "#b91c1c"


def score_color(score: float) -> str:
    """Return the display color for one sequence-quality score."""
    if score >= 90.0:
        return "#0f766e"
    if score >= 75.0:
        return "#1d4ed8"
    if score >= 60.0:
        return "#d97706"
    return "#b91c1c"


def region_classification(
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> dict[tuple[int, int], str]:
    """Classify windows into reviewer-facing over, under, or clear regions."""
    over = {(region.start, region.end) for region in over_regions}
    under = {(region.start, region.end) for region in under_regions}
    return {
        (window.start, window.end): (
            "over_aligned"
            if (window.start, window.end) in over
            else "under_aligned"
            if (window.start, window.end) in under
            else "clear"
        )
        for window in windows
    }


def line_path(points: list[tuple[float, float]]) -> str:
    """Serialize one SVG line path from chart points."""
    if not points:
        return ""
    return "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points)


def write_missingness_heatmap(
    path: Path,
    *,
    cells: list[AlignmentHeatmapCell],
    ranking_rows: list[SequenceQualityRankingRow],
    heatmap_bin_count: int,
) -> tuple[int, int]:
    """Write the reviewer-facing alignment missingness heatmap SVG."""
    row_count = len(ranking_rows)
    if heatmap_bin_count <= 0 or row_count <= 0:
        svg = "\n".join(
            [
                '<svg xmlns="http://www.w3.org/2000/svg" width="960" height="160" viewBox="0 0 960 160">',
                '<rect width="100%" height="100%" fill="#f8fafc"/>',
                '<text x="24" y="38" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">Alignment Missingness Heatmap</text>',
                '<text x="24" y="88" font-family="Avenir Next, Segoe UI, sans-serif" font-size="16" fill="#334155">No aligned sequence-by-site bins were available for heatmap review.</text>',
                "</svg>",
            ]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(svg + "\n", encoding="utf-8")
        return row_count, heatmap_bin_count

    cell_width = 10 if heatmap_bin_count <= 80 else 6
    cell_height = 22 if row_count <= 18 else 16
    left_margin = 190
    top_margin = 72
    width = left_margin + heatmap_bin_count * cell_width + 60
    height = top_margin + row_count * cell_height + 92
    ranked_ids = [row.identifier for row in ranking_rows]
    cell_lookup = {
        (cell.identifier, cell.bin_start, cell.bin_end): cell for cell in cells
    }
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">Alignment Missingness Heatmap</text>',
        '<text x="24" y="58" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#475569">Rows follow sequence-quality burden order; darker cells mean higher combined gap, missing, and ambiguity burden.</text>',
    ]
    bins = sorted({(cell.bin_start, cell.bin_end) for cell in cells})
    label_stride = max(1, len(bins) // 8)
    for row_index, identifier in enumerate(ranked_ids):
        y = top_margin + row_index * cell_height
        lines.append(
            f'<text x="24" y="{y + cell_height - 6}" font-family="SFMono-Regular, Consolas, monospace" font-size="12" fill="#1f2937">{escape(identifier)}</text>'
        )
        for column_index, site_bin in enumerate(bins):
            cell = cell_lookup[(identifier, site_bin[0], site_bin[1])]
            x = left_margin + column_index * cell_width
            lines.append(
                f'<rect x="{x}" y="{y}" width="{cell_width - 1}" height="{cell_height - 1}" fill="{heatmap_color(cell.uncertainty_fraction)}" stroke="#e2e8f0" stroke-width="0.5"/>'
            )
    for column_index, (start, end) in enumerate(bins):
        if column_index % label_stride != 0 and column_index != len(bins) - 1:
            continue
        x = left_margin + column_index * cell_width
        label = str(start) if start == end else f"{start}-{end}"
        lines.append(
            f'<text x="{x}" y="{top_margin - 10}" transform="rotate(-45 {x},{top_margin - 10})" font-family="SFMono-Regular, Consolas, monospace" font-size="10" fill="#475569">{escape(label)}</text>'
        )
    legend_y = height - 34
    legend_x = 24
    for offset, (label, color) in enumerate(
        [
            ("low", "#ecfeff"),
            ("mild", "#67e8f9"),
            ("moderate", "#22c55e"),
            ("high", "#f59e0b"),
            ("severe", "#b91c1c"),
        ]
    ):
        x = legend_x + offset * 120
        lines.extend(
            [
                f'<rect x="{x}" y="{legend_y}" width="18" height="18" rx="4" fill="{color}"/>',
                f'<text x="{x + 26}" y="{legend_y + 13}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#1f2937">{escape(label)}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return row_count, heatmap_bin_count


def write_site_quality_summary(
    path: Path,
    *,
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> int:
    """Write the reviewer-facing site-quality summary SVG."""
    width = 980
    height = 380
    chart_left = 74
    chart_top = 72
    chart_width = 860
    chart_height = 240
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">Alignment Site-Quality Summary</text>',
        '<text x="24" y="58" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#475569">Missingness, ambiguity, variability, and disagreement are summarized across sliding windows and suspicious windows are shaded directly on the panel.</text>',
    ]
    if not windows:
        lines.extend(
            [
                '<text x="24" y="108" font-family="Avenir Next, Segoe UI, sans-serif" font-size="16" fill="#334155">No alignment windows were available for plotting.</text>',
                "</svg>",
            ]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0

    lines.append(
        f'<rect x="{chart_left}" y="{chart_top}" width="{chart_width}" height="{chart_height}" fill="#ffffff" stroke="#cbd5e1" stroke-width="1"/>'
    )
    if len(windows) == 1:
        x_positions = [chart_left + chart_width / 2]
    else:
        x_positions = [
            chart_left + (chart_width * index / (len(windows) - 1))
            for index in range(len(windows))
        ]
    classification = region_classification(windows, over_regions, under_regions)
    for x, window in zip(x_positions, windows, strict=True):
        kind = classification[(window.start, window.end)]
        if kind == "clear":
            continue
        fill = "#fef3c7" if kind == "over_aligned" else "#fee2e2"
        lines.append(
            f'<rect x="{x - 18}" y="{chart_top}" width="36" height="{chart_height}" fill="{fill}" opacity="0.75"/>'
        )
    for tick in range(6):
        y = chart_top + chart_height * tick / 5
        value = 1.0 - tick / 5
        lines.extend(
            [
                f'<line x1="{chart_left}" y1="{y:.2f}" x2="{chart_left + chart_width}" y2="{y:.2f}" stroke="#e2e8f0" stroke-width="1"/>',
                f'<text x="24" y="{y + 4:.2f}" font-family="SFMono-Regular, Consolas, monospace" font-size="11" fill="#475569">{value:.1f}</text>',
            ]
        )

    def point(metric: float, x: float) -> tuple[float, float]:
        return (x, chart_top + (1.0 - metric) * chart_height)

    metric_specs = [
        ("missing", "#c2410c", [window.missing_fraction for window in windows]),
        ("ambiguity", "#b91c1c", [window.ambiguity_fraction for window in windows]),
        ("variable", "#0f766e", [window.variable_fraction for window in windows]),
        (
            "disagreement",
            "#1d4ed8",
            [window.disagreement_fraction for window in windows],
        ),
    ]
    for label, color, values in metric_specs:
        points = [point(value, x) for x, value in zip(x_positions, values, strict=True)]
        lines.append(
            f'<path d="{line_path(points)}" fill="none" stroke="{color}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}"/>')
        legend_y = chart_top + chart_height + 32
        legend_x = 28 + metric_specs.index((label, color, values)) * 170
        lines.extend(
            [
                f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 24}" y2="{legend_y}" stroke="{color}" stroke-width="3"/>',
                f'<text x="{legend_x + 32}" y="{legend_y + 4}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#1f2937">{escape(label)}</text>',
            ]
        )
    for x, window in zip(x_positions, windows, strict=True):
        label = (
            str(window.start)
            if window.start == window.end
            else f"{window.start}-{window.end}"
        )
        lines.append(
            f'<text x="{x:.2f}" y="{chart_top + chart_height + 14}" text-anchor="middle" font-family="SFMono-Regular, Consolas, monospace" font-size="10" fill="#475569">{escape(label)}</text>'
        )
    lines.extend(
        [
            '<rect x="712" y="324" width="14" height="14" fill="#fef3c7" opacity="0.75"/>',
            '<text x="734" y="335" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#1f2937">over-aligned window</text>',
            '<rect x="848" y="324" width="14" height="14" fill="#fee2e2" opacity="0.75"/>',
            '<text x="870" y="335" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#1f2937">under-aligned window</text>',
            "</svg>",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(windows)


def write_sequence_quality_panel(
    path: Path,
    *,
    ranking_rows: list[SequenceQualityRankingRow],
    maximum_rows: int,
) -> int:
    """Write the reviewer-facing sequence-quality burden SVG."""
    plotted_rows = ranking_rows[: max(maximum_rows, 1)]
    width = 980
    height = 120 + 46 * max(len(plotted_rows), 1)
    bar_left = 320
    bar_width = 520
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">Sequence-Quality Panel</text>',
        '<text x="24" y="58" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#475569">Lower-ranked sequences appear first so missingness, ambiguity, duplicate burden, and composition anomalies stay explicit on the figure instead of only in notes.</text>',
    ]
    if not plotted_rows:
        lines.extend(
            [
                '<text x="24" y="106" font-family="Avenir Next, Segoe UI, sans-serif" font-size="16" fill="#334155">No aligned sequences were available for ranking.</text>',
                "</svg>",
            ]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0
    for index, row in enumerate(plotted_rows):
        y = 76 + index * 46
        lines.extend(
            [
                f'<text x="24" y="{y + 15}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{escape(row.identifier)}</text>',
                f'<text x="24" y="{y + 31}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#475569">{escape(row.note)}</text>',
                f'<text x="250" y="{y + 16}" font-family="SFMono-Regular, Consolas, monospace" font-size="12" fill="#475569">rank {row.rank}</text>',
                f'<rect x="{bar_left}" y="{y}" width="{bar_width}" height="22" rx="8" fill="#e2e8f0"/>',
                f'<rect x="{bar_left}" y="{y}" width="{round(bar_width * row.score / 100.0, 3)}" height="22" rx="8" fill="{score_color(row.score)}"/>',
                f'<text x="{bar_left + bar_width + 12}" y="{y + 16}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{row.score:.3f}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(plotted_rows)


def write_heatmap_table(path: Path, cells: list[AlignmentHeatmapCell]) -> Path:
    """Write the heatmap support table artifact."""
    return write_taxon_rows(
        path,
        columns=[
            "identifier",
            "bin_start",
            "bin_end",
            "uncertainty_fraction",
            "gap_fraction",
            "missing_fraction",
            "ambiguity_fraction",
        ],
        rows=[
            {
                "identifier": cell.identifier,
                "bin_start": cell.bin_start,
                "bin_end": cell.bin_end,
                "uncertainty_fraction": format(cell.uncertainty_fraction, ".15g"),
                "gap_fraction": format(cell.gap_fraction, ".15g"),
                "missing_fraction": format(cell.missing_fraction, ".15g"),
                "ambiguity_fraction": format(cell.ambiguity_fraction, ".15g"),
            }
            for cell in cells
        ],
    )


def write_window_table(
    path: Path,
    *,
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> Path:
    """Write the sliding-window quality table artifact."""
    classification = region_classification(windows, over_regions, under_regions)
    return write_taxon_rows(
        path,
        columns=[
            "start",
            "end",
            "site_count",
            "gap_fraction",
            "missing_fraction",
            "ambiguity_fraction",
            "variable_fraction",
            "disagreement_fraction",
            "comparable_fraction",
            "classification",
        ],
        rows=[
            {
                "start": row.start,
                "end": row.end,
                "site_count": row.site_count,
                "gap_fraction": format(row.gap_fraction, ".15g"),
                "missing_fraction": format(row.missing_fraction, ".15g"),
                "ambiguity_fraction": format(row.ambiguity_fraction, ".15g"),
                "variable_fraction": format(row.variable_fraction, ".15g"),
                "disagreement_fraction": format(row.disagreement_fraction, ".15g"),
                "comparable_fraction": format(row.comparable_fraction, ".15g"),
                "classification": classification[(row.start, row.end)],
            }
            for row in windows
        ],
    )


def write_ranking_table(
    path: Path, ranking_rows: list[SequenceQualityRankingRow]
) -> Path:
    """Write the ranked sequence-quality table artifact."""
    return write_taxon_rows(
        path,
        columns=[
            "identifier",
            "rank",
            "score",
            "missing_fraction",
            "gap_fraction",
            "ambiguity_fraction",
            "composition_outlier",
            "duplicate_status",
            "note",
        ],
        rows=[
            {
                "identifier": row.identifier,
                "rank": row.rank,
                "score": format(row.score, ".15g"),
                "missing_fraction": format(row.missing_fraction, ".15g"),
                "gap_fraction": format(row.gap_fraction, ".15g"),
                "ambiguity_fraction": format(row.ambiguity_fraction, ".15g"),
                "composition_outlier": str(row.composition_outlier).lower(),
                "duplicate_status": row.duplicate_status,
                "note": row.note,
            }
            for row in ranking_rows
        ],
    )


def build_legend_entries() -> list[AlignmentFigureLegendEntry]:
    """Build the explicit figure legend entries for alignment review."""
    return [
        AlignmentFigureLegendEntry(
            surface="missingness-heatmap",
            label="sequence-by-site missingness burden",
            swatch="#67e8f9",
            detail="darker cells mean higher combined gap, explicit missing-data, and ambiguity burden within one site bin",
        ),
        AlignmentFigureLegendEntry(
            surface="site-quality-summary",
            label="windowed site-quality metrics",
            swatch="#1d4ed8",
            detail="the site summary keeps missingness, ambiguity, variability, and disagreement visible across sliding windows and highlights suspicious regions directly on the figure",
        ),
        AlignmentFigureLegendEntry(
            surface="sequence-quality-panel",
            label="reviewer-ranked sequence burden",
            swatch="#0f766e",
            detail="lower scores indicate higher missingness, ambiguity, duplicate burden, or composition risk",
        ),
    ]


def write_legend_table(path: Path, entries: list[AlignmentFigureLegendEntry]) -> Path:
    """Write the legend artifact for the alignment figure package."""
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
