from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.core.alignment import (
    AlignmentForensicReport,
    AlignmentRecord,
    AlignmentSummary,
    AlignmentWindowSummary,
    SequenceQualityRankingRow,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_forensic_report,
    summarize_alignment_windows,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.reports.review import (
    build_reviewer_audit_checklist,
    ReviewerAuditChecklist,
    write_reviewer_audit_checklist,
)


@dataclass(frozen=True, slots=True)
class AlignmentHeatmapCell:
    """One sequence-by-site-bin uncertainty cell for the alignment heatmap."""

    identifier: str
    bin_start: int
    bin_end: int
    uncertainty_fraction: float
    gap_fraction: float
    missing_fraction: float
    ambiguity_fraction: float


@dataclass(frozen=True, slots=True)
class AlignmentFigureLegendEntry:
    """One explicit legend entry for the alignment figure package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class AlignmentFigureCaptionDraft:
    """Structured caption draft for the alignment figure package."""

    title: str
    lead_sentence: str
    heatmap_sentence: str
    site_summary_sentence: str
    sequence_panel_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class AlignmentFigureAudit:
    """Publication-oriented audit for alignment figure completeness and risk."""

    publication_ready: bool
    heatmap_visible: bool
    site_summary_visible: bool
    sequence_panel_visible: bool
    legend_complete: bool
    caption_ready: bool
    suspicious_alignment: bool
    quality_score: float
    heatmap_row_count: int
    heatmap_bin_count: int
    plotted_window_count: int
    plotted_sequence_count: int
    invalid_character_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class AlignmentFigurePackageResult:
    output_dir: Path
    heatmap_figure_path: Path
    site_summary_figure_path: Path
    sequence_panel_figure_path: Path
    heatmap_table_path: Path
    window_table_path: Path
    ranking_table_path: Path
    legend_path: Path
    caption_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    reviewer_audit_checklist_path: Path
    summary: AlignmentSummary
    forensic: AlignmentForensicReport
    windows: list[AlignmentWindowSummary]
    heatmap_cells: list[AlignmentHeatmapCell]
    legend_entries: list[AlignmentFigureLegendEntry]
    caption_draft: AlignmentFigureCaptionDraft
    audit: AlignmentFigureAudit
    reviewer_audit_checklist: ReviewerAuditChecklist
    machine_manifest: dict[str, object]


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _site_bins(alignment_length: int, *, maximum_bins: int) -> list[tuple[int, int]]:
    if alignment_length <= 0:
        return []
    if maximum_bins <= 0:
        raise ValueError(f"maximum_bins must be positive, got {maximum_bins}")
    if alignment_length <= maximum_bins:
        return [(position, position) for position in range(1, alignment_length + 1)]
    bin_width = -(-alignment_length // maximum_bins)
    bins: list[tuple[int, int]] = []
    start = 1
    while start <= alignment_length:
        end = min(start + bin_width - 1, alignment_length)
        bins.append((start, end))
        start = end + 1
    return bins


def _mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return round(sum(values) / len(values), 15)


def _classify_residue(
    residue: str,
    *,
    alphabet: str,
) -> tuple[float, float, float]:
    if residue == "-":
        return 1.0, 0.0, 0.0
    if residue == "?":
        return 0.0, 1.0, 0.0
    canonical = (
        {"A", "C", "G", "T"}
        if alphabet == "dna"
        else {"A", "C", "G", "U"}
        if alphabet == "rna"
        else set("ACDEFGHIKLMNPQRSTVWY")
    )
    return (0.0, 0.0, 0.0) if residue.upper() in canonical else (0.0, 0.0, 1.0)


def _build_heatmap_cells(
    summary: AlignmentSummary,
    records: list[AlignmentRecord],
    ranking_rows: list[SequenceQualityRankingRow],
    *,
    maximum_bins: int,
) -> tuple[list[AlignmentHeatmapCell], int, int]:
    bins = _site_bins(summary.alignment_length, maximum_bins=maximum_bins)
    records_by_id = {record.identifier: record for record in records}
    cells: list[AlignmentHeatmapCell] = []
    for row in ranking_rows:
        record = records_by_id[row.identifier]
        for start, end in bins:
            gap_values: list[float] = []
            missing_values: list[float] = []
            ambiguity_values: list[float] = []
            for position in range(start - 1, end):
                gap_fraction, missing_fraction, ambiguity_fraction = _classify_residue(
                    record.sequence[position],
                    alphabet=summary.inferred_alphabet,
                )
                gap_values.append(gap_fraction)
                missing_values.append(missing_fraction)
                ambiguity_values.append(ambiguity_fraction)
            cells.append(
                AlignmentHeatmapCell(
                    identifier=row.identifier,
                    bin_start=start,
                    bin_end=end,
                    uncertainty_fraction=_mean(
                        [
                            gap_fraction + missing_fraction + ambiguity_fraction
                            for gap_fraction, missing_fraction, ambiguity_fraction in zip(
                                gap_values,
                                missing_values,
                                ambiguity_values,
                                strict=True,
                            )
                        ]
                    ),
                    gap_fraction=_mean(gap_values),
                    missing_fraction=_mean(missing_values),
                    ambiguity_fraction=_mean(ambiguity_values),
                )
            )
        if not bins:
            cells.append(
                AlignmentHeatmapCell(
                    identifier=row.identifier,
                    bin_start=1,
                    bin_end=0,
                    uncertainty_fraction=0.0,
                    gap_fraction=0.0,
                    missing_fraction=0.0,
                    ambiguity_fraction=0.0,
                )
            )
    return cells, len(ranking_rows), len(bins)


def _heatmap_color(fraction: float) -> str:
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


def _score_color(score: float) -> str:
    if score >= 90.0:
        return "#0f766e"
    if score >= 75.0:
        return "#1d4ed8"
    if score >= 60.0:
        return "#d97706"
    return "#b91c1c"


def _write_missingness_heatmap(
    path: Path,
    *,
    cells: list[AlignmentHeatmapCell],
    ranking_rows: list[SequenceQualityRankingRow],
    heatmap_bin_count: int,
) -> tuple[int, int]:
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
                f'<rect x="{x}" y="{y}" width="{cell_width - 1}" height="{cell_height - 1}" fill="{_heatmap_color(cell.uncertainty_fraction)}" stroke="#e2e8f0" stroke-width="0.5"/>'
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


def _region_classification(
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> dict[tuple[int, int], str]:
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


def _line_path(points: list[tuple[float, float]]) -> str:
    if not points:
        return ""
    return "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points)


def _write_site_quality_summary(
    path: Path,
    *,
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> int:
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
    classification = _region_classification(windows, over_regions, under_regions)
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
        ("disagreement", "#1d4ed8", [window.disagreement_fraction for window in windows]),
    ]
    for label, color, values in metric_specs:
        points = [point(value, x) for x, value in zip(x_positions, values, strict=True)]
        lines.append(
            f'<path d="{_line_path(points)}" fill="none" stroke="{color}" stroke-width="2.5"/>'
        )
        for x, y in points:
            lines.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{color}"/>'
            )
        legend_y = chart_top + chart_height + 32
        legend_x = 28 + metric_specs.index((label, color, values)) * 170
        lines.extend(
            [
                f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 24}" y2="{legend_y}" stroke="{color}" stroke-width="3"/>',
                f'<text x="{legend_x + 32}" y="{legend_y + 4}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#1f2937">{escape(label)}</text>',
            ]
        )
    for x, window in zip(x_positions, windows, strict=True):
        label = str(window.start) if window.start == window.end else f"{window.start}-{window.end}"
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


def _write_sequence_quality_panel(
    path: Path,
    *,
    ranking_rows: list[SequenceQualityRankingRow],
    maximum_rows: int,
) -> int:
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
                f'<rect x="{bar_left}" y="{y}" width="{round(bar_width * row.score / 100.0, 3)}" height="22" rx="8" fill="{_score_color(row.score)}"/>',
                f'<text x="{bar_left + bar_width + 12}" y="{y + 16}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{row.score:.3f}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(plotted_rows)


def _write_heatmap_table(path: Path, cells: list[AlignmentHeatmapCell]) -> Path:
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


def _write_window_table(
    path: Path,
    *,
    windows: list[AlignmentWindowSummary],
    over_regions,
    under_regions,
) -> Path:
    classification = _region_classification(windows, over_regions, under_regions)
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


def _write_ranking_table(
    path: Path, ranking_rows: list[SequenceQualityRankingRow]
) -> Path:
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


def _build_legend_entries() -> list[AlignmentFigureLegendEntry]:
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


def _write_legend_table(
    path: Path, entries: list[AlignmentFigureLegendEntry]
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


def _build_audit(
    *,
    summary: AlignmentSummary,
    forensic: AlignmentForensicReport,
    heatmap_row_count: int,
    heatmap_bin_count: int,
    plotted_window_count: int,
    plotted_sequence_count: int,
    legend_entries: list[AlignmentFigureLegendEntry],
) -> AlignmentFigureAudit:
    heatmap_visible = heatmap_row_count > 0 and heatmap_bin_count > 0
    site_summary_visible = plotted_window_count > 0
    sequence_panel_visible = plotted_sequence_count > 0
    legend_complete = {entry.surface for entry in legend_entries} == {
        "missingness-heatmap",
        "site-quality-summary",
        "sequence-quality-panel",
    }
    caption_ready = heatmap_visible and site_summary_visible and sequence_panel_visible
    suspicious_alignment = forensic.quality.suspicious_alignment
    publication_ready = (
        caption_ready
        and legend_complete
        and forensic.quality.quality_score >= 75.0
        and not suspicious_alignment
        and not summary.invalid_characters
    )
    reviewer_summary = [
        f"alignment quality score: {forensic.quality.quality_score}",
        f"heatmap rows and bins: {heatmap_row_count} x {heatmap_bin_count}",
        f"window summaries plotted: {plotted_window_count}",
        f"sequence rows plotted: {plotted_sequence_count}",
    ]
    limitations: list[str] = []
    if suspicious_alignment:
        limitations.extend(forensic.quality.suspicious_reasons)
    if summary.invalid_characters:
        limitations.append(
            "alignment contains invalid characters for the inferred alphabet"
        )
    if forensic.quality.quality_score < 75.0:
        limitations.append("alignment quality score remains below the reviewer threshold")
    if not summary.near_duplicate_scan_performed:
        limitations.append(
            "near-duplicate sequence review was skipped because the alignment exceeded the governed pairwise scan threshold"
        )
    if not heatmap_visible:
        limitations.append("the package does not currently render the missingness heatmap")
    if not site_summary_visible:
        limitations.append("the package does not currently render the site-quality summary")
    if not sequence_panel_visible:
        limitations.append("the package does not currently render the sequence-quality panel")
    if not limitations:
        limitations.append(
            "the current package keeps the key alignment quality figures explicit enough for publication-oriented review"
        )
    return AlignmentFigureAudit(
        publication_ready=publication_ready,
        heatmap_visible=heatmap_visible,
        site_summary_visible=site_summary_visible,
        sequence_panel_visible=sequence_panel_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        suspicious_alignment=suspicious_alignment,
        quality_score=forensic.quality.quality_score,
        heatmap_row_count=heatmap_row_count,
        heatmap_bin_count=heatmap_bin_count,
        plotted_window_count=plotted_window_count,
        plotted_sequence_count=plotted_sequence_count,
        invalid_character_count=len(summary.invalid_characters),
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def _build_caption_draft(
    *,
    summary: AlignmentSummary,
    audit: AlignmentFigureAudit,
) -> AlignmentFigureCaptionDraft:
    return AlignmentFigureCaptionDraft(
        title="Alignment quality review across missingness, site windows, and sequence burden",
        lead_sentence=(
            f"This package summarizes one {summary.inferred_alphabet} alignment with {summary.sequence_count} sequences and {summary.alignment_length} aligned sites through three explicit reviewer figures rather than leaving quality evidence buried in tables alone."
        ),
        heatmap_sentence=(
            f"The missingness heatmap keeps {audit.heatmap_row_count} sequence rows and {audit.heatmap_bin_count} site bins visible so missing-data concentration can be reviewed directly."
        ),
        site_summary_sentence=(
            f"The site-quality summary renders {audit.plotted_window_count} sliding windows and highlights suspicious over- or under-aligned regions on the figure itself."
        ),
        sequence_panel_sentence=(
            f"The sequence-quality panel keeps {audit.plotted_sequence_count} ranked sequence burdens explicit, with lower scores reflecting missingness, ambiguity, duplicate burden, or composition risk."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )


def _write_caption(path: Path, draft: AlignmentFigureCaptionDraft) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.heatmap_sentence,
                draft.site_summary_sentence,
                draft.sequence_panel_sentence,
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
    heatmap_figure_path: Path,
    site_summary_figure_path: Path,
    sequence_panel_figure_path: Path,
    heatmap_table_path: Path,
    window_table_path: Path,
    ranking_table_path: Path,
    legend_path: Path,
    caption_path: Path,
    reviewer_audit_checklist_path: Path,
    audit: AlignmentFigureAudit,
    reviewer_audit_checklist: ReviewerAuditChecklist,
) -> str:
    figures = {
        "heatmap": heatmap_figure_path.read_text(encoding="utf-8"),
        "site_summary": site_summary_figure_path.read_text(encoding="utf-8"),
        "sequence_panel": sequence_panel_figure_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>"
        + escape(label)
        + "</th><td>"
        + escape(value)
        + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("quality_score", format(audit.quality_score, ".15g")),
            ("suspicious_alignment", str(audit.suspicious_alignment).lower()),
            ("heatmap_visible", str(audit.heatmap_visible).lower()),
            ("site_summary_visible", str(audit.site_summary_visible).lower()),
            ("sequence_panel_visible", str(audit.sequence_panel_visible).lower()),
        ]
    )
    limitation_items = "".join(
        f"<li>{escape(item)}</li>" for item in audit.limitations
    )
    checklist_rows = "".join(
        "<tr><td>"
        + escape(item.section)
        + "</td><td>"
        + escape(item.status)
        + "</td><td>"
        + escape(item.summary)
        + "</td><td>"
        + escape("; ".join(item.evidence))
        + "</td></tr>"
        for item in reviewer_audit_checklist.items
    )
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Alignment Quality Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #eef6f4 0%, #f8fafc 100%); color: #1b1f24; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1220px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { color: #0f766e; margin-top: 0; }",
            "    .grid { display: grid; grid-template-columns: 1fr; gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.84); border: 1px solid rgba(15,118,110,0.14); border-radius: 18px; padding: 18px; box-shadow: 0 18px 42px rgba(15,118,110,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(15,118,110,0.12); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #0f766e; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Alignment Quality Review</h1>",
            "  <p>Reviewer-facing alignment figure package with one missingness heatmap, one site-quality summary, and one sequence-quality panel, backed by explicit ledgers and publication-oriented audit fields.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Reviewer Audit Checklist</h2>",
            "    <table><thead><tr><th>section</th><th>status</th><th>summary</th><th>evidence</th></tr></thead><tbody>"
            + checklist_rows
            + "</tbody></table>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Missingness Heatmap</h2><div class="figure-shell">' + figures["heatmap"] + "</div></section>",
            '    <section class="panel"><h2>Site-Quality Summary</h2><div class="figure-shell">' + figures["site_summary"] + "</div></section>",
            '    <section class="panel"><h2>Sequence-Quality Panel</h2><div class="figure-shell">' + figures["sequence_panel"] + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(heatmap_table_path.name)}">{escape(heatmap_table_path.name)}</a></li>',
            f'      <li><a href="{escape(window_table_path.name)}">{escape(window_table_path.name)}</a></li>',
            f'      <li><a href="{escape(ranking_table_path.name)}">{escape(ranking_table_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            f'      <li><a href="{escape(reviewer_audit_checklist_path.name)}">{escape(reviewer_audit_checklist_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def build_alignment_figure_package(
    alignment_path: Path,
    *,
    out_dir: Path,
    maximum_site_bins: int = 120,
    window_size: int = 30,
    step_size: int = 10,
    panel_row_limit: int = 12,
) -> AlignmentFigurePackageResult:
    """Build a reviewer-facing alignment figure package with explicit quality figures."""
    out_dir.mkdir(parents=True, exist_ok=True)
    heatmap_figure_path = out_dir / "alignment-missingness-heatmap.svg"
    site_summary_figure_path = out_dir / "alignment-site-quality-summary.svg"
    sequence_panel_figure_path = out_dir / "alignment-sequence-quality-panel.svg"
    heatmap_table_path = out_dir / "alignment-missingness-heatmap.tsv"
    window_table_path = out_dir / "alignment-site-quality-windows.tsv"
    ranking_table_path = out_dir / "alignment-sequence-quality-ranking.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "alignment-quality-review.html"
    manifest_path = out_dir / "alignment-quality-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"

    summary = summarise_fasta(alignment_path)
    records = load_fasta_alignment(alignment_path)
    forensic = build_alignment_forensic_report(alignment_path)
    windows = summarize_alignment_windows(
        alignment_path,
        window_size=window_size,
        step_size=step_size,
    )
    heatmap_cells, heatmap_row_count, heatmap_bin_count = _build_heatmap_cells(
        summary,
        records,
        forensic.sequence_ranking.rows,
        maximum_bins=maximum_site_bins,
    )
    _write_heatmap_table(heatmap_table_path, heatmap_cells)
    _write_window_table(
        window_table_path,
        windows=windows,
        over_regions=forensic.over_aligned_regions,
        under_regions=forensic.under_aligned_regions,
    )
    _write_ranking_table(ranking_table_path, forensic.sequence_ranking.rows)
    heatmap_row_count, heatmap_bin_count = _write_missingness_heatmap(
        heatmap_figure_path,
        cells=heatmap_cells,
        ranking_rows=forensic.sequence_ranking.rows,
        heatmap_bin_count=heatmap_bin_count,
    )
    plotted_window_count = _write_site_quality_summary(
        site_summary_figure_path,
        windows=windows,
        over_regions=forensic.over_aligned_regions,
        under_regions=forensic.under_aligned_regions,
    )
    plotted_sequence_count = _write_sequence_quality_panel(
        sequence_panel_figure_path,
        ranking_rows=forensic.sequence_ranking.rows,
        maximum_rows=panel_row_limit,
    )
    legend_entries = _build_legend_entries()
    _write_legend_table(legend_path, legend_entries)
    audit = _build_audit(
        summary=summary,
        forensic=forensic,
        heatmap_row_count=heatmap_row_count,
        heatmap_bin_count=heatmap_bin_count,
        plotted_window_count=plotted_window_count,
        plotted_sequence_count=plotted_sequence_count,
        legend_entries=legend_entries,
    )
    caption_draft = _build_caption_draft(summary=summary, audit=audit)
    _write_caption(caption_path, caption_draft)
    artifact_paths = [
        heatmap_figure_path,
        site_summary_figure_path,
        sequence_panel_figure_path,
        heatmap_table_path,
        window_table_path,
        ranking_table_path,
        legend_path,
        caption_path,
        review_path,
    ]
    existing_artifact_paths = artifact_paths[:-1]
    pre_review_manifest = {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": _checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {
            str(path): _checksum(path) for path in existing_artifact_paths
        },
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": _json_ready(asdict(summary)),
        "alignment_quality": _json_ready(asdict(forensic.quality)),
        "alignment_readiness": _json_ready(asdict(forensic.readiness)),
        "alignment_low_information": _json_ready(asdict(forensic.low_information)),
        "audit": _json_ready(asdict(audit)),
    }
    reviewer_audit_checklist = build_reviewer_audit_checklist(pre_review_manifest)
    review_path.write_text(
        _build_review_html(
            heatmap_figure_path=heatmap_figure_path,
            site_summary_figure_path=site_summary_figure_path,
            sequence_panel_figure_path=sequence_panel_figure_path,
            heatmap_table_path=heatmap_table_path,
            window_table_path=window_table_path,
            ranking_table_path=ranking_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            reviewer_audit_checklist_path=reviewer_audit_checklist_path,
            audit=audit,
            reviewer_audit_checklist=reviewer_audit_checklist,
        ),
        encoding="utf-8",
    )
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="alignment_quality_figure_package",
        input_files=[("alignment", alignment_path)],
        generated_figures=[
            ("missingness_heatmap", heatmap_figure_path),
            ("site_quality_summary", site_summary_figure_path),
            ("sequence_quality_panel", sequence_panel_figure_path),
        ],
        generated_tables=[
            ("missingness_heatmap", heatmap_table_path),
            ("site_quality_windows", window_table_path),
            ("sequence_quality_ranking", ranking_table_path),
            ("legend", legend_path),
        ],
        filters=None,
        model={
            "kind": "alignment_quality",
            "name": "summary-and-forensic-review",
        },
        settings={
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
            "alignment_length": summary.alignment_length,
            "sequence_count": summary.sequence_count,
        },
        linked_artifacts=[
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": _checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(
            reproducibility_manifest_path
        ),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": _json_ready(asdict(summary)),
        "alignment_quality": _json_ready(asdict(forensic.quality)),
        "alignment_readiness": _json_ready(asdict(forensic.readiness)),
        "alignment_low_information": _json_ready(asdict(forensic.low_information)),
        "audit": _json_ready(asdict(audit)),
    }
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["output_paths"].append(str(reviewer_audit_checklist_path))
    machine_manifest["output_checksums"][str(reviewer_audit_checklist_path)] = _checksum(
        reviewer_audit_checklist_path
    )
    machine_manifest["reviewer_audit_checklist_path"] = str(
        reviewer_audit_checklist_path
    )
    machine_manifest["reviewer_audit_checklist"] = _json_ready(
        asdict(reviewer_audit_checklist)
    )
    review_path.write_text(
        _build_review_html(
            heatmap_figure_path=heatmap_figure_path,
            site_summary_figure_path=site_summary_figure_path,
            sequence_panel_figure_path=sequence_panel_figure_path,
            heatmap_table_path=heatmap_table_path,
            window_table_path=window_table_path,
            ranking_table_path=ranking_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            reviewer_audit_checklist_path=reviewer_audit_checklist_path,
            audit=audit,
            reviewer_audit_checklist=reviewer_audit_checklist,
        ),
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return AlignmentFigurePackageResult(
        output_dir=out_dir,
        heatmap_figure_path=heatmap_figure_path,
        site_summary_figure_path=site_summary_figure_path,
        sequence_panel_figure_path=sequence_panel_figure_path,
        heatmap_table_path=heatmap_table_path,
        window_table_path=window_table_path,
        ranking_table_path=ranking_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        summary=summary,
        forensic=forensic,
        windows=windows,
        heatmap_cells=heatmap_cells,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )
