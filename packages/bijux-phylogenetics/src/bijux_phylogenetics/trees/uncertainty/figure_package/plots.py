from __future__ import annotations

from html import escape
from pathlib import Path


def write_horizontal_bar_plot(
    path: Path,
    *,
    title: str,
    rows: list[tuple[str, float, str]],
    empty_message: str,
    fill_color: str,
    maximum_rows: int,
) -> int:
    """Write one reviewer-facing horizontal bar plot SVG artifact."""
    visible_rows = rows[: max(maximum_rows, 1)]
    width = 920
    height = 110 + 44 * max(len(visible_rows), 1)
    bar_left = 320
    bar_width = 520
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        f'<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">{escape(title)}</text>',
    ]
    if not visible_rows:
        lines.extend(
            [
                f'<text x="24" y="86" font-family="Avenir Next, Segoe UI, sans-serif" font-size="16" fill="#1f2937">{escape(empty_message)}</text>',
                "</svg>",
            ]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0
    for index, (label, fraction, note) in enumerate(visible_rows):
        y = 70 + index * 44
        filled_width = round(bar_width * max(0.0, min(fraction, 1.0)), 3)
        lines.extend(
            [
                f'<text x="24" y="{y + 16}" font-family="SFMono-Regular, Consolas, monospace" font-size="14" fill="#1f2937">{escape(label)}</text>',
                f'<text x="24" y="{y + 32}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#476b67">{escape(note)}</text>',
                f'<rect x="{bar_left}" y="{y}" width="{bar_width}" height="22" rx="8" fill="#e2e8f0"/>',
                f'<rect x="{bar_left}" y="{y}" width="{filled_width}" height="22" rx="8" fill="{fill_color}"/>',
                f'<text x="{bar_left + bar_width + 12}" y="{y + 17}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{format(fraction, ".3f")}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(visible_rows)


def write_clade_support_plot(
    path: Path,
    *,
    clade_frequencies,
    maximum_rows: int,
) -> int:
    """Write one clade-support plot for the uncertainty figure package."""
    return write_horizontal_bar_plot(
        path,
        title="Tree-Set Clade Support",
        rows=[
            (
                row.clade,
                row.frequency,
                f"{row.tree_count} / {clade_frequencies.tree_count} trees",
            )
            for row in sorted(
                clade_frequencies.clade_frequencies,
                key=lambda row: (-row.frequency, row.clade),
            )
        ],
        empty_message="No informative clades were available for support plotting.",
        fill_color="#0f766e",
        maximum_rows=maximum_rows,
    )


def write_unstable_taxa_plot(
    path: Path,
    *,
    unstable_taxa,
    maximum_rows: int,
) -> int:
    """Write one unstable-taxa plot for the uncertainty figure package."""
    return write_horizontal_bar_plot(
        path,
        title="Unstable Taxa Across Tree Set",
        rows=[
            (
                row.taxon,
                row.instability_score,
                f"{row.unique_placements} placements; dominant frequency {format(row.dominant_frequency, '.3f')}",
            )
            for row in unstable_taxa.taxa
        ],
        empty_message="No unstable taxa were detected; every shared taxon kept one placement signature across the tree set.",
        fill_color="#b45309",
        maximum_rows=maximum_rows,
    )


def write_topology_cluster_plot(
    path: Path,
    *,
    topology_clusters,
    maximum_rows: int,
) -> int:
    """Write one topology-cluster plot for the uncertainty figure package."""
    return write_horizontal_bar_plot(
        path,
        title="Topology Clusters Across Tree Set",
        rows=[
            (
                cluster.rooted_topology_id[:12],
                cluster.frequency,
                f"{cluster.tree_count} trees; representative index {cluster.representative_index}",
            )
            for cluster in topology_clusters.clusters
        ],
        empty_message="No topology clusters were available for plotting.",
        fill_color="#1d4ed8",
        maximum_rows=maximum_rows,
    )
