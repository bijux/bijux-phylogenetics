from __future__ import annotations

import csv
from pathlib import Path

from .contracts import (
    CladeCompatibilityEdgeRow,
    CladeCompatibilityGraphReport,
    CladeCompatibilityNodeRow,
)
from .inventory import _analyze_tree_set, _require_exact_taxa, _TreeSetAnalysis
from .topology import _clades_conflict, _format_clade


def _compatibility_reason(
    left: frozenset[str],
    right: frozenset[str],
) -> tuple[str, str]:
    if _clades_conflict(left, right):
        return ("conflict", "overlap-without-containment")
    if left <= right or right <= left:
        return ("compatible", "nested")
    return ("compatible", "disjoint")


def _build_clade_compatibility_graph_report(
    analysis: _TreeSetAnalysis,
) -> CladeCompatibilityGraphReport:
    exact_taxa = _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    total_tree_count = len(analysis.trees)
    ordered_clades = sorted(counts, key=_format_clade)
    node_rows: list[CladeCompatibilityNodeRow] = []
    edge_rows: list[CladeCompatibilityEdgeRow] = []
    compatible_neighbor_counts = dict.fromkeys(ordered_clades, 0)
    conflict_neighbor_counts = dict.fromkeys(ordered_clades, 0)
    compatible_edge_count = 0
    conflict_edge_count = 0

    for left_index, left_clade in enumerate(ordered_clades):
        for right_clade in ordered_clades[left_index + 1 :]:
            relation, reason = _compatibility_reason(left_clade, right_clade)
            if relation == "compatible":
                compatible_edge_count += 1
                compatible_neighbor_counts[left_clade] += 1
                compatible_neighbor_counts[right_clade] += 1
            else:
                conflict_edge_count += 1
                conflict_neighbor_counts[left_clade] += 1
                conflict_neighbor_counts[right_clade] += 1
            edge_rows.append(
                CladeCompatibilityEdgeRow(
                    left_clade=_format_clade(left_clade),
                    right_clade=_format_clade(right_clade),
                    compatibility_relation=relation,
                    compatibility_reason=reason,
                    left_tree_count=counts[left_clade],
                    right_tree_count=counts[right_clade],
                    left_frequency=round(counts[left_clade] / total_tree_count, 15),
                    right_frequency=round(counts[right_clade] / total_tree_count, 15),
                )
            )

    for clade in ordered_clades:
        node_rows.append(
            CladeCompatibilityNodeRow(
                clade=_format_clade(clade),
                tree_count=counts[clade],
                frequency=round(counts[clade] / total_tree_count, 15),
                compatible_neighbor_count=compatible_neighbor_counts[clade],
                conflict_neighbor_count=conflict_neighbor_counts[clade],
            )
        )

    return CladeCompatibilityGraphReport(
        path=analysis.path,
        tree_count=total_tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        node_count=len(node_rows),
        edge_count=len(edge_rows),
        compatible_edge_count=compatible_edge_count,
        conflict_edge_count=conflict_edge_count,
        nodes=node_rows,
        edges=edge_rows,
    )


def compute_clade_compatibility_graph(path: Path) -> CladeCompatibilityGraphReport:
    """Build one exact clade compatibility graph from the informative clades in a tree set."""
    return _build_clade_compatibility_graph_report(_analyze_tree_set(path))


def write_clade_compatibility_node_table(
    path: Path,
    report: CladeCompatibilityGraphReport,
) -> Path:
    """Write one node row per informative clade in the compatibility graph."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "tree_count",
                "frequency",
                "compatible_neighbor_count",
                "conflict_neighbor_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.nodes:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "compatible_neighbor_count": row.compatible_neighbor_count,
                    "conflict_neighbor_count": row.conflict_neighbor_count,
                }
            )
    return path


def write_clade_compatibility_edge_table(
    path: Path,
    report: CladeCompatibilityGraphReport,
) -> Path:
    """Write one compatibility or conflict edge row per clade pair."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_clade",
                "right_clade",
                "compatibility_relation",
                "compatibility_reason",
                "left_tree_count",
                "right_tree_count",
                "left_frequency",
                "right_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.edges:
            writer.writerow(
                {
                    "left_clade": row.left_clade,
                    "right_clade": row.right_clade,
                    "compatibility_relation": row.compatibility_relation,
                    "compatibility_reason": row.compatibility_reason,
                    "left_tree_count": row.left_tree_count,
                    "right_tree_count": row.right_tree_count,
                    "left_frequency": format(row.left_frequency, ".15g"),
                    "right_frequency": format(row.right_frequency, ".15g"),
                }
            )
    return path


def write_clade_compatibility_graph_dot(
    path: Path,
    report: CladeCompatibilityGraphReport,
) -> Path:
    """Write one DOT graph with clade nodes and compatibility or conflict edges."""
    path.parent.mkdir(parents=True, exist_ok=True)
    node_ids = {
        row.clade: f"n{index}" for index, row in enumerate(report.nodes, start=1)
    }
    lines = ["graph clade_compatibility {"]
    lines.append("  graph [overlap=false, splines=true];")
    lines.append('  node [shape=box, style="rounded"];')
    for row in report.nodes:
        label = (
            f"{row.clade}\\n"
            f"frequency={format(row.frequency, '.15g')}\\n"
            f"compatible={row.compatible_neighbor_count}\\n"
            f"conflict={row.conflict_neighbor_count}"
        )
        lines.append(f'  {node_ids[row.clade]} [label="{label}"];')
    for row in report.edges:
        color = (
            "darkgreen" if row.compatibility_relation == "compatible" else "firebrick"
        )
        style = "solid" if row.compatibility_relation == "compatible" else "dashed"
        label = row.compatibility_reason
        lines.append(
            "  "
            f"{node_ids[row.left_clade]} -- {node_ids[row.right_clade]} "
            f'[color="{color}", style="{style}", label="{label}"];'
        )
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
