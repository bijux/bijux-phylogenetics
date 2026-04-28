from __future__ import annotations

from dataclasses import dataclass
from html import escape
import math
from pathlib import Path

from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree


@dataclass(slots=True)
class TreeRenderResult:
    output_path: Path
    format: str
    layout: str
    tip_count: int
    label_count: int
    has_scale_bar: bool
    rendered_support_count: int
    rendered_categorical_trait_count: int
    rendered_continuous_trait_count: int
    rendered_metadata_strip_count: int
    rendered_heatmap_column_count: int
    collapsed_clade_count: int
    missing_metadata_labels: list[str]


@dataclass(frozen=True, slots=True)
class _Point:
    x: float
    y: float


def _count_render_leaves(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_count_render_leaves(child) for child in node.children)


def _max_depth(node: TreeNode, depth: int = 0) -> int:
    if node.is_leaf():
        return depth
    return max(_max_depth(child, depth + 1) for child in node.children)


def _max_distance(node: TreeNode, distance: float = 0.0) -> float:
    next_distance = distance + float(node.branch_length or 0.0)
    if node.is_leaf():
        return next_distance
    return max(_max_distance(child, next_distance) for child in node.children)


def _nice_scale_bar_length(max_distance: float) -> float:
    if max_distance <= 0:
        return 0.0
    exponent = math.floor(math.log10(max_distance))
    base = 10**exponent
    for factor in (1.0, 0.5, 0.2, 0.1):
        candidate = base * factor
        if candidate <= max_distance / 3:
            return candidate
    return base / 10


def _format_branch_value(value: float) -> str:
    return format(round(value, 15), ".15g")


def render_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    labels: dict[str, str] | None = None,
    layout: str = "cladogram",
) -> TreeRenderResult:
    """Render a deterministic SVG tree as a cladogram or phylogram."""
    if layout not in {"cladogram", "phylogram"}:
        raise ValueError(f"unsupported tree layout: {layout}")

    tree = _load_tree(tree_path)
    labels = labels or {}
    row_height = 56
    left_margin = 48
    right_margin = 320
    top_margin = 40
    bottom_margin = 72
    horizontal_step = 150
    scale_width = 520
    leaf_count = _count_render_leaves(tree.root)
    max_depth = max(_max_depth(tree.root), 1)
    max_distance = _max_distance(tree.root, 0.0) if layout == "phylogram" else 0.0

    if layout == "phylogram" and max_distance > 0:
        tree_width = scale_width
        width = left_margin + tree_width + right_margin
    else:
        tree_width = horizontal_step * (max_depth + 1)
        width = left_margin + tree_width + right_margin
    height = top_margin + bottom_margin + row_height * max(leaf_count, 1)

    lines: list[str] = []
    texts: list[str] = []
    missing_labels: list[str] = []
    next_leaf_index = 0

    def node_x(depth: int, distance: float) -> float:
        if layout == "phylogram" and max_distance > 0:
            return left_margin + (distance / max_distance) * tree_width
        return left_margin + depth * horizontal_step

    def visit(node: TreeNode, depth: int, distance: float) -> _Point:
        nonlocal next_leaf_index
        branch_distance = distance + float(node.branch_length or 0.0)
        x = node_x(depth, branch_distance if node is not tree.root else distance)
        if node.is_leaf():
            y = top_margin + next_leaf_index * row_height + row_height / 2
            next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            if node.name and node.name not in labels and labels:
                missing_labels.append(node.name)
            texts.append(
                f'<text x="{x + 18:.1f}" y="{y + 5:.1f}" class="tip-label">{escape(label)}</text>'
            )
            return _Point(x=x, y=y)

        child_points = [visit(child, depth + 1, branch_distance) for child in node.children]
        y = sum(point.y for point in child_points) / len(child_points)
        min_y = min(point.y for point in child_points)
        max_y = max(point.y for point in child_points)
        lines.append(
            f'<line x1="{x:.1f}" y1="{min_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" class="branch"/>'
        )
        for child_point in child_points:
            lines.append(
                f'<line x1="{x:.1f}" y1="{child_point.y:.1f}" x2="{child_point.x:.1f}" y2="{child_point.y:.1f}" class="branch"/>'
            )
        return _Point(x=x, y=y)

    visit(tree.root, 0, 0.0)

    scale_bar = ""
    has_scale_bar = layout == "phylogram" and max_distance > 0
    if has_scale_bar:
        scale_length = _nice_scale_bar_length(max_distance)
        scale_start = left_margin
        scale_end = left_margin + (scale_length / max_distance) * tree_width
        scale_y = height - 28
        scale_bar = (
            f'<line x1="{scale_start:.1f}" y1="{scale_y:.1f}" x2="{scale_end:.1f}" y2="{scale_y:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_start:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_start:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_end:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_end:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<text x="{(scale_start + scale_end) / 2:.1f}" y="{scale_y - 10:.1f}" class="scale-label">{escape(_format_branch_value(scale_length))}</text>'
        )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="phylogenetic tree">
  <style>
    .panel {{ fill: #f7fbfa; stroke: #d7e3e1; stroke-width: 1; rx: 18; ry: 18; }}
    .branch {{ stroke: #0f172a; stroke-width: 2.2; stroke-linecap: round; fill: none; }}
    .scale-bar {{ stroke: #0f172a; stroke-width: 2; stroke-linecap: round; }}
    .tip-label {{ fill: #0f172a; font: 16px "Avenir Next", "Segoe UI", sans-serif; }}
    .scale-label {{ fill: #334155; text-anchor: middle; font: 13px "Avenir Next", "Segoe UI", sans-serif; }}
  </style>
  <rect x="1" y="1" width="{width - 2}" height="{height - 2}" class="panel" />
  {''.join(lines)}
  {''.join(texts)}
  {scale_bar}
</svg>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return TreeRenderResult(
        output_path=out_path,
        format="svg",
        layout=layout,
        tip_count=tree.tip_count,
        label_count=len(texts),
        has_scale_bar=has_scale_bar,
        rendered_support_count=0,
        rendered_categorical_trait_count=0,
        rendered_continuous_trait_count=0,
        rendered_metadata_strip_count=0,
        rendered_heatmap_column_count=0,
        collapsed_clade_count=0,
        missing_metadata_labels=sorted(set(missing_labels)),
    )
