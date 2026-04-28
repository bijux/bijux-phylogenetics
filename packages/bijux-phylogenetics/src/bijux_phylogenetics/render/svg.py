from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree


@dataclass(slots=True)
class TreeRenderResult:
    output_path: Path
    format: str
    tip_count: int
    label_count: int
    missing_metadata_labels: list[str]


def _count_leaves(node) -> int:
    if node.is_leaf():
        return 1
    return sum(_count_leaves(child) for child in node.children)


def _max_depth(node, depth: int = 0) -> int:
    if node.is_leaf():
        return depth
    return max(_max_depth(child, depth + 1) for child in node.children)


def render_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    labels: dict[str, str] | None = None,
) -> TreeRenderResult:
    """Render a deterministic SVG cladogram for a tree."""
    tree = _load_tree(tree_path)
    labels = labels or {}
    row_height = 56
    left_margin = 40
    top_margin = 32
    horizontal_step = 160
    leaf_count = _count_leaves(tree.root)
    max_depth = max(_max_depth(tree.root), 1)
    width = left_margin + horizontal_step * (max_depth + 1) + 260
    height = top_margin * 2 + row_height * max(leaf_count - 1, 1)

    lines: list[str] = []
    texts: list[str] = []
    missing_labels: list[str] = []
    next_leaf_index = 0

    def visit(node, depth: int) -> tuple[float, float]:
        nonlocal next_leaf_index
        x = left_margin + depth * horizontal_step
        if node.is_leaf():
            y = top_margin + next_leaf_index * row_height
            next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            if node.name and node.name not in labels and labels:
                missing_labels.append(node.name)
            texts.append(f'<text x="{x + 18:.1f}" y="{y + 5:.1f}" class="tip-label">{escape(label)}</text>')
            return x, y

        child_positions = [visit(child, depth + 1) for child in node.children]
        y = sum(position[1] for position in child_positions) / len(child_positions)
        min_y = min(position[1] for position in child_positions)
        max_y = max(position[1] for position in child_positions)
        lines.append(f'<line x1="{x:.1f}" y1="{min_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" class="branch"/>')
        for child_x, child_y in child_positions:
            lines.append(f'<line x1="{x:.1f}" y1="{child_y:.1f}" x2="{child_x:.1f}" y2="{child_y:.1f}" class="branch"/>')
        return x, y

    visit(tree.root, 0)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="phylogenetic tree">
  <style>
    .panel {{ fill: #f7fbfa; stroke: #d7e3e1; stroke-width: 1; rx: 18; ry: 18; }}
    .branch {{ stroke: #0f172a; stroke-width: 2.2; stroke-linecap: round; }}
    .tip-label {{ fill: #0f172a; font: 16px "Avenir Next", "Segoe UI", sans-serif; }}
  </style>
  <rect x="1" y="1" width="{width - 2}" height="{height - 2}" class="panel" />
  {''.join(lines)}
  {''.join(texts)}
</svg>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return TreeRenderResult(
        output_path=out_path,
        format="svg",
        tip_count=tree.tip_count,
        label_count=len(texts),
        missing_metadata_labels=sorted(set(missing_labels)),
    )
