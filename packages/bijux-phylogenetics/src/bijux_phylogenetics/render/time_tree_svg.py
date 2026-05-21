from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import TreeNode


@dataclass(frozen=True, slots=True)
class TimeTreeNodeInterval:
    """Posterior age summary rendered for one internal time-tree node."""

    clade: str
    node_kind: str
    mean_age: float
    median_age: float
    minimum_age: float
    maximum_age: float
    lower_95_credible_interval: float
    upper_95_credible_interval: float
    tree_count: int


@dataclass(slots=True)
class TimeTreeRenderResult:
    output_path: Path
    format: str
    tip_count: int
    internal_node_count: int
    rendered_interval_count: int
    rendered_age_label_count: int
    root_age: float
    axis_tick_count: int
    ultrametric: bool
    longest_tip_label_length: int
    warnings: list[str]


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _clade_signature(node: TreeNode) -> str:
    return "|".join(_descendant_taxa(node))


def _count_internal_nodes(node: TreeNode) -> int:
    if node.is_leaf():
        return 0
    return 1 + sum(_count_internal_nodes(child) for child in node.children)


def _count_leaves(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_count_leaves(child) for child in node.children)


def _tip_depths(node: TreeNode, current_depth: float) -> list[float]:
    if node.is_leaf():
        return [current_depth]
    depths: list[float] = []
    for child in node.children:
        depths.extend(
            _tip_depths(child, current_depth + float(child.branch_length or 0.0))
        )
    return depths


def _tick_values(root_age: float) -> list[float]:
    if root_age <= 0.0:
        return [0.0]
    fractions = (1.0, 0.75, 0.5, 0.25, 0.0)
    ticks = [round(root_age * fraction, 6) for fraction in fractions]
    ordered: list[float] = []
    for value in ticks:
        if value not in ordered:
            ordered.append(value)
    return ordered


def render_time_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    node_intervals: list[TimeTreeNodeInterval],
    labels: dict[str, str] | None = None,
    title: str = "Bijux Time Tree Figure",
) -> TimeTreeRenderResult:
    """Render one rooted time tree with node-age labels and 95% HPD intervals."""
    tree = load_tree(tree_path)
    labels = labels or {taxon: taxon for taxon in tree.tip_names}
    interval_by_clade = {row.clade: row for row in node_intervals}
    tip_depths = _tip_depths(tree.root, 0.0)
    root_age = max(tip_depths, default=0.0)
    ultrametric = max(tip_depths, default=0.0) - min(tip_depths, default=0.0) <= 1e-9

    row_height = 72
    left_margin = 80
    right_margin = 340
    top_margin = 48
    bottom_margin = 96
    tree_width = 620
    tip_count = _count_leaves(tree.root)
    width = left_margin + tree_width + right_margin
    height = top_margin + bottom_margin + row_height * max(tip_count, 1)
    next_leaf_index = 0
    lines: list[str] = []
    texts: list[str] = []
    overlays: list[str] = []

    def node_x(current_depth: float) -> float:
        if root_age <= 0.0:
            return left_margin
        return left_margin + (current_depth / root_age) * tree_width

    def visit(node: TreeNode, current_depth: float) -> tuple[float, float]:
        nonlocal next_leaf_index
        x = node_x(current_depth)
        if node.is_leaf():
            y = top_margin + next_leaf_index * row_height + row_height / 2
            next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            texts.append(
                f'<text x="{x + 16:.1f}" y="{y + 5:.1f}" class="tip-label">{escape(label)}</text>'
            )
            return x, y

        child_points = []
        for child in node.children:
            branch_depth = current_depth + float(child.branch_length or 0.0)
            child_points.append(visit(child, branch_depth))
        y = sum(point[1] for point in child_points) / len(child_points)
        min_y = min(point[1] for point in child_points)
        max_y = max(point[1] for point in child_points)
        lines.append(
            f'<line x1="{x:.1f}" y1="{min_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" class="branch"/>'
        )
        for child_x, child_y in child_points:
            lines.append(
                f'<line x1="{x:.1f}" y1="{child_y:.1f}" x2="{child_x:.1f}" y2="{child_y:.1f}" class="branch"/>'
            )

        interval = interval_by_clade.get(_clade_signature(node))
        if interval is not None:
            older_depth = max(0.0, root_age - interval.upper_95_credible_interval)
            younger_depth = max(0.0, root_age - interval.lower_95_credible_interval)
            interval_x_left = node_x(older_depth)
            interval_x_right = node_x(younger_depth)
            interval_y = y - 15
            overlays.extend(
                [
                    f'<line x1="{interval_x_left:.1f}" y1="{interval_y:.1f}" x2="{interval_x_right:.1f}" y2="{interval_y:.1f}" class="hpd-interval"/>',
                    f'<line x1="{interval_x_left:.1f}" y1="{interval_y - 6:.1f}" x2="{interval_x_left:.1f}" y2="{interval_y + 6:.1f}" class="hpd-whisker"/>',
                    f'<line x1="{interval_x_right:.1f}" y1="{interval_y - 6:.1f}" x2="{interval_x_right:.1f}" y2="{interval_y + 6:.1f}" class="hpd-whisker"/>',
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" class="node-age-dot"/>',
                ]
            )
            texts.append(
                f'<text x="{x + 8:.1f}" y="{y - 20:.1f}" class="node-age-label">{escape(format(interval.median_age, ".3g"))}</text>'
            )
        return x, y

    visit(tree.root, 0.0)

    tick_values = _tick_values(root_age)
    axis_y = height - 48
    axis_markup = [
        f'<line x1="{left_margin:.1f}" y1="{axis_y:.1f}" x2="{left_margin + tree_width:.1f}" y2="{axis_y:.1f}" class="axis"/>'
    ]
    for value in tick_values:
        tick_depth = max(0.0, root_age - value)
        tick_x = node_x(tick_depth)
        axis_markup.extend(
            [
                f'<line x1="{tick_x:.1f}" y1="{axis_y - 5:.1f}" x2="{tick_x:.1f}" y2="{axis_y + 5:.1f}" class="axis"/>',
                f'<text x="{tick_x:.1f}" y="{axis_y + 22:.1f}" text-anchor="middle" class="axis-label">{escape(format(value, ".3g"))}</text>',
            ]
        )
    axis_markup.append(
        f'<text x="{left_margin + tree_width / 2:.1f}" y="{axis_y + 46:.1f}" text-anchor="middle" class="axis-title">node age before present</text>'
    )

    warnings: list[str] = []
    longest_tip_label_length = max((len(label) for label in labels.values()), default=0)
    if longest_tip_label_length > 28:
        warnings.append(
            "one or more tip labels are long enough that the publication figure may require abbreviated labels or a separate legend"
        )
    if tip_count > 40:
        warnings.append(
            "tip density is high enough that the time-tree figure may require panel subdivision before journal use"
        )
    if not ultrametric:
        warnings.append(
            "tree is not ultrametric, so age-axis interpretation is unsafe until branch lengths are corrected"
        )
    if len(node_intervals) < _count_internal_nodes(tree.root):
        warnings.append(
            "one or more internal nodes are missing posterior age intervals"
        )

    svg = "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            "<style>"
            ".title{font:700 20px 'Avenir Next','Segoe UI',sans-serif;fill:#0f172a}"
            ".branch{stroke:#334155;stroke-width:2.2;fill:none}"
            ".tip-label{font:500 14px 'Iowan Old Style','Palatino Linotype',serif;fill:#0f172a}"
            ".node-age-label{font:600 12px 'SFMono-Regular','SF Mono',Consolas,monospace;fill:#0f766e}"
            ".node-age-dot{fill:#0f766e;stroke:#ffffff;stroke-width:1.4}"
            ".hpd-interval{stroke:#b45309;stroke-width:3.2;stroke-linecap:round}"
            ".hpd-whisker{stroke:#b45309;stroke-width:1.8}"
            ".axis{stroke:#475569;stroke-width:1.6}"
            ".axis-label{font:500 11px 'SFMono-Regular','SF Mono',Consolas,monospace;fill:#334155}"
            ".axis-title{font:600 12px 'Avenir Next','Segoe UI',sans-serif;fill:#334155;letter-spacing:0.03em;text-transform:uppercase}"
            "</style>",
            '<rect width="100%" height="100%" fill="#f8fafc" />',
            f'<text x="{left_margin:.1f}" y="28" class="title">{escape(title)}</text>',
            *lines,
            *overlays,
            *texts,
            *axis_markup,
            "</svg>",
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return TimeTreeRenderResult(
        output_path=out_path,
        format=out_path.suffix.lstrip(".") or "svg",
        tip_count=tip_count,
        internal_node_count=_count_internal_nodes(tree.root),
        rendered_interval_count=sum(
            1
            for node in tree.iter_internal_nodes(order="preorder")
            if _clade_signature(node) in interval_by_clade
        ),
        rendered_age_label_count=sum(
            1
            for node in tree.iter_internal_nodes(order="preorder")
            if _clade_signature(node) in interval_by_clade
        ),
        root_age=round(root_age, 15),
        axis_tick_count=len(tick_values),
        ultrametric=ultrametric,
        longest_tip_label_length=longest_tip_label_length,
        warnings=warnings,
    )
