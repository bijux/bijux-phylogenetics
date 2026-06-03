from __future__ import annotations

from html import escape
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.render.tree_svg import render_tree_svg

from ..palette import _DEFAULT_STATE_COLORS
from .models import (
    StochasticMapBranchProbabilityRow,
    StochasticMapCollectionReport,
    StochasticMapDensityArtifactResult,
    StochasticMapDensityBranchRow,
    StochasticMapDensityReport,
    StochasticMapDensitySliceRow,
)
from .summary import summarize_discrete_stochastic_maps


def _binary_entropy(probability: float) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return float(
        format(
            -(probability * math.log2(probability))
            - ((1.0 - probability) * math.log2(1.0 - probability)),
            ".15g",
        )
    )


def _normalize_probability_interval(value: float, branch_length: float) -> float:
    if branch_length <= 0.0:
        return 0.0
    return float(format(value / branch_length, ".15g"))


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    token = color.removeprefix("#")
    return (
        int(token[0:2], 16),
        int(token[2:4], 16),
        int(token[4:6], 16),
    )


def _blend_density_color(
    probability: float,
    *,
    low_color: str = _DEFAULT_STATE_COLORS[0],
    high_color: str = _DEFAULT_STATE_COLORS[1],
) -> str:
    bounded = min(max(probability, 0.0), 1.0)
    low_red, low_green, low_blue = _hex_to_rgb(low_color)
    high_red, high_green, high_blue = _hex_to_rgb(high_color)
    red = round(low_red + (high_red - low_red) * bounded)
    green = round(low_green + (high_green - low_green) * bounded)
    blue = round(low_blue + (high_blue - low_blue) * bounded)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _tree_branch_geometry(
    report: StochasticMapCollectionReport,
) -> list[tuple[int, str, str, float, float, float]]:
    tree = load_tree(report.tree_path)
    root_candidates = {
        history.parent_node for history in report.maps[0].branch_histories
    } - {history.child_node for history in report.maps[0].branch_histories}
    if len(root_candidates) == 1:
        kept_taxa = next(iter(root_candidates)).split("|")
        if sorted(kept_taxa) != sorted(tree.tip_names):
            tree, _ = prune_tree_to_requested_taxa(report.tree_path, kept_taxa)
    expected_branch_indexes = {
        (
            history.parent_node,
            history.child_node,
        ): history.branch_index
        for history in report.maps[0].branch_histories
    }
    rows: list[tuple[int, str, str, float, float, float]] = []

    def visit(node, parent_depth: float) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            parent_node = node_signature(node)
            child_node = node_signature(child)
            branch_index = expected_branch_indexes.get((parent_node, child_node))
            if branch_index is not None:
                rows.append(
                    (
                        branch_index,
                        parent_node,
                        child_node,
                        branch_length,
                        float(format(parent_depth, ".15g")),
                        float(format(parent_depth + branch_length, ".15g")),
                    )
                )
            visit(child, parent_depth + branch_length)

    visit(tree.root, 0.0)
    return rows


def summarize_discrete_stochastic_map_density(
    report: StochasticMapCollectionReport,
    *,
    resolution: int = 100,
    focal_state: str | None = None,
) -> StochasticMapDensityReport:
    """Summarize one stochastic-map collection as branch probability density."""
    if resolution < 1:
        raise ValueError(f"resolution must be at least 1, got {resolution}")
    summary = summarize_discrete_stochastic_maps(report)
    branch_state_rows = [
        StochasticMapBranchProbabilityRow(
            branch_index=row.branch_index,
            parent_node=row.parent_node,
            child_node=row.child_node,
            state=row.state,
            branch_length=row.branch_length,
            mean_probability=_normalize_probability_interval(
                row.mean_time, row.branch_length
            ),
            lower_95_probability=_normalize_probability_interval(
                row.lower_95_interval, row.branch_length
            ),
            upper_95_probability=_normalize_probability_interval(
                row.upper_95_interval, row.branch_length
            ),
            minimum_probability=_normalize_probability_interval(
                row.minimum_time, row.branch_length
            ),
            maximum_probability=_normalize_probability_interval(
                row.maximum_time, row.branch_length
            ),
            presence_fraction=row.presence_fraction,
        )
        for row in summary.branch_occupancy_rows
    ]
    observed_state_order = list(dict.fromkeys(row.state for row in branch_state_rows))
    declared_state_order = list(report.fit_audit.state_order)
    if declared_state_order:
        state_order = [
            state for state in declared_state_order if state in observed_state_order
        ]
    else:
        state_order = []
    if not state_order:
        state_order = observed_state_order
    warnings = list(summary.warnings)
    resolved_focal_state = focal_state
    baseline_state: str | None = None
    if resolved_focal_state is None:
        if len(state_order) == 2:
            baseline_state = state_order[0]
            resolved_focal_state = state_order[1]
        elif len(state_order) > 2:
            warnings.append(
                "branch-state probability summaries are available for multistate collections, but density slices require one explicit focal state"
            )
    else:
        allowed_focal_states = declared_state_order or state_order
        if resolved_focal_state not in allowed_focal_states:
            raise ValueError(
                f"focal_state '{resolved_focal_state}' is not present in the stochastic-map state order"
            )
        if len(state_order) == 2:
            baseline_state = next(
                state for state in state_order if state != resolved_focal_state
            )
    branch_rows: list[StochasticMapDensityBranchRow] = []
    density_rows: list[StochasticMapDensitySliceRow] = []
    branch_geometry = _tree_branch_geometry(report)
    total_tree_depth = max(
        (end_depth for _, _, _, _, _, end_depth in branch_geometry),
        default=0.0,
    )
    if resolved_focal_state is None:
        return StochasticMapDensityReport(
            replicate_count=summary.replicate_count,
            resolution=resolution,
            total_tree_depth=total_tree_depth,
            state_order=state_order,
            focal_state=None,
            baseline_state=baseline_state,
            branch_state_rows=branch_state_rows,
            density_rows=density_rows,
            branch_rows=branch_rows,
            warnings=warnings,
        )
    branch_histories_by_replicate = [
        {history.branch_index: history for history in replicate.branch_histories}
        for replicate in report.maps
    ]
    steps = [
        (float(step) / float(resolution)) * total_tree_depth
        for step in range(resolution + 1)
    ]
    for (
        branch_index,
        parent_node,
        child_node,
        branch_length,
        start_depth,
        end_depth,
    ) in branch_geometry:
        boundaries = [start_depth]
        boundaries.extend(step for step in steps if start_depth < step < end_depth)
        boundaries.append(end_depth)
        branch_slice_rows: list[StochasticMapDensitySliceRow] = []
        if branch_length <= 0.0:
            branch_rows.append(
                StochasticMapDensityBranchRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    focal_state=resolved_focal_state,
                    baseline_state=baseline_state,
                    mean_posterior_probability=0.0,
                    minimum_posterior_probability=0.0,
                    maximum_posterior_probability=0.0,
                    uncertainty=0.0,
                    slice_count=0,
                )
            )
            continue
        for slice_index, (slice_start, slice_end) in enumerate(
            zip(boundaries, boundaries[1:], strict=False)
        ):
            slice_length = slice_end - slice_start
            if slice_length <= 0.0:
                continue
            slice_start_local = slice_start - start_depth
            slice_end_local = slice_end - start_depth
            replicate_probabilities: list[float] = []
            for branch_histories in branch_histories_by_replicate:
                history = branch_histories[branch_index]
                focal_duration = 0.0
                for segment in history.segments:
                    segment_start = segment.start_time_fraction * branch_length
                    segment_end = segment.end_time_fraction * branch_length
                    overlap = min(segment_end, slice_end_local) - max(
                        segment_start, slice_start_local
                    )
                    if overlap <= 0.0 or segment.state != resolved_focal_state:
                        continue
                    focal_duration += overlap
                replicate_probabilities.append(focal_duration / slice_length)
            posterior_probability = float(
                format(
                    sum(replicate_probabilities) / max(len(replicate_probabilities), 1),
                    ".15g",
                )
            )
            branch_slice_rows.append(
                StochasticMapDensitySliceRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    slice_index=slice_index,
                    start_depth=float(format(slice_start, ".15g")),
                    end_depth=float(format(slice_end, ".15g")),
                    start_time_fraction=float(
                        format(slice_start_local / branch_length, ".15g")
                    ),
                    end_time_fraction=float(
                        format(slice_end_local / branch_length, ".15g")
                    ),
                    posterior_probability=posterior_probability,
                    posterior_uncertainty=_binary_entropy(posterior_probability),
                )
            )
        density_rows.extend(branch_slice_rows)
        if not branch_slice_rows:
            branch_rows.append(
                StochasticMapDensityBranchRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    focal_state=resolved_focal_state,
                    baseline_state=baseline_state,
                    mean_posterior_probability=0.0,
                    minimum_posterior_probability=0.0,
                    maximum_posterior_probability=0.0,
                    uncertainty=0.0,
                    slice_count=0,
                )
            )
            continue
        weighted_probability_total = 0.0
        weighted_uncertainty_total = 0.0
        probability_values: list[float] = []
        for row in branch_slice_rows:
            slice_length = row.end_depth - row.start_depth
            weighted_probability_total += row.posterior_probability * slice_length
            weighted_uncertainty_total += row.posterior_uncertainty * slice_length
            probability_values.append(row.posterior_probability)
        branch_rows.append(
            StochasticMapDensityBranchRow(
                branch_index=branch_index,
                parent_node=parent_node,
                child_node=child_node,
                branch_length=branch_length,
                focal_state=resolved_focal_state,
                baseline_state=baseline_state,
                mean_posterior_probability=float(
                    format(weighted_probability_total / branch_length, ".15g")
                ),
                minimum_posterior_probability=min(probability_values, default=0.0),
                maximum_posterior_probability=max(probability_values, default=0.0),
                uncertainty=float(
                    format(weighted_uncertainty_total / branch_length, ".15g")
                ),
                slice_count=len(branch_slice_rows),
            )
        )
    return StochasticMapDensityReport(
        replicate_count=summary.replicate_count,
        resolution=resolution,
        total_tree_depth=total_tree_depth,
        state_order=state_order,
        focal_state=resolved_focal_state,
        baseline_state=baseline_state,
        branch_state_rows=branch_state_rows,
        density_rows=density_rows,
        branch_rows=branch_rows,
        warnings=warnings,
    )


def _write_stochastic_map_density_html(
    *,
    report: StochasticMapDensityReport,
    svg_path: Path,
    out_path: Path,
    layout: str,
) -> None:
    svg_markup = svg_path.read_text(encoding="utf-8")
    low_color = _DEFAULT_STATE_COLORS[0]
    high_color = _DEFAULT_STATE_COLORS[1]
    summary_cards = [
        ("focal state", report.focal_state or ""),
        ("baseline state", report.baseline_state or "complement"),
        ("resolution", str(report.resolution)),
        ("replicates", str(report.replicate_count)),
        ("branch density rows", str(len(report.branch_rows))),
    ]
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(f"Bijux Stochastic Density Map: {report.focal_state or 'state density'}")}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f8fafc;
      --panel: #ffffff;
      --rule: #d6dee8;
      --accent: #0f766e;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
    body {{
      margin: 0;
      padding: 2rem;
      background: linear-gradient(180deg, #eef6f4 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 1200px;
      margin: 0 auto;
      display: grid;
      gap: 1.5rem;
    }}
    section {{
      background: var(--panel);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1.5rem;
      box-shadow: 0 20px 60px rgba(15, 118, 110, 0.08);
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
    }}
    h1 {{
      margin: 0;
      color: var(--accent);
    }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 0.75rem;
    }}
    .summary-card {{
      background: #f6fffd;
      border: 1px solid rgba(15, 118, 110, 0.18);
      border-radius: 14px;
      padding: 0.9rem 1rem;
    }}
    .summary-card dt {{
      margin: 0;
      font: 600 0.82rem/1.2 "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      color: #476b67;
    }}
    .summary-card dd {{
      margin: 0.35rem 0 0;
      font: 700 1.1rem/1.2 var(--mono);
    }}
    .figure-shell svg {{
      width: 100%;
      height: auto;
      display: block;
    }}
    .legend-bar {{
      width: min(420px, 100%);
      height: 18px;
      border-radius: 999px;
      border: 1px solid rgba(15, 23, 42, 0.14);
      background: linear-gradient(90deg, {low_color} 0%, {high_color} 100%);
    }}
    .legend-labels {{
      display: flex;
      justify-content: space-between;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.9rem;
      margin-top: 0.35rem;
    }}
  </style>
</head>
<body>
  <main>
    <section>
      <h1>{escape(report.focal_state or "stochastic density map")}</h1>
      <p>Reviewer-facing stochastic-map density artifact with branch colors scaled by mean posterior probability.</p>
    </section>
    <section>
      <h2>Summary</h2>
      <div class="summary-grid">
        {"".join(f'<dl class="summary-card"><dt>{escape(label)}</dt><dd>{escape(value)}</dd></dl>' for label, value in summary_cards)}
      </div>
    </section>
    <section>
      <h2>Legend</h2>
      <div class="legend-bar"></div>
      <div class="legend-labels">
        <span>{escape(report.baseline_state or "0.0")}</span>
        <span>{escape(report.focal_state or "1.0")}</span>
      </div>
    </section>
    <section>
      <h2>Figure</h2>
      <div class="figure-shell">{svg_markup}</div>
      <p>Layout: <code>{escape(layout)}</code></p>
    </section>
  </main>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")


def render_stochastic_map_density_artifact(
    report: StochasticMapDensityReport,
    *,
    tree_path: Path,
    out_path: Path,
    layout: str = "phylogram",
) -> StochasticMapDensityArtifactResult:
    """Render one reviewer-facing stochastic-map density artifact."""
    if report.focal_state is None:
        raise ValueError(
            "stochastic-map density rendering requires one resolved focal state"
        )
    output_format = out_path.suffix.lower().lstrip(".")
    if output_format not in {"svg", "html"}:
        raise ValueError("stochastic-map density output must end in .svg or .html")
    branch_colors = {
        row.child_node: _blend_density_color(row.mean_posterior_probability)
        for row in report.branch_rows
    }
    svg_path = out_path if output_format == "svg" else out_path.with_suffix(".svg")
    render_result = render_tree_svg(
        tree_path,
        out_path=svg_path,
        layout=layout,
        branch_colors=branch_colors,
    )
    if output_format == "html":
        _write_stochastic_map_density_html(
            report=report,
            svg_path=svg_path,
            out_path=out_path,
            layout=layout,
        )
    return StochasticMapDensityArtifactResult(
        output_path=out_path,
        svg_path=svg_path,
        format=output_format,
        layout=layout,
        focal_state=report.focal_state,
        baseline_state=report.baseline_state,
        branch_count=len(report.branch_rows),
        rendered_branch_color_count=render_result.rendered_branch_color_count,
    )
