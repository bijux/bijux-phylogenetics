from __future__ import annotations

from dataclasses import dataclass
from html import escape
import json
from pathlib import Path
import shutil

# Controlled local renderer execution only.
import subprocess  # nosec B404
import tempfile

from bijux_phylogenetics.ancestral.common import reconstruction_manifest
from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.render.tree_svg import TreeRenderResult, render_tree_svg

try:
    import cairosvg  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional renderer may miss native Cairo.
    cairosvg = None

_CATEGORICAL_PALETTE = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


@dataclass(slots=True)
class AncestralVisualizationResult:
    """Reviewer-facing ancestral visualization artifact."""

    output_path: Path
    format: str
    svg_path: Path
    layout: str
    reconstruction_kind: str
    model: str
    discrete_node_style: str
    branch_coloring: str
    tree_render: TreeRenderResult


def render_ancestral_state_visualization(
    tree_path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    *,
    out_path: Path,
    layout: str = "phylogram",
    discrete_node_style: str = "labels",
    branch_coloring: str = "none",
    publication_annotations: bool = False,
) -> AncestralVisualizationResult:
    """Render a governed ancestral-state visualization as SVG, PNG, or HTML."""
    output_format = out_path.suffix.lower().lstrip(".")
    if output_format not in {"svg", "png", "html"}:
        raise ValueError(
            "ancestral visualization output must end in .svg, .png, or .html"
        )
    if isinstance(report, ContinuousAncestralReport):
        if branch_coloring not in {"none", "regime"}:
            raise ValueError(
                "continuous ancestral visualization supports branch coloring 'none' or 'regime'"
            )
    else:
        if discrete_node_style not in {"labels", "pies"}:
            raise ValueError(
                "discrete ancestral visualization supports node style 'labels' or 'pies'"
            )
        if branch_coloring not in {"none", "state"}:
            raise ValueError(
                "discrete ancestral visualization supports branch coloring 'none' or 'state'"
            )

    svg_path = out_path if output_format == "svg" else out_path.with_suffix(".svg")
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="bijux-ancestral-visualization-"
    ) as temp_dir:
        analysis_tree_path = Path(temp_dir) / "analysis-tree.nwk"
        analysis_tree_path.write_text(
            f"{report.analysis_tree_newick}\n", encoding="utf-8"
        )
        tree_render = _render_ancestral_tree_svg(
            analysis_tree_path=analysis_tree_path,
            report=report,
            svg_path=svg_path,
            layout=layout,
            discrete_node_style=discrete_node_style,
            branch_coloring=branch_coloring,
            publication_annotations=publication_annotations,
        )

    if output_format == "png":
        _convert_svg_to_png(svg_path, out_path)
    elif output_format == "html":
        _write_ancestral_visualization_html(
            tree_path=tree_path,
            report=report,
            out_path=out_path,
            svg_path=svg_path,
            layout=layout,
            discrete_node_style=discrete_node_style,
            branch_coloring=branch_coloring,
            publication_annotations=publication_annotations,
        )

    return AncestralVisualizationResult(
        output_path=out_path,
        format=output_format,
        svg_path=svg_path,
        layout=layout,
        reconstruction_kind=(
            "continuous"
            if isinstance(report, ContinuousAncestralReport)
            else "discrete"
        ),
        model=report.model,
        discrete_node_style=discrete_node_style,
        branch_coloring=branch_coloring,
        tree_render=tree_render,
    )


def _categorical_palette(states: list[str]) -> dict[str, str]:
    return {
        state: _CATEGORICAL_PALETTE[index % len(_CATEGORICAL_PALETTE)]
        for index, state in enumerate(sorted(states))
    }


def _continuous_color(value: float, minimum: float, maximum: float) -> str:
    fraction = 0.5 if maximum <= minimum else (value - minimum) / (maximum - minimum)
    red = round(219 + (15 - 219) * fraction)
    green = round(234 + (118 - 234) * fraction)
    blue = round(254 + (110 - 254) * fraction)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _render_ancestral_tree_svg(
    *,
    analysis_tree_path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    svg_path: Path,
    layout: str,
    discrete_node_style: str,
    branch_coloring: str,
    publication_annotations: bool,
) -> TreeRenderResult:
    if isinstance(report, ContinuousAncestralReport):
        return _render_continuous_ancestral_tree_svg(
            analysis_tree_path=analysis_tree_path,
            report=report,
            svg_path=svg_path,
            layout=layout,
            branch_coloring=branch_coloring,
            publication_annotations=publication_annotations,
        )

    return _render_discrete_ancestral_tree_svg(
        analysis_tree_path=analysis_tree_path,
        report=report,
        svg_path=svg_path,
        layout=layout,
        discrete_node_style=discrete_node_style,
        branch_coloring=branch_coloring,
        publication_annotations=publication_annotations,
    )


def _render_continuous_ancestral_tree_svg(
    *,
    analysis_tree_path: Path,
    report: ContinuousAncestralReport,
    svg_path: Path,
    layout: str,
    branch_coloring: str,
    publication_annotations: bool,
) -> TreeRenderResult:
    (
        continuous_traits,
        internal_annotations,
        internal_annotation_colors,
        branch_colors,
    ) = _build_continuous_render_data(
        report,
        branch_coloring=branch_coloring,
        publication_annotations=publication_annotations,
    )
    return render_tree_svg(
        analysis_tree_path,
        out_path=svg_path,
        layout=layout,
        continuous_traits=continuous_traits,
        internal_annotations=internal_annotations,
        internal_annotation_colors=internal_annotation_colors,
        branch_colors=branch_colors,
    )


def _build_continuous_render_data(
    report: ContinuousAncestralReport,
    *,
    branch_coloring: str,
    publication_annotations: bool,
) -> tuple[dict[str, float], dict[str, str], dict[str, str], dict[str, str]]:
    value_by_node = {estimate.node: estimate.estimate for estimate in report.estimates}
    minimum = min(value_by_node.values()) if value_by_node else 0.0
    maximum = max(value_by_node.values()) if value_by_node else 0.0
    continuous_traits = {
        estimate.node_name: estimate.estimate
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    internal_annotations = {
        estimate.node: (
            f"{format(estimate.estimate, '.3g')} +/-{format(estimate.standard_error, '.2g')}"
            if publication_annotations
            else format(estimate.estimate, ".3g")
        )
        for estimate in report.estimates
        if not estimate.is_tip
    }
    internal_annotation_colors = {
        estimate.node: _continuous_color(estimate.estimate, minimum, maximum)
        for estimate in report.estimates
        if not estimate.is_tip
    }
    branch_colors = (
        {
            node: _continuous_color(value, minimum, maximum)
            for node, value in value_by_node.items()
        }
        if branch_coloring == "regime"
        else {}
    )
    return (
        continuous_traits,
        internal_annotations,
        internal_annotation_colors,
        branch_colors,
    )


def _render_discrete_ancestral_tree_svg(
    *,
    analysis_tree_path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    svg_path: Path,
    layout: str,
    discrete_node_style: str,
    branch_coloring: str,
    publication_annotations: bool,
) -> TreeRenderResult:
    (
        categorical_traits,
        internal_annotations,
        internal_annotation_colors,
        internal_pies,
        internal_pie_colors,
        branch_colors,
    ) = _build_discrete_render_data(
        report,
        discrete_node_style=discrete_node_style,
        branch_coloring=branch_coloring,
        publication_annotations=publication_annotations,
    )
    return render_tree_svg(
        analysis_tree_path,
        out_path=svg_path,
        layout=layout,
        categorical_traits=categorical_traits,
        internal_annotations=internal_annotations,
        internal_annotation_colors=internal_annotation_colors,
        branch_colors=branch_colors,
        internal_pies=internal_pies,
        internal_pie_colors=internal_pie_colors,
    )


def _build_discrete_render_data(
    report: DiscreteAncestralReport,
    *,
    discrete_node_style: str,
    branch_coloring: str,
    publication_annotations: bool,
) -> tuple[
    dict[str, str],
    dict[str, str],
    dict[str, str],
    dict[str, dict[str, float]],
    dict[str, str],
    dict[str, str],
]:
    palette = _categorical_palette(report.observed_states)
    internal_annotations = {
        estimate.node: (
            (
                f"{estimate.most_likely_state} {format(estimate.confidence, '.2f')}"
                if not estimate.ambiguous
                else f"{'/'.join(estimate.state_set)} {format(estimate.confidence, '.2f')}"
            )
            if publication_annotations
            else (
                estimate.most_likely_state
                if not estimate.ambiguous
                else "/".join(estimate.state_set)
            )
        )
        for estimate in report.estimates
        if not estimate.is_tip
        and (publication_annotations or discrete_node_style == "labels")
    }
    internal_pies = (
        {
            estimate.node: {
                state: probability
                for state, probability in estimate.state_probabilities.items()
                if probability > 0.0
            }
            for estimate in report.estimates
            if not estimate.is_tip
        }
        if discrete_node_style == "pies"
        else {}
    )
    categorical_traits = {
        estimate.node_name: estimate.most_likely_state
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    internal_annotation_colors = {
        estimate.node: palette.get(estimate.most_likely_state, "#6d28d9")
        for estimate in report.estimates
        if not estimate.is_tip
    }
    branch_colors = (
        {
            estimate.node: palette.get(estimate.most_likely_state, "#6d28d9")
            for estimate in report.estimates
            if estimate.node
        }
        if branch_coloring == "state"
        else {}
    )
    return (
        categorical_traits,
        internal_annotations,
        internal_annotation_colors,
        internal_pies,
        palette,
        branch_colors,
    )


def _convert_svg_to_png(svg_path: Path, png_path: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    if cairosvg is not None:
        try:
            cairosvg.svg2png(url=str(svg_path), write_to=str(png_path))
        except Exception as exc:  # pragma: no cover - defensive error capture
            failures.append(f"cairosvg: {exc}")
        else:
            if png_path.exists():
                return
    candidates = (
        ["rsvg-convert", str(svg_path), "-o", str(png_path)],
        ["sips", "-s", "format", "png", str(svg_path), "--out", str(png_path)],
    )
    for command in candidates:
        executable = command[0]
        if shutil.which(executable) is None:
            continue
        # The executable is checked with shutil.which before invocation.
        completed = subprocess.run(  # nosec B603
            command,
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0 and png_path.exists():
            return
        failures.append(
            f"{executable}: {(completed.stderr or completed.stdout).strip()}"
        )
    failure_text = "; ".join(failure for failure in failures if failure)
    if not failure_text:
        failure_text = "no supported SVG-to-PNG converter is available"
    raise RuntimeError(f"failed to export ancestral visualization PNG: {failure_text}")


def _write_ancestral_visualization_html(
    *,
    tree_path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    out_path: Path,
    svg_path: Path,
    layout: str,
    discrete_node_style: str,
    branch_coloring: str,
    publication_annotations: bool,
) -> None:
    svg_markup = svg_path.read_text(encoding="utf-8")
    title = f"Bijux Ancestral Visualization: {report.trait}"
    manifest = reconstruction_manifest(
        report_kind="ancestral-visualization",
        title=title,
        tree_path=tree_path,
        traits_path=report.traits_path,
        trait=report.trait,
        model=report.model,
        rendered_tree=str(svg_path),
    )
    manifest["layout"] = layout
    manifest["branch_coloring"] = branch_coloring
    manifest["discrete_node_style"] = discrete_node_style
    manifest["publication_annotations"] = publication_annotations
    manifest_script = json.dumps(manifest, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    kind = "continuous" if isinstance(report, ContinuousAncestralReport) else "discrete"
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f6faf9;
      --panel: #ffffff;
      --rule: #d6dee8;
      --accent: #0f766e;
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
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin: 0;
      color: var(--accent);
    }}
    figure {{
      margin: 0;
      overflow-x: auto;
    }}
    dl {{
      display: grid;
      grid-template-columns: max-content 1fr;
      gap: 0.5rem 1rem;
      margin: 0;
    }}
    dt {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      color: #475569;
    }}
    dd {{
      margin: 0;
    }}
    code {{
      font-family: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
  </style>
</head>
<body>
  <main>
    <section>
      <h1>{escape(title)}</h1>
      <p>Reviewer-facing ancestral {escape(kind)} visualization with explicit node styling and branch coloring choices.</p>
      <script id="bijux-report-manifest" type="application/json">{manifest_script}</script>
    </section>
    <section>
      <h2>Figure</h2>
      <figure>{svg_markup}</figure>
    </section>
    <section>
      <h2>Visualization Contract</h2>
      <dl>
        <dt>Trait</dt><dd><code>{escape(report.trait)}</code></dd>
        <dt>Model</dt><dd><code>{escape(report.model)}</code></dd>
        <dt>Layout</dt><dd><code>{escape(layout)}</code></dd>
        <dt>Branch coloring</dt><dd><code>{escape(branch_coloring)}</code></dd>
        <dt>Discrete node style</dt><dd><code>{escape(discrete_node_style)}</code></dd>
        <dt>Analyzed taxa</dt><dd><code>{report.taxon_count}</code></dd>
        <dt>Warnings</dt><dd><code>{len(report.warnings)}</code></dd>
      </dl>
    </section>
  </main>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
