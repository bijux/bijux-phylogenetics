from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
    write_discrete_ancestral_probability_table,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    write_ancestral_state_table,
)
from bijux_phylogenetics.ancestral.presentation.visualization import (
    render_ancestral_state_visualization,
)
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)


@dataclass(slots=True)
class AncestralFigurePublicationAudit:
    """Publication-readiness audit for one ancestral-state figure package."""

    publication_ready: bool
    internal_state_visible: bool
    uncertainty_visible: bool
    legend_complete: bool
    caption_ready: bool
    internal_node_count: int
    rendered_internal_annotation_count: int
    rendered_internal_pie_count: int
    ambiguous_internal_node_count: int
    unstable_internal_node_count: int
    weak_support_node_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class AncestralFigurePackageResult:
    output_dir: Path
    figure_path: Path
    figure_png_path: Path
    figure_html_path: Path
    review_path: Path
    node_table_path: Path
    uncertainty_table_path: Path
    node_review_path: Path
    legend_path: Path
    model_description_path: Path
    caption_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    audit: AncestralFigurePublicationAudit
    reconstruction: ContinuousAncestralReport | DiscreteAncestralReport


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_ancestral_figure_package(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    out_dir: Path,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    layout: str = "phylogram",
) -> AncestralFigurePackageResult:
    """Build a publication-ready package for one ancestral-state reconstruction."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if reconstruction_kind == "continuous":
        report: ContinuousAncestralReport | DiscreteAncestralReport = (
            reconstruct_continuous_ancestral_states(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                alpha=alpha,
            )
        )
    else:
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )

    figure_path = out_dir / "ancestral-figure.svg"
    figure_png_path = out_dir / "ancestral-figure.png"
    figure_html_path = out_dir / "ancestral-figure.html"
    review_path = out_dir / "ancestral-figure-review.html"
    node_table_path = out_dir / "node-states.tsv"
    uncertainty_table_path = out_dir / "uncertainty.tsv"
    node_review_path = out_dir / "node-uncertainty-review.tsv"
    legend_path = out_dir / "legend.md"
    model_description_path = out_dir / "model-description.md"
    caption_path = out_dir / "figure-caption.md"
    manifest_path = out_dir / "figure-manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    figure_render = render_ancestral_state_visualization(
        tree_path,
        report,
        out_path=figure_path,
        layout=layout,
        discrete_node_style="pies",
        branch_coloring="regime"
        if isinstance(report, ContinuousAncestralReport)
        else "state",
        publication_annotations=True,
    )
    render_ancestral_state_visualization(
        tree_path,
        report,
        out_path=figure_png_path,
        layout=layout,
        discrete_node_style="pies",
        branch_coloring="regime"
        if isinstance(report, ContinuousAncestralReport)
        else "state",
        publication_annotations=True,
    )
    render_ancestral_state_visualization(
        tree_path,
        report,
        out_path=figure_html_path,
        layout=layout,
        discrete_node_style="pies",
        branch_coloring="regime"
        if isinstance(report, ContinuousAncestralReport)
        else "state",
        publication_annotations=True,
    )
    write_ancestral_state_table(node_table_path, report)
    if isinstance(report, ContinuousAncestralReport):
        write_continuous_ancestral_uncertainty_table(uncertainty_table_path, report)
    else:
        write_discrete_ancestral_probability_table(uncertainty_table_path, report)
    node_review_rows = _build_node_review_rows(report)
    _write_node_review_table(node_review_path, node_review_rows)
    legend_path.write_text(_legend_markdown(report), encoding="utf-8")
    model_description_path.write_text(_model_description(report), encoding="utf-8")
    audit = _build_publication_audit(
        report=report,
        render=figure_render,
        legend_text=legend_path.read_text(encoding="utf-8"),
    )
    caption_path.write_text(
        _caption_markdown(report, layout=layout, audit=audit),
        encoding="utf-8",
    )
    _write_publication_review_html(
        path=review_path,
        report=report,
        audit=audit,
        figure_svg_path=figure_path,
        node_review_rows=node_review_rows,
        uncertainty_table_path=uncertainty_table_path,
    )
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="ancestral_figure_package",
        input_files=[
            ("tree", tree_path),
            ("traits", traits_path),
        ],
        generated_figures=[
            ("figure_svg", figure_path),
            ("figure_png", figure_png_path),
            ("figure_html", figure_html_path),
        ],
        generated_tables=[
            ("node_states", node_table_path),
            ("uncertainty", uncertainty_table_path),
            ("node_review", node_review_path),
        ],
        filters=None,
        model={
            "kind": reconstruction_kind,
            "name": model,
            "alpha": alpha if reconstruction_kind == "continuous" else None,
            "state_ordering": (
                state_ordering if reconstruction_kind == "discrete" else None
            ),
            "ordered_states": (
                ordered_states if reconstruction_kind == "discrete" else None
            ),
        },
        settings={
            "trait": trait,
            "taxon_column": taxon_column,
            "layout": layout,
            "internal_node_count": audit.internal_node_count,
            "rendered_internal_annotation_count": (
                audit.rendered_internal_annotation_count
            ),
            "rendered_internal_pie_count": audit.rendered_internal_pie_count,
        },
        linked_artifacts=[
            ("legend", legend_path),
            ("model_description", model_description_path),
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    manifest_path.write_text(
        json.dumps(
            {
                "report_kind": "ancestral_figure_package",
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "reconstruction_kind": reconstruction_kind,
                "model": model,
                "layout": layout,
                "input_checksums": {
                    str(tree_path): _sha256(tree_path),
                    str(traits_path): _sha256(traits_path),
                },
                "artifacts": {
                    "figure": str(figure_path),
                    "figure_png": str(figure_png_path),
                    "figure_html": str(figure_html_path),
                    "review": str(review_path),
                    "node_table": str(node_table_path),
                    "uncertainty_table": str(uncertainty_table_path),
                    "node_review": str(node_review_path),
                    "legend": str(legend_path),
                    "model_description": str(model_description_path),
                    "caption": str(caption_path),
                    "reproducibility_manifest": str(reproducibility_manifest_path),
                },
                "artifact_checksums": {
                    str(path): _sha256(path)
                    for path in (
                        figure_path,
                        figure_png_path,
                        figure_html_path,
                        review_path,
                        node_table_path,
                        uncertainty_table_path,
                        node_review_path,
                        legend_path,
                        model_description_path,
                        caption_path,
                        reproducibility_manifest_path,
                    )
                },
                "reproducibility_manifest_path": str(reproducibility_manifest_path),
                "reproducibility_manifest_checksum": _sha256(
                    reproducibility_manifest_path
                ),
                "reproducibility_manifest": reproducibility_manifest,
                "audit": asdict(audit),
                "report": asdict(report),
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )
    return AncestralFigurePackageResult(
        output_dir=out_dir,
        figure_path=figure_path,
        figure_png_path=figure_png_path,
        figure_html_path=figure_html_path,
        review_path=review_path,
        node_table_path=node_table_path,
        uncertainty_table_path=uncertainty_table_path,
        node_review_path=node_review_path,
        legend_path=legend_path,
        model_description_path=model_description_path,
        caption_path=caption_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        audit=audit,
        reconstruction=report,
    )


def _build_node_review_rows(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for estimate in report.estimates:
        if estimate.is_tip:
            continue
        if isinstance(report, ContinuousAncestralReport):
            rows.append(
                {
                    "node": estimate.node,
                    "descendant_taxa": ",".join(estimate.descendant_taxa),
                    "state_label": format(estimate.estimate, ".3g"),
                    "uncertainty_label": f"+/-{format(estimate.standard_error, '.2g')}",
                    "interpretation": estimate.interpretation,
                    "ambiguous": "false",
                    "unstable": str(estimate.unstable).lower(),
                }
            )
            continue
        rows.append(
            {
                "node": estimate.node,
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "state_label": (
                    estimate.most_likely_state
                    if not estimate.ambiguous
                    else "/".join(estimate.state_set)
                ),
                "uncertainty_label": format(estimate.confidence, ".2f"),
                "interpretation": estimate.interpretation,
                "ambiguous": str(estimate.ambiguous).lower(),
                "unstable": str(estimate.unstable).lower(),
            }
        )
    return rows


def _write_node_review_table(path: Path, rows: list[dict[str, object]]) -> Path:
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "state_label",
            "uncertainty_label",
            "interpretation",
            "ambiguous",
            "unstable",
        ],
        rows=rows,
    )


def _build_publication_audit(
    *,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    render,
    legend_text: str,
) -> AncestralFigurePublicationAudit:
    internal_node_count = sum(1 for estimate in report.estimates if not estimate.is_tip)
    ambiguous_internal_node_count = (
        0
        if isinstance(report, ContinuousAncestralReport)
        else sum(
            1
            for estimate in report.estimates
            if not estimate.is_tip and estimate.ambiguous
        )
    )
    unstable_internal_node_count = len(report.unstable_nodes)
    weak_support_node_count = len(report.weak_support_nodes)
    legend_complete = (
        "uncertainty" in legend_text.lower()
        and "internal" in legend_text.lower()
        and (
            "pie" in legend_text.lower()
            or isinstance(report, ContinuousAncestralReport)
        )
    )
    caption_ready = True
    internal_state_visible = (
        render.tree_render.rendered_internal_annotation_count == internal_node_count
    )
    uncertainty_visible = (
        render.tree_render.rendered_internal_annotation_count == internal_node_count
        if isinstance(report, ContinuousAncestralReport)
        else (
            render.tree_render.rendered_internal_pie_count == internal_node_count
            and render.tree_render.rendered_internal_annotation_count
            == internal_node_count
        )
    )
    limitations = list(report.warnings)
    if not internal_state_visible:
        limitations.append(
            "one or more internal nodes are missing rendered ancestral-state labels"
        )
    if not uncertainty_visible:
        limitations.append(
            "one or more internal nodes are missing rendered uncertainty cues"
        )
    reviewer_summary = [
        f"rendered {internal_node_count} internal ancestral nodes with reviewer-facing labels",
        (
            "continuous node labels include estimates and standard errors"
            if isinstance(report, ContinuousAncestralReport)
            else "discrete node pies are paired with top-state confidence labels"
        ),
        f"unstable internal nodes: {unstable_internal_node_count}",
    ]
    return AncestralFigurePublicationAudit(
        publication_ready=internal_state_visible
        and uncertainty_visible
        and legend_complete
        and caption_ready,
        internal_state_visible=internal_state_visible,
        uncertainty_visible=uncertainty_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        internal_node_count=internal_node_count,
        rendered_internal_annotation_count=render.tree_render.rendered_internal_annotation_count,
        rendered_internal_pie_count=render.tree_render.rendered_internal_pie_count,
        ambiguous_internal_node_count=ambiguous_internal_node_count,
        unstable_internal_node_count=unstable_internal_node_count,
        weak_support_node_count=weak_support_node_count,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def _write_publication_review_html(
    *,
    path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    audit: AncestralFigurePublicationAudit,
    figure_svg_path: Path,
    node_review_rows: list[dict[str, object]],
    uncertainty_table_path: Path,
) -> Path:
    svg_markup = figure_svg_path.read_text(encoding="utf-8")
    table_rows = "".join(
        "<tr>"
        f"<td>{row['node']}</td>"
        f"<td>{row['state_label']}</td>"
        f"<td>{row['uncertainty_label']}</td>"
        f"<td>{row['interpretation']}</td>"
        "</tr>"
        for row in node_review_rows
    )
    path.write_text(
        "\n".join(
            [
                "<!doctype html>",
                '<html lang="en">',
                "<head>",
                '  <meta charset="utf-8">',
                '  <meta name="viewport" content="width=device-width, initial-scale=1">',
                f"  <title>Bijux Ancestral Figure Review: {report.trait}</title>",
                "  <style>",
                '    body { font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif; color: #172033; background: #f7f9fc; margin: 0; padding: 2rem; }',
                "    main { max-width: 1180px; margin: 0 auto; background: rgba(255,255,255,0.96); border: 1px solid #d7e0ea; border-radius: 24px; padding: 2rem; box-shadow: 0 28px 80px rgba(15,23,42,0.08); }",
                '    h1, h2 { font-family: "Avenir Next", "Segoe UI", sans-serif; }',
                "    table { width: 100%; border-collapse: collapse; margin-top: 1rem; }",
                "    th, td { border-bottom: 1px solid #d7e0ea; padding: 0.55rem; text-align: left; vertical-align: top; }",
                "    ul { margin: 0.5rem 0 0; }",
                '    .metric { display: inline-block; margin: 0.25rem 0.75rem 0.25rem 0; padding: 0.45rem 0.7rem; border-radius: 999px; background: #e8f0fb; font-family: "Avenir Next", "Segoe UI", sans-serif; }',
                "    .figure-shell { margin-top: 1.5rem; overflow-x: auto; }",
                "  </style>",
                "</head>",
                "<body>",
                "<main>",
                f"<h1>Bijux Ancestral Figure Review: {report.trait}</h1>",
                "<p>This publication review checks whether internal ancestral states and their uncertainty remain visible and interpretable directly on the figure, not only in supplementary tables.</p>",
                "<div>",
                f'<span class="metric">publication_ready={str(audit.publication_ready).lower()}</span>',
                f'<span class="metric">internal_state_visible={str(audit.internal_state_visible).lower()}</span>',
                f'<span class="metric">uncertainty_visible={str(audit.uncertainty_visible).lower()}</span>',
                f'<span class="metric">ambiguous_internal_nodes={audit.ambiguous_internal_node_count}</span>',
                f'<span class="metric">unstable_internal_nodes={audit.unstable_internal_node_count}</span>',
                "</div>",
                "<h2>Reviewer Summary</h2>",
                "<ul>",
                *[f"<li>{line}</li>" for line in audit.reviewer_summary],
                "</ul>",
                "<h2>Publication Limitations</h2>",
                (
                    "<p>none</p>"
                    if not audit.limitations
                    else "<ul>"
                    + "".join(f"<li>{line}</li>" for line in audit.limitations)
                    + "</ul>"
                ),
                '<div class="figure-shell">',
                svg_markup,
                "</div>",
                "<h2>Node Review</h2>",
                "<table>",
                "<thead><tr><th>node</th><th>state label</th><th>uncertainty label</th><th>interpretation</th></tr></thead>",
                f"<tbody>{table_rows}</tbody>",
                "</table>",
                f"<p>Supplementary uncertainty ledger: <code>{uncertainty_table_path.name}</code></p>",
                "</main>",
                "</body>",
                "</html>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def _write_uncertainty_table(
    path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> Path:
    if isinstance(report, ContinuousAncestralReport):
        return write_ancestral_rows(
            path,
            columns=[
                "node",
                "descendant_taxa",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
                "interpretation",
            ],
            rows=[
                {
                    "node": estimate.node,
                    "descendant_taxa": ",".join(estimate.descendant_taxa),
                    "standard_error": str(estimate.standard_error),
                    "lower_95_interval": str(estimate.lower_95_interval),
                    "upper_95_interval": str(estimate.upper_95_interval),
                    "interpretation": estimate.interpretation,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
        )
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "most_likely_state",
            "confidence",
            "unstable",
            "interpretation",
        ],
        rows=[
            {
                "node": estimate.node,
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "most_likely_state": estimate.most_likely_state,
                "confidence": str(estimate.confidence),
                "unstable": str(estimate.unstable).lower(),
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def _legend_markdown(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> str:
    if isinstance(report, ContinuousAncestralReport):
        return (
            "# Legend\n\n"
            "- Internal labels: reconstructed continuous ancestral values followed by `+/-` one standard error.\n"
            "- Uncertainty visibility: each internal node label carries the same uncertainty surface shown in the supplementary uncertainty ledger.\n"
            "- Uncertainty table: standard error and 95% interval for each internal node.\n"
        )
    states = "\n".join(f"- `{state}`" for state in report.observed_states)
    return (
        "# Legend\n\n"
        "- Internal labels: reconstructed discrete ancestral states with top-state confidence shown as a probability.\n"
        "- Internal pies: posterior or parsimony-supported state mixtures at each internal node.\n"
        "- Observed states:\n"
        f"{states}\n"
        "- Uncertainty table: most-likely state confidence and instability flag per internal node.\n"
    )


def _model_description(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> str:
    if isinstance(report, ContinuousAncestralReport):
        return (
            "# Model Description\n\n"
            f"- Kind: continuous ancestral reconstruction\n"
            f"- Model: `{report.model}`\n"
            f"- Alpha: `{report.alpha}`\n"
            f"- Trait: `{report.trait}`\n"
        )
    return (
        "# Model Description\n\n"
        f"- Kind: discrete ancestral reconstruction\n"
        f"- Model: `{report.model}`\n"
        f"- Trait: `{report.trait}`\n"
        f"- Observed states: `{', '.join(report.observed_states)}`\n"
    )


def _caption_markdown(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
    *,
    layout: str,
    audit: AncestralFigurePublicationAudit,
) -> str:
    return (
        f"# Ancestral State Figure: {report.trait}\n\n"
        f"- Layout: `{layout}`\n"
        f"- Model: `{report.model}`\n"
        f"- Taxa analyzed: `{report.taxon_count}`\n"
        f"- Internal nodes summarized: `{sum(1 for estimate in report.estimates if not estimate.is_tip)}`\n"
        + (
            "- Figure uncertainty surface: internal labels show estimate plus standard error at each node.\n"
            if isinstance(report, ContinuousAncestralReport)
            else "- Figure uncertainty surface: internal pies show state mixtures and internal labels report top-state confidence.\n"
        )
        + f"- Publication ready: `{audit.publication_ready}`\n"
    )
