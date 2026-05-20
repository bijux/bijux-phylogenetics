from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
import math
from pathlib import Path

from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)

from ...continuous import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from .contracts import (
    ComparativeModelCriteriaRow,
    ComparativeModelFigureAudit,
    ComparativeModelFigureCaptionDraft,
    ComparativeModelFigureLegendEntry,
    ComparativeModelFigurePackageResult,
    ComparativeModelFitRow,
    ComparativeModelLikelihoodRow,
    ComparativeModelParameterRow,
)
from .summaries import (
    build_model_criteria_rows as _criteria_rows,
    build_model_figure_audit as _build_audit,
    build_model_figure_caption_draft as _build_caption_draft,
    build_model_figure_legend as _legend_entries,
    build_model_fit_rows as _fit_rows,
    build_model_likelihood_rows as _likelihood_rows,
    build_model_parameter_rows as _parameter_rows,
    format_model_figure_number as _format_number,
    model_figure_color as _model_color,
    parameter_label as _parameter_label,
)
from .outputs import (
    write_model_figure_caption as _write_caption,
    write_model_figure_legend_table as _write_legend_table,
    write_model_figure_table as _write_table,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))




def _write_information_criteria_svg(
    path: Path,
    rows: list[ComparativeModelCriteriaRow],
    *,
    selected_model: str,
) -> int:
    width = 1080
    row_height = 52
    top = 84
    height = top + len(rows) * row_height + 48
    columns = [
        ("model", 56),
        ("logL", 264),
        ("AIC", 430),
        ("AICc", 566),
        ("delta AIC", 710),
        ("delta AICc", 866),
        ("selected", 992),
    ]
    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model information criteria">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        '<text x="56" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Comparative model information criteria</text>',
        '<text x="56" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">AIC and AICc remain explicit so the selected model is reviewable from the package itself.</text>',
        '<rect x="40" y="64" width="1000" height="36" rx="14" fill="#e2e8f0" />',
    ]
    for label, x in columns:
        anchor = "middle" if label != "model" else "start"
        segments.append(
            f'<text x="{x}" y="87" text-anchor="{anchor}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#334155">{escape(label)}</text>'
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        fill = "#ecfeff" if row.model == selected_model else "#ffffff"
        stroke = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="40" y="{y - 24}" width="1000" height="40" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="56" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<text x="264" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.log_likelihood))}</text>',
                f'<text x="430" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.aic))}</text>',
                f'<text x="566" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.aicc))}</text>',
                f'<text x="710" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.delta_aic))}</text>',
                f'<text x="866" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#1e293b">{escape(_format_number(row.delta_aicc))}</text>',
                f'<text x="992" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{"yes" if row.selected else "no"}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_likelihood_svg(
    path: Path,
    rows: list[ComparativeModelLikelihoodRow],
    *,
    selected_model: str,
) -> int:
    width = 980
    height = max(250, 130 + len(rows) * 84)
    left = 240
    right = 52
    top = 52
    bottom = 48
    plot_width = width - left - right
    max_log_likelihood = max(row.log_likelihood for row in rows)
    min_log_likelihood = min(row.log_likelihood for row in rows)
    span = max(max_log_likelihood - min_log_likelihood, 1e-9)

    def x_position(value: float) -> float:
        return left + ((value - min_log_likelihood) / span) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model likelihood comparison">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        '<text x="48" y="26" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Comparative model likelihood comparison</text>',
        '<text x="48" y="50" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">Higher log-likelihood bars indicate better likelihood fit before the information-criterion penalty is applied.</text>',
        f'<line x1="{left}" y1="{height - bottom}" x2="{left + plot_width}" y2="{height - bottom}" stroke="#cbd5e1" stroke-width="1.5" />',
    ]
    for tick_value in (min_log_likelihood, (min_log_likelihood + max_log_likelihood) / 2.0, max_log_likelihood):
        x = x_position(tick_value)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{height - bottom}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{x}" y="{height - 12}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{escape(_format_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(sorted(rows, key=lambda item: item.log_likelihood)):
        y = top + index * 84 + 24
        bar_end = x_position(row.log_likelihood)
        color = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<text x="{left - 16}" y="{y + 4}" text-anchor="end" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">{escape(row.model)}</text>',
                f'<rect x="{left}" y="{y - 10}" width="{max(bar_end - left, 4):.2f}" height="20" rx="10" fill="{color}" opacity="0.9" />',
                f'<text x="{bar_end + 12}" y="{y + 4}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">logL={escape(_format_number(row.log_likelihood))}; gap={escape(_format_number(row.delta_log_likelihood))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_parameter_svg(
    path: Path,
    rows: list[ComparativeModelParameterRow],
    *,
    selected_model: str,
) -> int:
    width = 1120
    row_height = 74
    left = 310
    right = 64
    top = 72
    bottom = 48
    height = top + len(rows) * row_height + bottom
    plot_width = width - left - right
    minimum = min(min(row.lower_95, row.estimate) for row in rows)
    maximum = max(max(row.upper_95, row.estimate) for row in rows)
    if math.isclose(minimum, maximum, abs_tol=1e-12):
        minimum -= 1.0
        maximum += 1.0
    padding = max((maximum - minimum) * 0.08, 1e-6)
    minimum -= padding
    maximum += padding
    span = maximum - minimum

    def x_position(value: float) -> float:
        return left + ((value - minimum) / span) * plot_width

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model parameter intervals">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#faf5ff" />',
        '<text x="56" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#3b0764">Comparative model parameter summaries</text>',
        '<text x="56" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#6b21a8">Each parameter keeps its fitted estimate and interval visible so the model comparison stays interpretable beyond the winner call.</text>',
    ]
    for tick_value in (minimum, (minimum + maximum) / 2.0, maximum):
        x = x_position(tick_value)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top - 8}" x2="{x}" y2="{height - bottom}" stroke="#e9d5ff" stroke-width="1" />',
                f'<text x="{x}" y="{height - 14}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#6b21a8">{escape(_format_number(tick_value))}</text>',
            ]
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        color = _model_color(row.model, selected_model=selected_model)
        lower = x_position(row.lower_95)
        estimate = x_position(row.estimate)
        upper = x_position(row.upper_95)
        segments.extend(
            [
                f'<text x="56" y="{y + 6}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#3b0764">{escape(row.model)}</text>',
                f'<text x="56" y="{y + 28}" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#581c87">{escape(_parameter_label(row.parameter))}</text>',
                f'<text x="56" y="{y + 48}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#7e22ce">{escape(row.interval_method)}</text>',
                f'<line x1="{lower:.2f}" y1="{y + 22}" x2="{upper:.2f}" y2="{y + 22}" stroke="{color}" stroke-width="4" stroke-linecap="round" />',
                f'<line x1="{lower:.2f}" y1="{y + 12}" x2="{lower:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<line x1="{upper:.2f}" y1="{y + 12}" x2="{upper:.2f}" y2="{y + 32}" stroke="{color}" stroke-width="2" />',
                f'<circle cx="{estimate:.2f}" cy="{y + 22}" r="7" fill="{color}" stroke="#ffffff" stroke-width="1.5" />',
                f'<text x="{estimate:.2f}" y="{y + 54}" text-anchor="middle" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#334155">{escape(_format_number(row.estimate))}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _write_fit_summary_svg(
    path: Path,
    rows: list[ComparativeModelFitRow],
    *,
    selected_model: str,
) -> int:
    width = 1140
    row_height = 56
    top = 84
    height = top + len(rows) * row_height + 48
    columns = [
        ("model", 60),
        ("residual variance", 286),
        ("max |z|", 474),
        ("residual lambda", 644),
        ("outliers", 816),
        ("warnings", 934),
        ("status", 1042),
    ]
    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="comparative model fit summary">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fffaf0" />',
        '<text x="60" y="28" font-size="22" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">Comparative model fit summaries</text>',
        '<text x="60" y="52" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#9a3412">Residual scale, standardized outlier burden, and warning counts remain reviewable alongside the information criteria.</text>',
        '<rect x="44" y="64" width="1050" height="36" rx="14" fill="#fed7aa" />',
    ]
    for label, x in columns:
        anchor = "middle" if label != "model" else "start"
        segments.append(
            f'<text x="{x}" y="87" text-anchor="{anchor}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#7c2d12">{escape(label)}</text>'
        )
    for index, row in enumerate(rows):
        y = top + index * row_height
        fill = "#fffbeb" if row.model == selected_model else "#ffffff"
        stroke = _model_color(row.model, selected_model=selected_model)
        segments.extend(
            [
                f'<rect x="44" y="{y - 26}" width="1050" height="42" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.25" />',
                f'<text x="60" y="{y}" font-size="14" font-family="Avenir Next, Segoe UI, sans-serif" fill="#451a03">{escape(row.model)}</text>',
                f'<text x="286" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.residual_variance))}</text>',
                f'<text x="474" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.max_abs_standardized_residual))}</text>',
                f'<text x="644" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{escape(_format_number(row.phylogenetic_residual_lambda))}</text>',
                f'<text x="816" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.outlier_taxon_count}</text>',
                f'<text x="934" y="{y}" text-anchor="middle" font-size="13" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#451a03">{row.warning_count}</text>',
                f'<text x="1042" y="{y}" text-anchor="middle" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="{stroke}">{escape(row.convergence_status)}</text>',
            ]
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(rows)


def _build_review_html(
    *,
    criteria_figure_path: Path,
    likelihood_figure_path: Path,
    parameter_figure_path: Path,
    fit_figure_path: Path,
    criteria_table_path: Path,
    likelihood_table_path: Path,
    parameter_table_path: Path,
    fit_table_path: Path,
    legend_path: Path,
    caption_path: Path,
    audit: ComparativeModelFigureAudit,
) -> str:
    figures = {
        "criteria": criteria_figure_path.read_text(encoding="utf-8"),
        "likelihood": likelihood_figure_path.read_text(encoding="utf-8"),
        "parameters": parameter_figure_path.read_text(encoding="utf-8"),
        "fit": fit_figure_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>"
        + escape(label)
        + "</th><td>"
        + escape(value)
        + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("selected_model", audit.selected_model),
            ("support_distinct", str(audit.support_distinct).lower()),
            ("aicc_delta", _format_number(audit.aicc_delta)),
            ("finite_aicc_model_count", str(audit.finite_aicc_model_count)),
            ("warning_count", str(audit.warning_count)),
        ]
    )
    limitation_items = "".join(f"<li>{escape(item)}</li>" for item in audit.limitations)
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Comparative Model Comparison Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #f8fafc 0%, #fffaf0 100%); color: #0f172a; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1280px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { margin-top: 0; color: #0f172a; }",
            "    .grid { display: grid; grid-template-columns: 1fr; gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.92); border: 1px solid rgba(15,23,42,0.08); border-radius: 18px; padding: 18px; box-shadow: 0 18px 40px rgba(15,23,42,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(148,163,184,0.35); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #0f766e; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Comparative Model Comparison Review</h1>",
            "  <p>Reviewer-facing figure package for continuous-trait model support across information criteria, likelihood, parameters, and fit diagnostics.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Information Criteria</h2><div class="figure-shell">' + figures["criteria"] + "</div></section>",
            '    <section class="panel"><h2>Likelihood</h2><div class="figure-shell">' + figures["likelihood"] + "</div></section>",
            '    <section class="panel"><h2>Parameters</h2><div class="figure-shell">' + figures["parameters"] + "</div></section>",
            '    <section class="panel"><h2>Fit Summary</h2><div class="figure-shell">' + figures["fit"] + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(criteria_table_path.name)}">{escape(criteria_table_path.name)}</a></li>',
            f'      <li><a href="{escape(likelihood_table_path.name)}">{escape(likelihood_table_path.name)}</a></li>',
            f'      <li><a href="{escape(parameter_table_path.name)}">{escape(parameter_table_path.name)}</a></li>',
            f'      <li><a href="{escape(fit_table_path.name)}">{escape(fit_table_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def build_comparative_model_figure_package(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    out_dir: Path,
    taxon_column: str | None = None,
) -> ComparativeModelFigurePackageResult:
    """Build a reviewer-facing figure package for BM versus OU continuous-trait model comparison."""
    out_dir.mkdir(parents=True, exist_ok=True)
    criteria_figure_path = out_dir / "model-comparison-criteria.svg"
    likelihood_figure_path = out_dir / "model-comparison-likelihood.svg"
    parameter_figure_path = out_dir / "model-comparison-parameters.svg"
    fit_figure_path = out_dir / "model-comparison-fit-summary.svg"
    criteria_table_path = out_dir / "model-comparison-criteria.tsv"
    likelihood_table_path = out_dir / "model-comparison-likelihood.tsv"
    parameter_table_path = out_dir / "model-comparison-parameters.tsv"
    fit_table_path = out_dir / "model-comparison-fit-summary.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "model-comparison-review.html"
    manifest_path = out_dir / "model-comparison-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    comparison_report = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    brownian_report = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    ou_report = fit_ornstein_uhlenbeck_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    criteria_rows = _criteria_rows(comparison_report)
    likelihood_rows = _likelihood_rows(comparison_report)
    parameter_rows = _parameter_rows(brownian_report, ou_report)
    fit_rows = _fit_rows(
        brownian_report,
        ou_report,
        selected_model=comparison_report.better_model,
    )
    _write_table(criteria_table_path, criteria_rows)
    _write_table(likelihood_table_path, likelihood_rows)
    _write_table(parameter_table_path, parameter_rows)
    _write_table(fit_table_path, fit_rows)
    _write_information_criteria_svg(
        criteria_figure_path,
        criteria_rows,
        selected_model=comparison_report.better_model,
    )
    _write_likelihood_svg(
        likelihood_figure_path,
        likelihood_rows,
        selected_model=comparison_report.better_model,
    )
    _write_parameter_svg(
        parameter_figure_path,
        parameter_rows,
        selected_model=comparison_report.better_model,
    )
    _write_fit_summary_svg(
        fit_figure_path,
        fit_rows,
        selected_model=comparison_report.better_model,
    )
    legend_entries = _legend_entries()
    _write_legend_table(legend_path, legend_entries)
    audit = _build_audit(
        comparison_report=comparison_report,
        brownian=brownian_report,
        ou=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
    )
    caption_draft = _build_caption_draft(
        comparison_report=comparison_report,
        audit=audit,
        fit_rows=fit_rows,
    )
    _write_caption(caption_path, caption_draft)
    review_path.write_text(
        _build_review_html(
            criteria_figure_path=criteria_figure_path,
            likelihood_figure_path=likelihood_figure_path,
            parameter_figure_path=parameter_figure_path,
            fit_figure_path=fit_figure_path,
            criteria_table_path=criteria_table_path,
            likelihood_table_path=likelihood_table_path,
            parameter_table_path=parameter_table_path,
            fit_table_path=fit_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            audit=audit,
        ),
        encoding="utf-8",
    )
    artifact_paths = [
        criteria_figure_path,
        likelihood_figure_path,
        parameter_figure_path,
        fit_figure_path,
        criteria_table_path,
        likelihood_table_path,
        parameter_table_path,
        fit_table_path,
        legend_path,
        caption_path,
        review_path,
    ]
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="comparative_model_figure_package",
        input_files=[
            ("tree", tree_path),
            ("traits", traits_path),
        ],
        generated_figures=[
            ("information_criteria", criteria_figure_path),
            ("likelihood", likelihood_figure_path),
            ("parameter_intervals", parameter_figure_path),
            ("fit_diagnostics", fit_figure_path),
        ],
        generated_tables=[
            ("information_criteria", criteria_table_path),
            ("likelihood", likelihood_table_path),
            ("parameter_intervals", parameter_table_path),
            ("fit_diagnostics", fit_table_path),
        ],
        filters=None,
        model={
            "kind": "comparative_model_comparison",
            "name": comparison_report.better_model,
            "candidate_models": [row.model for row in comparison_report.rows],
            "aicc_delta": audit.aicc_delta,
        },
        settings={
            "trait": trait,
            "taxon_column": taxon_column,
            "taxon_count": comparison_report.taxon_count,
        },
        linked_artifacts=[
            ("legend", legend_path),
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "comparative_model_figure_package",
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(
            reproducibility_manifest_path
        ),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "trait": trait,
            "taxon_column": taxon_column,
        },
        "metrics": {
            "taxon_count": comparison_report.taxon_count,
            "selected_model": audit.selected_model,
            "publication_ready": audit.publication_ready,
            "support_distinct": audit.support_distinct,
            "aicc_delta": audit.aicc_delta,
            "finite_aicc_model_count": audit.finite_aicc_model_count,
            "plotted_model_count": audit.plotted_model_count,
            "rendered_parameter_count": audit.rendered_parameter_count,
            "rendered_fit_row_count": audit.rendered_fit_row_count,
            "warning_count": audit.warning_count,
        },
        "comparison_report": _json_ready(asdict(comparison_report)),
        "brownian_report": _json_ready(asdict(brownian_report)),
        "ou_report": _json_ready(asdict(ou_report)),
        "criteria_rows": _json_ready([asdict(row) for row in criteria_rows]),
        "likelihood_rows": _json_ready([asdict(row) for row in likelihood_rows]),
        "parameter_rows": _json_ready([asdict(row) for row in parameter_rows]),
        "fit_rows": _json_ready([asdict(row) for row in fit_rows]),
        "audit": _json_ready(asdict(audit)),
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return ComparativeModelFigurePackageResult(
        output_dir=out_dir,
        criteria_figure_path=criteria_figure_path,
        likelihood_figure_path=likelihood_figure_path,
        parameter_figure_path=parameter_figure_path,
        fit_figure_path=fit_figure_path,
        criteria_table_path=criteria_table_path,
        likelihood_table_path=likelihood_table_path,
        parameter_table_path=parameter_table_path,
        fit_table_path=fit_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        comparison_report=comparison_report,
        brownian_report=brownian_report,
        ou_report=ou_report,
        criteria_rows=criteria_rows,
        likelihood_rows=likelihood_rows,
        parameter_rows=parameter_rows,
        fit_rows=fit_rows,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        machine_manifest=machine_manifest,
    )
