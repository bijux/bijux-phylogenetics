from __future__ import annotations

from html import escape
from pathlib import Path

from .contracts import ComparativeModelFigureAudit
from .summaries import format_model_figure_number


def build_model_figure_review_html(
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
        "<tr><th>" + escape(label) + "</th><td>" + escape(value) + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("selected_model", audit.selected_model),
            ("support_distinct", str(audit.support_distinct).lower()),
            ("aicc_delta", format_model_figure_number(audit.aicc_delta)),
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
            '    <section class="panel"><h2>Information Criteria</h2><div class="figure-shell">'
            + figures["criteria"]
            + "</div></section>",
            '    <section class="panel"><h2>Likelihood</h2><div class="figure-shell">'
            + figures["likelihood"]
            + "</div></section>",
            '    <section class="panel"><h2>Parameters</h2><div class="figure-shell">'
            + figures["parameters"]
            + "</div></section>",
            '    <section class="panel"><h2>Fit Summary</h2><div class="figure-shell">'
            + figures["fit"]
            + "</div></section>",
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
