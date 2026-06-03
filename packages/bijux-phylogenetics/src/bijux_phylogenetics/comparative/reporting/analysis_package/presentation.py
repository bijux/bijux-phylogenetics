from __future__ import annotations

from dataclasses import asdict
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import ComparativeMethodReport
from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist

from .contracts import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, default=str, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    return (
        '<script id="bijux-comparative-report-manifest" type="application/json">'
        f"{serialized}</script>"
    )


def _list(items: list[str]) -> str:
    if not items:
        return "<p>none</p>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        )
    return (
        f"<table><thead><tr>{head}</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table>"
    )


def _method_tier_note(method_tier: MethodTierAssessment) -> str:
    basis = (
        "<p><strong>validation basis:</strong> "
        + escape("; ".join(method_tier.validation_basis))
        + "</p>"
        if method_tier.validation_basis
        else ""
    )
    warning = (
        f"<p><strong>warning:</strong> {escape(method_tier.warning)}</p>"
        if method_tier.warning is not None
        else ""
    )
    approximation = (
        "<p><strong>approximation:</strong> "
        + escape(method_tier.approximation)
        + "</p>"
        if method_tier.approximation is not None
        else ""
    )
    return (
        '<div class="note">'
        f"<p><strong>{escape(method_tier.tier)}</strong> ({escape(method_tier.inference_mode)})</p>"
        f"<p>{escape(method_tier.summary)}</p>"
        f"{basis}{approximation}{warning}"
        "</div>"
    )


def _reviewer_audit_table(checklist: ReviewerAuditChecklist) -> str:
    return _table(
        ["section", "status", "summary", "evidence"],
        [
            [
                item.section,
                item.status,
                item.summary,
                "; ".join(item.evidence),
            ]
            for item in checklist.items
        ],
    )


def write_comparative_report_html(
    *,
    path: Path,
    report: ComparativeMethodReport,
    methods_summary_text: str,
    summary_row: ComparativeAnalysisSummaryRow,
    coefficient_rows: list[ComparativeCoefficientTableRow],
    residual_rows: list[ComparativeResidualTableRow],
    signal_row: ComparativeSignalTableRow,
    interpretation_rows: list[ComparativeInterpretationRow],
    method_tier: MethodTierAssessment,
    reviewer_audit_checklist: ReviewerAuditChecklist,
    manifest: dict[str, object],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Bijux Comparative Analysis Report</title>
  <style>
    :root {{
      --ink: #142132;
      --muted: #5a667d;
      --bg: #f5f9fb;
      --panel: #ffffff;
      --rule: #d7e3ec;
      --accent: #1d4ed8;
      --accent-soft: #dbeafe;
    }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top right, #e3f0ff 0, transparent 30rem),
        linear-gradient(180deg, #eef5fa 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 2rem;
    }}
    .shell {{
      background: rgba(255, 255, 255, 0.94);
      border: 1px solid var(--rule);
      border-radius: 24px;
      box-shadow: 0 28px 80px rgba(15, 23, 42, 0.08);
      padding: 2rem;
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin: 0 0 0.4rem;
      color: var(--accent);
    }}
    .lead {{
      margin: 0 0 1.5rem;
      color: var(--muted);
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
      gap: 0.75rem;
      margin: 1.5rem 0 2rem;
    }}
    .card {{
      background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1rem;
    }}
    .label {{
      color: var(--muted);
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .value {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 1.7rem;
      margin-top: 0.2rem;
    }}
    section + section {{
      margin-top: 2rem;
      padding-top: 2rem;
      border-top: 1px solid var(--rule);
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.92rem;
    }}
    th, td {{
      text-align: left;
      padding: 0.65rem 0.55rem;
      border-bottom: 1px solid var(--rule);
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-size: 0.76rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    details {{
      margin-top: 1rem;
      border: 1px solid var(--rule);
      border-radius: 14px;
      background: #fbfdff;
      padding: 0.85rem 1rem;
    }}
    summary {{
      cursor: pointer;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-weight: 600;
    }}
    pre {{
      margin: 1rem 0 0;
      white-space: pre-wrap;
      word-break: break-word;
      background: #f6f8fa;
      border-radius: 12px;
      padding: 1rem;
      font-family: "SFMono-Regular", "SF Mono", Consolas, monospace;
      font-size: 0.86rem;
    }}
    .note {{
      background: var(--accent-soft);
      border-left: 4px solid var(--accent);
      border-radius: 12px;
      padding: 0.9rem 1rem;
    }}
  </style>
</head>
<body>
  <main>
    <div class="shell">
      <h1>Bijux Comparative Analysis Report</h1>
      <p class="lead">Integrated reviewer-facing output for comparative regression, phylogenetic signal, contrasts, residual diagnostics, and process-model comparison.</p>
      {_json_script(manifest)}
      <div class="cards">
        <div class="card"><div class="label">Formula</div><div class="value">{
        escape(summary_row.formula)
    }</div></div>
        <div class="card"><div class="label">Analysis Taxa</div><div class="value">{
        summary_row.analysis_taxa
    }</div></div>
        <div class="card"><div class="label">Selected Model</div><div class="value">{
        escape(summary_row.selected_model)
    }</div></div>
        <div class="card"><div class="label">Method Tier</div><div class="value">{
        escape(method_tier.tier)
    }</div></div>
        <div class="card"><div class="label">PGLS R²</div><div class="value">{
        summary_row.pgls_r_squared:.3f}</div></div>
      </div>
      <section>
        <h2>Method Tier</h2>
        {_method_tier_note(method_tier)}
      </section>
      <section>
        <h2>Reviewer Summary</h2>
        {_list(report.snapshot.limitations[:4])}
        <div class="note">
          Comparative coefficients remain associational. Review residual diagnostics, phylogenetic signal, and model-comparison tables before treating direction or significance as robust.
        </div>
      </section>
      <section>
        <h2>Methods Summary</h2>
        <pre>{escape(methods_summary_text)}</pre>
      </section>
      <section>
        <h2>Reviewer Audit Checklist</h2>
        {_reviewer_audit_table(reviewer_audit_checklist)}
      </section>
      <section>
        <h2>Coefficient Table</h2>
        {
        _table(
            [
                "term",
                "estimate",
                "standard error",
                "test statistic",
                "p value",
                "95% CI",
                "significant",
            ],
            [
                [
                    row.term,
                    f"{row.estimate:.6f}",
                    f"{row.standard_error:.6f}",
                    f"{row.test_statistic:.6f}",
                    f"{row.p_value:.6f}",
                    f"[{row.lower_95_confidence_interval:.6f}, {row.upper_95_confidence_interval:.6f}]",
                    "true" if row.significant else "false",
                ]
                for row in coefficient_rows
            ],
        )
    }
      </section>
      <section>
        <h2>Residual Summary</h2>
        {
        _table(
            [
                "analysis",
                "residual variance",
                "max |standardized residual|",
                "residual lambda",
                "max leverage",
                "outlier taxa",
            ],
            [
                [
                    row.analysis,
                    f"{row.residual_variance:.6f}",
                    f"{row.max_abs_standardized_residual:.6f}",
                    ""
                    if row.phylogenetic_residual_lambda is None
                    else f"{row.phylogenetic_residual_lambda:.6f}",
                    "" if row.max_leverage is None else f"{row.max_leverage:.6f}",
                    ", ".join(row.outlier_taxa),
                ]
                for row in residual_rows
            ],
        )
    }
      </section>
      <section>
        <h2>Phylogenetic Signal</h2>
        {
        _table(
            [
                "trait",
                "blomberg's k",
                "pagel's lambda",
                "contrast count",
                "root estimate",
            ],
            [
                [
                    signal_row.trait,
                    f"{signal_row.blombergs_k:.6f}",
                    f"{signal_row.pagels_lambda:.6f}",
                    str(signal_row.independent_contrast_count),
                    f"{signal_row.independent_contrast_root_estimate:.6f}",
                ]
            ],
        )
    }
      </section>
      <section>
        <h2>Model Comparison</h2>
        {
        _table(
            ["model", "parameter count", "log likelihood", "AIC", "AICc", "selected"],
            [
                [
                    row.model,
                    str(row.parameter_count),
                    f"{row.log_likelihood:.6f}",
                    f"{row.aic:.6f}",
                    f"{row.aicc:.6f}",
                    "true" if row.selected else "false",
                ]
                for row in report.snapshot.model_comparison.rows
            ],
        )
    }
      </section>
      <section>
        <h2>Biological Interpretation</h2>
        {
        _table(
            ["topic", "claim", "evidence", "caution"],
            [
                [row.topic, row.claim, row.evidence, row.caution]
                for row in interpretation_rows
            ],
        )
    }
      </section>
      <section>
        <h2>Diagnostics</h2>
        <details open>
          <summary>PGLS inputs</summary>
          <pre>{
        escape(
            json.dumps(
                asdict(report.snapshot.pgls_inputs),
                default=str,
                indent=2,
                sort_keys=True,
            )
        )
    }</pre>
        </details>
        <details>
          <summary>PGLS model</summary>
          <pre>{
        escape(
            json.dumps(
                asdict(report.snapshot.pgls_model),
                default=str,
                indent=2,
                sort_keys=True,
            )
        )
    }</pre>
        </details>
        <details>
          <summary>Influence</summary>
          <pre>{
        escape(
            json.dumps(asdict(report.influence), default=str, indent=2, sort_keys=True)
        )
    }</pre>
        </details>
      </section>
    </div>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    return path
