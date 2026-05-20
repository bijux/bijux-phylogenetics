from __future__ import annotations

import csv
from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    ComparativeMethodsSummaryTextResult,
    build_comparative_method_report,
    write_comparative_methods_summary_text,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    comparative_report_method_tier,
)
from bijux_phylogenetics.reports.review import (
    ReviewerAuditChecklist,
    write_reviewer_audit_checklist,
)
from .contracts import (
    ComparativeAnalysisSummaryRow,
    ComparativeAuditTableRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeReportPackageResult,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


from .summaries import (
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
)


def _write_rows(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=_table_delimiter(path)
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_comparative_summary_table(
    path: Path, row: ComparativeAnalysisSummaryRow
) -> Path:
    return _write_rows(path, list(asdict(row).keys()), [asdict(row)])


def write_comparative_coefficient_table(
    path: Path, rows: list[ComparativeCoefficientTableRow]
) -> Path:
    return _write_rows(
        path, list(asdict(rows[0]).keys()), [asdict(row) for row in rows]
    )


def write_comparative_residual_table(
    path: Path, rows: list[ComparativeResidualTableRow]
) -> Path:
    rendered = []
    for row in rows:
        payload = asdict(row)
        payload["outlier_taxa"] = "|".join(row.outlier_taxa)
        payload["high_leverage_taxa"] = "|".join(row.high_leverage_taxa)
        payload["warnings"] = "|".join(row.warnings)
        rendered.append(payload)
    return _write_rows(path, list(rendered[0].keys()), rendered)


def write_comparative_signal_table(path: Path, row: ComparativeSignalTableRow) -> Path:
    return _write_rows(path, list(asdict(row).keys()), [asdict(row)])


def write_comparative_model_comparison_table(
    path: Path, report: ComparativeMethodReport
) -> Path:
    rows = [asdict(row) for row in report.snapshot.model_comparison.rows]
    return _write_rows(path, list(rows[0].keys()), rows)


def write_comparative_interpretation_table(
    path: Path, rows: list[ComparativeInterpretationRow]
) -> Path:
    return _write_rows(
        path, list(asdict(rows[0]).keys()), [asdict(row) for row in rows]
    )


def write_comparative_audit_table(
    path: Path, rows: list[ComparativeAuditTableRow]
) -> Path:
    rendered = []
    for row in rows:
        payload = asdict(row)
        payload["taxa_used"] = "|".join(row.taxa_used)
        payload["traits_used"] = "|".join(row.traits_used)
        payload["excluded_taxa"] = "|".join(row.excluded_taxa)
        payload["assumptions"] = "|".join(row.assumptions)
        payload["warnings"] = "|".join(row.warnings)
        rendered.append(payload)
    return _write_rows(path, list(rendered[0].keys()), rendered)


def write_comparative_contrast_table(
    path: Path, report: ComparativeMethodReport
) -> Path:
    rows = [
        {
            "node": row.node,
            "left_taxa": "|".join(row.left_taxa),
            "right_taxa": "|".join(row.right_taxa),
            "contrast": row.contrast,
            "expected_variance": row.expected_variance,
            "ancestral_value": row.ancestral_value,
        }
        for row in report.snapshot.contrasts.contrasts
    ]
    return _write_rows(path, list(rows[0].keys()), rows)


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


def _write_comparative_report_html(
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


def build_comparative_report_package(
    tree_path: Path,
    traits_path: Path,
    *,
    out_dir: Path,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeReportPackageResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    report = build_comparative_method_report(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    summary_row = summarize_comparative_analysis(report)
    coefficient_rows = summarize_comparative_coefficients(report)
    residual_rows = summarize_comparative_residuals(report)
    signal_row = summarize_comparative_signal(report)
    interpretation_rows = summarize_comparative_interpretation(report)
    audit_rows = summarize_comparative_audit(report)
    method_tier = comparative_report_method_tier()

    report_path = out_dir / "comparative-report.html"
    methods_summary_path = out_dir / "comparative-methods-summary.md"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"
    summary_table_path = out_dir / "comparative-summary.tsv"
    coefficient_table_path = out_dir / "coefficient-table.tsv"
    residual_table_path = out_dir / "residual-summary.tsv"
    signal_table_path = out_dir / "signal-summary.tsv"
    model_comparison_table_path = out_dir / "model-comparison.tsv"
    interpretation_table_path = out_dir / "interpretation-table.tsv"
    audit_table_path = out_dir / "audit-table.tsv"
    contrast_table_path = out_dir / "contrast-table.tsv"
    manifest_path = out_dir / "comparative-report.manifest.json"

    methods_summary = write_comparative_methods_summary_text(
        methods_summary_path, report
    )
    write_comparative_summary_table(summary_table_path, summary_row)
    write_comparative_coefficient_table(coefficient_table_path, coefficient_rows)
    write_comparative_residual_table(residual_table_path, residual_rows)
    write_comparative_signal_table(signal_table_path, signal_row)
    write_comparative_model_comparison_table(model_comparison_table_path, report)
    write_comparative_interpretation_table(
        interpretation_table_path, interpretation_rows
    )
    write_comparative_audit_table(audit_table_path, audit_rows)
    write_comparative_contrast_table(contrast_table_path, report)

    machine_manifest = {
        "report_kind": "comparative_package",
        "input_paths": [str(tree_path), str(traits_path)],
        "input_checksums": {
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
        },
        "outputs": {
            "report_path": str(report_path),
            "methods_summary_path": str(methods_summary_path),
            "reviewer_audit_checklist_path": str(reviewer_audit_checklist_path),
            "summary_table_path": str(summary_table_path),
            "coefficient_table_path": str(coefficient_table_path),
            "residual_table_path": str(residual_table_path),
            "signal_table_path": str(signal_table_path),
            "model_comparison_table_path": str(model_comparison_table_path),
            "interpretation_table_path": str(interpretation_table_path),
            "audit_table_path": str(audit_table_path),
            "contrast_table_path": str(contrast_table_path),
        },
        "metrics": {
            "analysis_taxa": summary_row.analysis_taxa,
            "selected_model": summary_row.selected_model,
            "coefficient_count": len(coefficient_rows),
            "contrast_count": signal_row.independent_contrast_count,
            "limitation_count": len(report.snapshot.limitations),
            "methods_summary_warning_count": methods_summary.warning_count,
        },
        "summary": asdict(summary_row),
        "limitations": report.snapshot.limitations,
    }
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["reviewer_audit_checklist"] = asdict(reviewer_audit_checklist)
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_comparative_report_html(
        path=report_path,
        report=report,
        methods_summary_text=methods_summary.text,
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        method_tier=method_tier,
        reviewer_audit_checklist=reviewer_audit_checklist,
        manifest=machine_manifest,
    )
    return ComparativeReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        methods_summary_path=methods_summary_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        summary_table_path=summary_table_path,
        coefficient_table_path=coefficient_table_path,
        residual_table_path=residual_table_path,
        signal_table_path=signal_table_path,
        model_comparison_table_path=model_comparison_table_path,
        interpretation_table_path=interpretation_table_path,
        audit_table_path=audit_table_path,
        contrast_table_path=contrast_table_path,
        manifest_path=manifest_path,
        report=report,
        methods_summary=methods_summary,
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        audit_rows=audit_rows,
        method_tier=method_tier,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )
