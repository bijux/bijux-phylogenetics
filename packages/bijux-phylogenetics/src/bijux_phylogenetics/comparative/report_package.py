from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)


@dataclass(frozen=True, slots=True)
class ComparativeAnalysisSummaryRow:
    response: str
    formula: str
    predictor_count: int
    analysis_taxa: int
    excluded_taxa: int
    selected_model: str
    pgls_lambda: float
    pgls_log_likelihood: float
    pgls_r_squared: float
    phylogenetic_signal_k: float
    phylogenetic_signal_lambda: float
    independent_contrast_count: int
    better_model_aicc_delta: float


@dataclass(frozen=True, slots=True)
class ComparativeCoefficientTableRow:
    term: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    degrees_of_freedom: int
    inference_distribution: str
    significant: bool


@dataclass(frozen=True, slots=True)
class ComparativeResidualTableRow:
    analysis: str
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float | None
    max_leverage: float | None
    outlier_taxa: tuple[str, ...]
    high_leverage_taxa: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ComparativeSignalTableRow:
    trait: str
    taxon_count: int
    blombergs_k: float
    pagels_lambda: float
    lambda_log_likelihood: float
    lambda_null_log_likelihood: float
    lambda_brownian_log_likelihood: float
    independent_contrast_count: int
    independent_contrast_root_estimate: float


@dataclass(frozen=True, slots=True)
class ComparativeInterpretationRow:
    topic: str
    claim: str
    evidence: str
    caution: str


@dataclass(frozen=True, slots=True)
class ComparativeAuditTableRow:
    analysis: str
    taxa_used: tuple[str, ...]
    traits_used: tuple[str, ...]
    excluded_taxa: tuple[str, ...]
    assumptions: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(slots=True)
class ComparativeReportPackageResult:
    output_dir: Path
    report_path: Path
    summary_table_path: Path
    coefficient_table_path: Path
    residual_table_path: Path
    signal_table_path: Path
    model_comparison_table_path: Path
    interpretation_table_path: Path
    audit_table_path: Path
    contrast_table_path: Path
    manifest_path: Path
    report: ComparativeMethodReport
    summary_row: ComparativeAnalysisSummaryRow
    coefficient_rows: list[ComparativeCoefficientTableRow]
    residual_rows: list[ComparativeResidualTableRow]
    signal_row: ComparativeSignalTableRow
    interpretation_rows: list[ComparativeInterpretationRow]
    audit_rows: list[ComparativeAuditTableRow]
    machine_manifest: dict[str, object]


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def summarize_comparative_analysis(
    report: ComparativeMethodReport,
) -> ComparativeAnalysisSummaryRow:
    model_rows = report.snapshot.model_comparison.rows
    selected_row = next(row for row in model_rows if row.selected)
    runner_up_aicc = min(
        row.aicc for row in model_rows if row.model != selected_row.model
    )
    return ComparativeAnalysisSummaryRow(
        response=report.snapshot.response,
        formula=report.snapshot.formula.formula,
        predictor_count=len(report.snapshot.formula.predictors),
        analysis_taxa=report.snapshot.pgls_model.taxon_count,
        excluded_taxa=len(report.snapshot.pgls_inputs.formula_audit.excluded_taxa),
        selected_model=report.snapshot.model_comparison.better_model,
        pgls_lambda=report.snapshot.pgls_model.lambda_value,
        pgls_log_likelihood=report.snapshot.pgls_model.log_likelihood,
        pgls_r_squared=report.snapshot.pgls_model.r_squared,
        phylogenetic_signal_k=report.snapshot.signal_k.k,
        phylogenetic_signal_lambda=report.snapshot.signal_lambda.lambda_value,
        independent_contrast_count=len(report.snapshot.contrasts.contrasts),
        better_model_aicc_delta=runner_up_aicc - selected_row.aicc,
    )


def summarize_comparative_coefficients(
    report: ComparativeMethodReport,
) -> list[ComparativeCoefficientTableRow]:
    rows: list[ComparativeCoefficientTableRow] = []
    for coefficient in report.snapshot.pgls_model.coefficients:
        rows.append(
            ComparativeCoefficientTableRow(
                term=coefficient.name,
                estimate=coefficient.estimate,
                standard_error=coefficient.standard_error,
                test_statistic=coefficient.test_statistic,
                p_value=coefficient.p_value,
                lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
                upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
                degrees_of_freedom=coefficient.degrees_of_freedom,
                inference_distribution=coefficient.inference_distribution,
                significant=coefficient.p_value <= 0.05,
            )
        )
    return rows


def summarize_comparative_residuals(
    report: ComparativeMethodReport,
) -> list[ComparativeResidualTableRow]:
    rows: list[ComparativeResidualTableRow] = []
    for surface in report.snapshot.maturity.residual_diagnostics:
        rows.append(
            ComparativeResidualTableRow(
                analysis=surface.analysis,
                residual_variance=surface.residual_variance,
                max_abs_standardized_residual=surface.max_abs_standardized_residual,
                phylogenetic_residual_lambda=surface.phylogenetic_residual_lambda,
                max_leverage=surface.max_leverage,
                outlier_taxa=tuple(surface.outlier_taxa),
                high_leverage_taxa=tuple(surface.high_leverage_taxa),
                warnings=tuple(surface.warnings),
            )
        )
    return rows


def summarize_comparative_signal(
    report: ComparativeMethodReport,
) -> ComparativeSignalTableRow:
    return ComparativeSignalTableRow(
        trait=report.snapshot.signal_k.trait,
        taxon_count=report.snapshot.signal_k.taxon_count,
        blombergs_k=report.snapshot.signal_k.k,
        pagels_lambda=report.snapshot.signal_lambda.lambda_value,
        lambda_log_likelihood=report.snapshot.signal_lambda.log_likelihood,
        lambda_null_log_likelihood=report.snapshot.signal_lambda.null_log_likelihood,
        lambda_brownian_log_likelihood=(
            report.snapshot.signal_lambda.brownian_log_likelihood
        ),
        independent_contrast_count=len(report.snapshot.contrasts.contrasts),
        independent_contrast_root_estimate=report.snapshot.contrasts.root_estimate,
    )


def summarize_comparative_interpretation(
    report: ComparativeMethodReport,
) -> list[ComparativeInterpretationRow]:
    rows: list[ComparativeInterpretationRow] = []
    summary = summarize_comparative_analysis(report)
    rows.append(
        ComparativeInterpretationRow(
            topic="formula",
            claim=f"comparative regression was fit as `{summary.formula}`",
            evidence=(
                f"analysis_taxa={summary.analysis_taxa}; predictors={summary.predictor_count}"
            ),
            caution=(
                "formula interpretation is conditional on the supplied tree, taxon overlap, and encoded predictor structure"
            ),
        )
    )
    signal_level = (
        "high"
        if report.snapshot.signal_lambda.lambda_value >= 0.8
        else "moderate"
        if report.snapshot.signal_lambda.lambda_value >= 0.3
        else "low"
    )
    rows.append(
        ComparativeInterpretationRow(
            topic="phylogenetic-signal",
            claim=f"the response trait retains {signal_level} phylogenetic structure",
            evidence=(
                f"blombergs_k={report.snapshot.signal_k.k:.6f}; "
                f"pagels_lambda={report.snapshot.signal_lambda.lambda_value:.6f}"
            ),
            caution=(
                "signal metrics describe covariance with the supplied phylogeny and are not by themselves a process explanation"
            ),
        )
    )
    rows.append(
        ComparativeInterpretationRow(
            topic="model-comparison",
            claim=(
                f"{report.snapshot.model_comparison.better_model} is the preferred continuous trait model on AICc"
            ),
            evidence=(
                f"aicc_delta={summary.better_model_aicc_delta:.6f}; "
                f"selected_from={len(report.snapshot.model_comparison.rows)} models"
            ),
            caution=(
                "small AICc differences should not be overstated as decisive biological process separation"
            ),
        )
    )
    for row in summarize_comparative_coefficients(report):
        if row.term == "intercept":
            continue
        direction = "positive" if row.estimate > 0.0 else "negative"
        strength = (
            "nominally supported" if row.significant else "not nominally supported"
        )
        rows.append(
            ComparativeInterpretationRow(
                topic="coefficient",
                claim=f"`{row.term}` shows a {direction} association and is {strength}",
                evidence=(
                    f"estimate={row.estimate:.6f}; p_value={row.p_value:.6f}; "
                    f"95%_ci=[{row.lower_95_confidence_interval:.6f}, {row.upper_95_confidence_interval:.6f}]"
                ),
                caution=(
                    "coefficient interpretation remains associational and depends on lambda choice, sample size, and residual diagnostics"
                ),
            )
        )
    for residual in summarize_comparative_residuals(report):
        if residual.warnings:
            rows.append(
                ComparativeInterpretationRow(
                    topic="diagnostics",
                    claim=f"{residual.analysis} residual diagnostics flagged review concerns",
                    evidence="; ".join(residual.warnings),
                    caution=(
                        "diagnostic warnings weaken confidence in coefficient interpretation and should be reviewed before publication use"
                    ),
                )
            )
    return rows


def summarize_comparative_audit(
    report: ComparativeMethodReport,
) -> list[ComparativeAuditTableRow]:
    return [
        ComparativeAuditTableRow(
            analysis=row.analysis,
            taxa_used=tuple(row.taxa_used),
            traits_used=tuple(row.traits_used),
            excluded_taxa=tuple(row.excluded_taxa),
            assumptions=tuple(row.assumptions),
            warnings=tuple(row.warnings),
        )
        for row in report.snapshot.audit_rows
    ]


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


def _write_comparative_report_html(
    *,
    path: Path,
    report: ComparativeMethodReport,
    summary_row: ComparativeAnalysisSummaryRow,
    coefficient_rows: list[ComparativeCoefficientTableRow],
    residual_rows: list[ComparativeResidualTableRow],
    signal_row: ComparativeSignalTableRow,
    interpretation_rows: list[ComparativeInterpretationRow],
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
        <div class="card"><div class="label">PGLS R²</div><div class="value">{
        summary_row.pgls_r_squared:.3f}</div></div>
      </div>
      <section>
        <h2>Reviewer Summary</h2>
        {_list(report.snapshot.limitations[:4])}
        <div class="note">
          Comparative coefficients remain associational. Review residual diagnostics, phylogenetic signal, and model-comparison tables before treating direction or significance as robust.
        </div>
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

    report_path = out_dir / "comparative-report.html"
    summary_table_path = out_dir / "comparative-summary.tsv"
    coefficient_table_path = out_dir / "coefficient-table.tsv"
    residual_table_path = out_dir / "residual-summary.tsv"
    signal_table_path = out_dir / "signal-summary.tsv"
    model_comparison_table_path = out_dir / "model-comparison.tsv"
    interpretation_table_path = out_dir / "interpretation-table.tsv"
    audit_table_path = out_dir / "audit-table.tsv"
    contrast_table_path = out_dir / "contrast-table.tsv"
    manifest_path = out_dir / "comparative-report.manifest.json"

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
        },
        "summary": asdict(summary_row),
        "limitations": report.snapshot.limitations,
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_comparative_report_html(
        path=report_path,
        report=report,
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        manifest=machine_manifest,
    )
    return ComparativeReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
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
        summary_row=summary_row,
        coefficient_rows=coefficient_rows,
        residual_rows=residual_rows,
        signal_row=signal_row,
        interpretation_rows=interpretation_rows,
        audit_rows=audit_rows,
        machine_manifest=machine_manifest,
    )
