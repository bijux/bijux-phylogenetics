from __future__ import annotations

from dataclasses import asdict
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist
from bijux_phylogenetics.trees import CladeTableRow

from .contracts import TreeBranchStatisticsRow, TreeSupportRow


def json_script(payload: dict[str, object]) -> str:
    """Render the embedded manifest JSON script block."""
    serialized = json.dumps(payload, default=str, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    return (
        '<script id="bijux-tree-report-package-manifest" type="application/json">'
        f"{serialized}</script>"
    )


def render_support_table(rows: list[TreeSupportRow]) -> str:
    """Render the reviewer-facing support summary table."""
    table_rows = []
    for row in rows:
        table_rows.append(
            "<tr>"
            f"<td>{escape(row.node_kind)}</td>"
            f"<td>{escape(row.node)}</td>"
            f"<td>{'' if row.support is None else row.support}</td>"
            f"<td>{'' if row.support_fraction is None else row.support_fraction}</td>"
            f"<td>{escape(row.support_class)}</td>"
            f"<td>{len(row.descendant_taxa)}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr>"
        "<th>node kind</th><th>node</th><th>support</th><th>support fraction</th>"
        "<th>support class</th><th>descendant taxa</th>"
        "</tr></thead><tbody>" + "".join(table_rows) + "</tbody></table>"
    )


def render_clade_table(rows: list[CladeTableRow]) -> str:
    """Render the reviewer-facing clade table."""
    internal_rows = [row for row in rows if row.node_kind != "tip"]
    table_rows = []
    for row in internal_rows:
        table_rows.append(
            "<tr>"
            f"<td>{escape(row.node_kind)}</td>"
            f"<td>{escape(row.clade_id)}</td>"
            f"<td>{row.taxon_count}</td>"
            f"<td>{'' if row.support is None else row.support}</td>"
            f"<td>{'' if row.branch_length is None else row.branch_length}</td>"
            f"<td>{'' if row.root_depth is None else row.root_depth}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr>"
        "<th>node kind</th><th>clade</th><th>taxa</th><th>support</th>"
        "<th>branch length</th><th>root depth</th>"
        "</tr></thead><tbody>" + "".join(table_rows) + "</tbody></table>"
    )


def render_branch_stats(row: TreeBranchStatisticsRow) -> str:
    """Render the reviewer-facing branch statistics table."""
    return (
        "<table><thead><tr>"
        "<th>branch count</th><th>defined</th><th>missing</th><th>zero</th>"
        "<th>negative</th><th>positive median</th><th>long outliers</th><th>short outliers</th>"
        "</tr></thead><tbody><tr>"
        f"<td>{row.branch_count}</td>"
        f"<td>{row.defined_branch_count}</td>"
        f"<td>{row.missing_branch_count}</td>"
        f"<td>{row.zero_length_branch_count}</td>"
        f"<td>{row.negative_branch_count}</td>"
        f"<td>{'' if row.positive_branch_median is None else row.positive_branch_median}</td>"
        f"<td>{row.long_outlier_count}</td>"
        f"<td>{row.short_outlier_count}</td>"
        "</tr></tbody></table>"
    )


def render_list(items: list[str]) -> str:
    """Render a reviewer-facing bullet list or explicit none block."""
    if not items:
        return "<p>none</p>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def render_method_tier(method_tier: MethodTierAssessment) -> str:
    """Render the method-tier narrative card."""
    basis = (
        "<p><strong>validation basis:</strong> "
        + escape("; ".join(method_tier.validation_basis))
        + "</p>"
        if method_tier.validation_basis
        else ""
    )
    warning = (
        f'<p class="warn"><strong>warning:</strong> {escape(method_tier.warning)}</p>'
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


def render_reviewer_audit_table(checklist: ReviewerAuditChecklist) -> str:
    """Render the reviewer-audit checklist table."""
    rows = []
    for item in checklist.items:
        rows.append(
            "<tr>"
            f"<td>{escape(item.section)}</td>"
            f"<td>{escape(item.status)}</td>"
            f"<td>{escape(item.summary)}</td>"
            f"<td>{escape('; '.join(item.evidence))}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr>"
        "<th>section</th><th>status</th><th>summary</th><th>evidence</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
    )


def write_tree_report_html(
    *,
    path: Path,
    title: str,
    figure_svg: str,
    methods_summary_text: str,
    reviewer_summary: list[str],
    limitations: list[str],
    support_rows: list[TreeSupportRow],
    clade_rows: list[CladeTableRow],
    branch_stats: TreeBranchStatisticsRow,
    method_tier: MethodTierAssessment,
    reviewer_audit_checklist: ReviewerAuditChecklist,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    manifest: dict[str, object],
) -> Path:
    """Write the reviewer-facing HTML package report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    :root {{
      --ink: #172033;
      --muted: #52607a;
      --bg: #f5f8fb;
      --panel: #ffffff;
      --rule: #d7e0ea;
      --accent: #0f766e;
      --accent-soft: #d9f2ee;
      --warn: #92400e;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top left, #e6f7f3 0, transparent 28rem),
        linear-gradient(180deg, #eef4f8 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 2rem;
    }}
    .shell {{
      background: rgba(255, 255, 255, 0.92);
      border: 1px solid var(--rule);
      border-radius: 24px;
      box-shadow: 0 28px 80px rgba(15, 23, 42, 0.08);
      padding: 2rem;
      backdrop-filter: blur(12px);
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin: 0 0 0.5rem;
      color: var(--accent);
    }}
    p.lead {{
      margin: 0 0 1.5rem;
      color: var(--muted);
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 0.75rem;
      margin: 1.5rem 0 2rem;
    }}
    .card {{
      background: linear-gradient(180deg, #ffffff 0%, #f8fbfd 100%);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1rem;
    }}
    .card .label {{
      color: var(--muted);
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .card .value {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 1.8rem;
      margin-top: 0.2rem;
    }}
    section + section {{
      margin-top: 2rem;
      padding-top: 2rem;
      border-top: 1px solid var(--rule);
    }}
    .figure-frame {{
      overflow-x: auto;
      background: #fbfdfe;
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1rem;
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
      font-family: var(--mono);
      font-size: 0.86rem;
    }}
    .note {{
      background: var(--accent-soft);
      border-left: 4px solid var(--accent);
      border-radius: 12px;
      padding: 0.9rem 1rem;
    }}
    .warn {{
      color: var(--warn);
    }}
  </style>
</head>
<body>
  <main>
    <div class="shell">
      <h1>{escape(title)}</h1>
      <p class="lead">Reviewer-facing tree evidence with figure, support, clade, and branch-length review in one package.</p>
      {json_script(manifest)}
      <div class="cards">
        <div class="card"><div class="label">Tip Count</div><div class="value">{inspection.tip_count}</div></div>
        <div class="card"><div class="label">Tree Quality</div><div class="value">{inspection.tree_quality_score}</div></div>
        <div class="card"><div class="label">Supported Branches</div><div class="value">{sum(1 for row in support_rows if row.support is not None)}</div></div>
        <div class="card"><div class="label">Method Tier</div><div class="value">{escape(method_tier.tier)}</div></div>
        <div class="card"><div class="label">Long Outliers</div><div class="value">{branch_stats.long_outlier_count}</div></div>
      </div>
      <section>
        <h2>Method Tier</h2>
        {render_method_tier(method_tier)}
      </section>
      <section>
        <h2>Reviewer Summary</h2>
        {render_list(reviewer_summary)}
        <div class="note">
          Support classes are reviewer-facing thresholds:
          <strong>strong</strong> at support fraction at least 0.95,
          <strong>moderate</strong> at least 0.80,
          and <strong>weak</strong> below 0.80.
        </div>
      </section>
      <section>
        <h2>Methods Summary</h2>
        <pre>{escape(methods_summary_text)}</pre>
      </section>
      <section>
        <h2>Reviewer Audit Checklist</h2>
        {render_reviewer_audit_table(reviewer_audit_checklist)}
      </section>
      <section>
        <h2>Tree Image</h2>
        <div class="figure-frame">{figure_svg}</div>
      </section>
      <section>
        <h2>Support Table</h2>
        {render_support_table(support_rows)}
      </section>
      <section>
        <h2>Clade Table</h2>
        {render_clade_table(clade_rows)}
      </section>
      <section>
        <h2>Branch-Length Stats</h2>
        {render_branch_stats(branch_stats)}
      </section>
      <section>
        <h2>Limitations</h2>
        {render_list(limitations)}
      </section>
      <section>
        <h2>Diagnostics</h2>
        <details open>
          <summary>Tree validation</summary>
          <pre>{escape(json.dumps(asdict(validation), default=str, indent=2, sort_keys=True))}</pre>
        </details>
        <details>
          <summary>Tree inspection</summary>
          <pre>{escape(json.dumps(asdict(inspection), default=str, indent=2, sort_keys=True))}</pre>
        </details>
        <details>
          <summary>Tree forensic</summary>
          <pre>{escape(json.dumps(asdict(forensic), default=str, indent=2, sort_keys=True))}</pre>
        </details>
      </section>
    </div>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    return path
