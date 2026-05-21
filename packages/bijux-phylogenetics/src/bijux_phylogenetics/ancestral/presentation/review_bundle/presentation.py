from __future__ import annotations

from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralSummary
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralSummary
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist


def write_package_html(
    *,
    path: Path,
    title: str,
    figure_svg: str,
    methods_summary_text: str,
    reconstruction_kind: str,
    model: str,
    summary: ContinuousAncestralSummary | DiscreteAncestralSummary,
    warnings: list[str],
    node_table_rows: list[list[str]],
    uncertainty_table_rows: list[list[str]],
    transition_count_rows: list[list[str]],
    transition_branch_rows: list[list[str]],
    limitations: list[str],
    reviewer_audit_checklist: ReviewerAuditChecklist,
    manifest: dict[str, object],
) -> Path:
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
      --bg: #f7f9fc;
      --panel: #ffffff;
      --rule: #d7e0ea;
      --accent: #1d4ed8;
      --accent-soft: #dbeafe;
    }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top right, #e4eefc 0, transparent 28rem),
        linear-gradient(180deg, #eef4fb 0%, var(--bg) 100%);
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
      font-size: 1.6rem;
      margin-top: 0.2rem;
    }}
    section + section {{
      margin-top: 2rem;
      padding-top: 2rem;
      border-top: 1px solid var(--rule);
    }}
    .figure-frame {{
      overflow-x: auto;
      background: #fbfdff;
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
      <h1>{escape(title)}</h1>
      <p class="lead">Reviewer-facing ancestral reconstruction bundle with node estimates, uncertainty, branchwise change evidence, and visualization in one package.</p>
      {_json_script(manifest)}
      <div class="cards">
        <div class="card"><div class="label">Kind</div><div class="value">{
        escape(reconstruction_kind)
    }</div></div>
        <div class="card"><div class="label">Model</div><div class="value">{
        escape(model)
    }</div></div>
        <div class="card"><div class="label">Analyzed Taxa</div><div class="value">{
        summary.analyzed_taxon_count
    }</div></div>
        <div class="card"><div class="label">Warnings</div><div class="value">{
        summary.warning_count
    }</div></div>
      </div>
      <section>
        <h2>Reviewer Summary</h2>
        {_list(warnings[:6])}
        <div class="note">
          Continuous packages preserve branch-change direction counts. Discrete packages preserve inferred state-transition counts.
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
        <h2>Tree Visualization</h2>
        <div class="figure-frame">{figure_svg}</div>
      </section>
      <section>
        <h2>Node Table</h2>
        {
        _table(
            [
                "node",
                "node name",
                "tip",
                "descendant taxa",
                "value or state",
                "uncertainty",
            ],
            node_table_rows,
        )
    }
      </section>
      <section>
        <h2>Uncertainty Review</h2>
        {
        _table(
            ["node", "descendant taxa", "uncertainty evidence", "interpretation"],
            uncertainty_table_rows,
        )
    }
      </section>
      <section>
        <h2>Transition Review</h2>
        {
        _table(
            ["label", "count", "fraction or certainty", "detail"],
            transition_count_rows,
        )
    }
        {
        _table(
            ["parent", "child", "descendant taxa", "branch length", "change", "detail"],
            transition_branch_rows,
        )
    }
      </section>
      <section>
        <h2>Limitations</h2>
        {_list(limitations)}
      </section>
    </div>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    return path


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, default=str, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    return (
        '<script id="bijux-ancestral-report-package-manifest" type="application/json">'
        f"{serialized}</script>"
    )


def _table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _list(items: list[str]) -> str:
    if not items:
        return "<p>none</p>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _reviewer_audit_table(checklist: ReviewerAuditChecklist) -> str:
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
