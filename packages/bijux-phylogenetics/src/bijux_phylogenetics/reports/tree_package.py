from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.branch_lengths import (
    BranchLengthDistributionReport,
    analyze_branch_length_distribution,
)
from bijux_phylogenetics.clades import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
)
from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.iqtree_support import support_fraction
from bijux_phylogenetics.render.svg import (
    SupportLabelRenderAudit,
    TreeRenderResult,
    audit_support_label_rendering,
    render_tree_svg,
)


@dataclass(frozen=True, slots=True)
class TreeSupportRow:
    node_kind: str
    node: str
    node_label: str | None
    descendant_taxa: tuple[str, ...]
    support: float | None
    support_fraction: float | None
    support_class: str
    branch_length: float | None
    root_depth: float | None


@dataclass(frozen=True, slots=True)
class TreeBranchStatisticsRow:
    branch_count: int
    defined_branch_count: int
    missing_branch_count: int
    zero_length_branch_count: int
    negative_branch_count: int
    positive_branch_count: int
    long_outlier_count: int
    short_outlier_count: int
    minimum_branch_length: float | None
    maximum_branch_length: float | None
    mean_branch_length: float | None
    median_branch_length: float | None
    positive_branch_median: float | None


@dataclass(slots=True)
class TreeReportPackageResult:
    output_dir: Path
    report_path: Path
    figure_path: Path
    support_table_path: Path
    clade_table_path: Path
    branch_stats_path: Path
    manifest_path: Path
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    forensic: TreeForensicReport
    figure: TreeRenderResult
    support_audit: SupportLabelRenderAudit
    clades: CladeTableReport
    branch_lengths: BranchLengthDistributionReport
    support_rows: list[TreeSupportRow]
    branch_stats: TreeBranchStatisticsRow
    reviewer_summary: list[str]
    limitations: list[str]
    machine_manifest: dict[str, object]


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _support_class(value: float | None) -> str:
    fraction = support_fraction(value)
    if fraction is None:
        return "missing"
    if fraction >= 0.95:
        return "strong"
    if fraction >= 0.80:
        return "moderate"
    return "weak"


def summarize_tree_support(clades: CladeTableReport) -> list[TreeSupportRow]:
    rows: list[TreeSupportRow] = []
    for row in clades.rows:
        if row.node_kind == "tip":
            continue
        rows.append(
            TreeSupportRow(
                node_kind=row.node_kind,
                node=row.clade_id,
                node_label=row.node_label,
                descendant_taxa=tuple(row.taxa),
                support=row.support,
                support_fraction=row.support_fraction,
                support_class=_support_class(row.support),
                branch_length=row.branch_length,
                root_depth=row.root_depth,
            )
        )
    return rows


def summarize_tree_branch_statistics(
    branch_lengths: BranchLengthDistributionReport,
) -> TreeBranchStatisticsRow:
    aggregate = branch_lengths.aggregate
    return TreeBranchStatisticsRow(
        branch_count=aggregate.branch_count,
        defined_branch_count=aggregate.defined_branch_count,
        missing_branch_count=aggregate.missing_branch_count,
        zero_length_branch_count=aggregate.zero_length_branch_count,
        negative_branch_count=aggregate.negative_branch_count,
        positive_branch_count=aggregate.positive_branch_count,
        long_outlier_count=aggregate.long_outlier_count,
        short_outlier_count=aggregate.short_outlier_count,
        minimum_branch_length=aggregate.minimum_branch_length,
        maximum_branch_length=aggregate.maximum_branch_length,
        mean_branch_length=aggregate.mean_branch_length,
        median_branch_length=aggregate.median_branch_length,
        positive_branch_median=aggregate.positive_branch_median,
    )


def write_tree_support_table(path: Path, rows: list[TreeSupportRow]) -> Path:
    fieldnames = [
        "node_kind",
        "node",
        "node_label",
        "descendant_taxa",
        "support",
        "support_fraction",
        "support_class",
        "branch_length",
        "root_depth",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=_table_delimiter(path)
        )
        writer.writeheader()
        for row in rows:
            payload = asdict(row)
            payload["node_label"] = "" if row.node_label is None else row.node_label
            payload["descendant_taxa"] = "|".join(row.descendant_taxa)
            payload["support"] = "" if row.support is None else row.support
            payload["support_fraction"] = (
                "" if row.support_fraction is None else row.support_fraction
            )
            payload["branch_length"] = (
                "" if row.branch_length is None else row.branch_length
            )
            payload["root_depth"] = "" if row.root_depth is None else row.root_depth
            writer.writerow(payload)
    return path


def write_tree_branch_statistics_table(
    path: Path, row: TreeBranchStatisticsRow
) -> Path:
    fieldnames = list(asdict(row).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=_table_delimiter(path)
        )
        writer.writeheader()
        writer.writerow(asdict(row))
    return path


def _reviewer_summary(
    *,
    inspection: TreeInspectionReport,
    support_rows: list[TreeSupportRow],
    branch_stats: TreeBranchStatisticsRow,
    support_audit: SupportLabelRenderAudit,
) -> tuple[list[str], list[str]]:
    supported_rows = [row for row in support_rows if row.support is not None]
    strong_rows = [row for row in support_rows if row.support_class == "strong"]
    limitations: list[str] = []
    if not support_audit.validated:
        limitations.append(
            "support labels were withheld from the rendered figure because the input support surface was not safe to standardize"
        )
    if branch_stats.missing_branch_count:
        limitations.append(
            "branch-length summaries include missing lengths, so weighted interpretation is incomplete"
        )
    if branch_stats.negative_branch_count:
        limitations.append(
            "negative branch lengths remain in the source tree and should be corrected before downstream weighted analysis"
        )
    summary = [
        f"tree quality score: {inspection.tree_quality_score}",
        f"tip count: {inspection.tip_count}",
        f"internal clade count: {sum(1 for row in support_rows if row.node_kind == 'internal')}",
        f"supported branch count: {len(supported_rows)}",
        f"strong-support branch count: {len(strong_rows)}",
        f"long-branch outlier count: {branch_stats.long_outlier_count}",
    ]
    if support_audit.validated and support_audit.warnings:
        summary.extend(support_audit.warnings)
    return summary, limitations


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, default=str, indent=2, sort_keys=True).replace(
        "</", "<\\/"
    )
    return (
        '<script id="bijux-tree-report-package-manifest" type="application/json">'
        f"{serialized}</script>"
    )


def _render_support_table(rows: list[TreeSupportRow]) -> str:
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


def _render_clade_table(rows: list[CladeTableRow]) -> str:
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


def _render_branch_stats(row: TreeBranchStatisticsRow) -> str:
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


def _render_list(items: list[str]) -> str:
    if not items:
        return "<p>none</p>"
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _write_tree_report_html(
    *,
    path: Path,
    title: str,
    figure_svg: str,
    reviewer_summary: list[str],
    limitations: list[str],
    support_rows: list[TreeSupportRow],
    clade_rows: list[CladeTableRow],
    branch_stats: TreeBranchStatisticsRow,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
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
      {_json_script(manifest)}
      <div class="cards">
        <div class="card"><div class="label">Tip Count</div><div class="value">{inspection.tip_count}</div></div>
        <div class="card"><div class="label">Tree Quality</div><div class="value">{inspection.tree_quality_score}</div></div>
        <div class="card"><div class="label">Supported Branches</div><div class="value">{sum(1 for row in support_rows if row.support is not None)}</div></div>
        <div class="card"><div class="label">Long Outliers</div><div class="value">{branch_stats.long_outlier_count}</div></div>
      </div>
      <section>
        <h2>Reviewer Summary</h2>
        {_render_list(reviewer_summary)}
        <div class="note">
          Support classes are reviewer-facing thresholds:
          <strong>strong</strong> at support fraction at least 0.95,
          <strong>moderate</strong> at least 0.80,
          and <strong>weak</strong> below 0.80.
        </div>
      </section>
      <section>
        <h2>Tree Image</h2>
        <div class="figure-frame">{figure_svg}</div>
      </section>
      <section>
        <h2>Support Table</h2>
        {_render_support_table(support_rows)}
      </section>
      <section>
        <h2>Clade Table</h2>
        {_render_clade_table(clade_rows)}
      </section>
      <section>
        <h2>Branch-Length Stats</h2>
        {_render_branch_stats(branch_stats)}
      </section>
      <section>
        <h2>Limitations</h2>
        {_render_list(limitations)}
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


def build_tree_report_package(
    tree_path: Path,
    *,
    out_dir: Path,
    title: str = "Bijux Full Tree Report",
) -> TreeReportPackageResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "tree-report.html"
    figure_path = out_dir / "tree-image.svg"
    support_table_path = out_dir / "support-table.tsv"
    clade_table_path = out_dir / "clade-table.tsv"
    branch_stats_path = out_dir / "branch-stats.tsv"
    manifest_path = out_dir / "tree-report.manifest.json"

    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    support_audit = audit_support_label_rendering(tree_path)
    figure = render_tree_svg(
        tree_path,
        out_path=figure_path,
        layout="phylogram",
        show_support_values=support_audit.validated,
        validated_support_labels=support_audit.labels_by_node,
        support_validation_warnings=support_audit.warnings,
    )
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    branch_lengths = analyze_branch_length_distribution(tree_path)
    branch_stats = summarize_tree_branch_statistics(branch_lengths)
    reviewer_summary, limitations = _reviewer_summary(
        inspection=inspection,
        support_rows=support_rows,
        branch_stats=branch_stats,
        support_audit=support_audit,
    )

    write_tree_support_table(support_table_path, support_rows)
    from bijux_phylogenetics.clades import write_clade_table

    write_clade_table(clade_table_path, clades)
    write_tree_branch_statistics_table(branch_stats_path, branch_stats)

    machine_manifest = {
        "report_kind": "tree_package",
        "title": title,
        "input_path": str(tree_path),
        "input_checksum": _checksum(tree_path),
        "outputs": {
            "report_path": str(report_path),
            "figure_path": str(figure_path),
            "support_table_path": str(support_table_path),
            "clade_table_path": str(clade_table_path),
            "branch_stats_path": str(branch_stats_path),
        },
        "metrics": {
            "tip_count": inspection.tip_count,
            "node_count": inspection.node_count,
            "clade_count": inspection.clade_count,
            "supported_branch_count": sum(
                1 for row in support_rows if row.support is not None
            ),
            "long_outlier_count": branch_stats.long_outlier_count,
            "short_outlier_count": branch_stats.short_outlier_count,
            "rendered_support_count": figure.rendered_support_count,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
        "validation": asdict(validation),
        "inspection": asdict(inspection),
        "forensic": asdict(forensic),
        "support_audit": asdict(support_audit),
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_tree_report_html(
        path=report_path,
        title=title,
        figure_svg=figure_path.read_text(encoding="utf-8"),
        reviewer_summary=reviewer_summary,
        limitations=limitations,
        support_rows=support_rows,
        clade_rows=clades.rows,
        branch_stats=branch_stats,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        manifest=machine_manifest,
    )
    return TreeReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        figure_path=figure_path,
        support_table_path=support_table_path,
        clade_table_path=clade_table_path,
        branch_stats_path=branch_stats_path,
        manifest_path=manifest_path,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        figure=figure,
        support_audit=support_audit,
        clades=clades,
        branch_lengths=branch_lengths,
        support_rows=support_rows,
        branch_stats=branch_stats,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
        machine_manifest=machine_manifest,
    )
