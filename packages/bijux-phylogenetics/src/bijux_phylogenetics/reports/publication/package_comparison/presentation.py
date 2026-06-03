from __future__ import annotations

from html import escape
from pathlib import Path

from .contracts import PublicationPackageComparisonResult


def write_html_report(
    path: Path,
    *,
    result: PublicationPackageComparisonResult,
) -> Path:
    """Write the reviewer-facing package comparison HTML report."""
    check_rows = "\n".join(
        (
            "      <tr>"
            f"<td>{escape(row.section)}</td>"
            f"<td>{escape(row.check_id)}</td>"
            f"<td>{escape(row.status)}</td>"
            f"<td>{escape(row.summary)}</td>"
            f"<td>{escape(row.evidence)}</td>"
            "</tr>"
        )
        for row in result.check_rows
    )
    artifact_rows = "\n".join(
        (
            "      <tr>"
            f"<td>{escape(row.relative_path)}</td>"
            f"<td>{escape(row.kind)}</td>"
            f"<td>{escape(row.status)}</td>"
            f"<td>{escape(row.detail)}</td>"
            "</tr>"
        )
        for row in result.artifact_rows
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Publication Package Comparison</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: #f5f1ea; color: #1f2a24; }",
            "    main { max-width: 1120px; margin: 0 auto; padding: 28px; }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }",
            "    .card { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 14px; }",
            "    .label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em; color: #676154; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    .panel { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 18px; margin-top: 18px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 12px; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid #e6dfd1; vertical-align: top; }",
            "    code { background: #f0eadf; padding: 1px 4px; border-radius: 4px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Publication Package Comparison</h1>",
            '  <div class="cards">',
            f'    <div class="card"><span class="label">dataset id</span><span class="value">{escape(result.dataset_id)}</span></div>',
            f'    <div class="card"><span class="label">overall status</span><span class="value">{escape(result.overall_comparison_status)}</span></div>',
            f'    <div class="card"><span class="label">changed artifacts</span><span class="value">{result.changed_artifact_count}</span></div>',
            f'    <div class="card"><span class="label">finding differences</span><span class="value">{result.scientific_finding_difference_count}</span></div>',
            "  </div>",
            '  <section class="panel">',
            "    <h2>Checks</h2>",
            "    <table>",
            "      <thead><tr><th>Section</th><th>Check</th><th>Status</th><th>Summary</th><th>Evidence</th></tr></thead>",
            "      <tbody>",
            check_rows,
            "      </tbody>",
            "    </table>",
            "  </section>",
            '  <section class="panel">',
            "    <h2>Artifacts</h2>",
            "    <table>",
            "      <thead><tr><th>Artifact</th><th>Kind</th><th>Status</th><th>Detail</th></tr></thead>",
            "      <tbody>",
            artifact_rows,
            "      </tbody>",
            "    </table>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path
