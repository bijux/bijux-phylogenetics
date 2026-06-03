from __future__ import annotations

from html import escape
from pathlib import Path

from .contracts import PublicationPackageRevalidationResult


def write_html_report(
    path: Path,
    *,
    result: PublicationPackageRevalidationResult,
    manifest_path: Path,
) -> Path:
    artifact_rows = "\n".join(
        (
            "      <tr>"
            f"<td>{escape(row.relative_path)}</td>"
            f"<td>{escape(row.artifact_scope)}</td>"
            f"<td>{escape(row.section)}</td>"
            f"<td>{escape(row.status)}</td>"
            f"<td>{escape(row.detail)}</td>"
            "</tr>"
        )
        for row in result.artifact_rows
    )
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
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Publication Package Revalidation</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: #f5f1ea; color: #1f2a24; }",
            "    main { max-width: 1120px; margin: 0 auto; padding: 28px; }",
            "    h1, h2 { margin: 0 0 12px; }",
            "    .panel { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 18px; margin-top: 18px; }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }",
            "    .card { background: #fffdf8; border: 1px solid #d7d0c3; border-radius: 14px; padding: 14px; }",
            "    .label { font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em; color: #676154; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 12px; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid #e6dfd1; vertical-align: top; }",
            "    code { background: #f0eadf; padding: 1px 4px; border-radius: 4px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Publication Package Revalidation</h1>",
            '  <div class="cards">',
            f'    <div class="card"><span class="label">report kind</span><span class="value">{escape(result.report_kind)}</span></div>',
            f'    <div class="card"><span class="label">overall status</span><span class="value">{escape(result.overall_revalidation_status)}</span></div>',
            f'    <div class="card"><span class="label">original artifacts match</span><span class="value">{str(result.all_original_artifacts_match).lower()}</span></div>',
            f'    <div class="card"><span class="label">unexpected files</span><span class="value">{result.unexpected_file_count}</span></div>',
            "  </div>",
            '  <section class="panel">',
            "    <h2>Package Root</h2>",
            f"    <p><code>{escape(str(result.package_root))}</code></p>",
            '    <h2 style="margin-top: 16px;">Manifest</h2>',
            f"    <p><code>{escape(str(manifest_path))}</code></p>",
            "  </section>",
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
            "      <thead><tr><th>Artifact</th><th>Scope</th><th>Section</th><th>Status</th><th>Detail</th></tr></thead>",
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
