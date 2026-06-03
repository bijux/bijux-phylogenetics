from __future__ import annotations

from collections.abc import Sequence
from html import escape
import json
from pathlib import Path

SummaryValue = str | int | float | bool
SummaryMetric = tuple[str, SummaryValue]
ArtifactLink = tuple[str, str, str | None]


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, indent=2, sort_keys=True).replace("</", "<\\/")
    return f'<script id="bijux-report-manifest" type="application/json">{serialized}</script>'


def write_html_report(
    *,
    title: str,
    sections: list[tuple[str, str]],
    out_path: Path,
    embedded_json: dict[str, object] | None = None,
    summary_metrics: Sequence[SummaryMetric] | None = None,
    artifact_links: Sequence[ArtifactLink] | None = None,
) -> Path:
    """Write a simple standalone HTML report."""
    body = "\n".join(
        f"<section><h2>{escape(name)}</h2><pre>{escape(content)}</pre></section>"
        for name, content in sections
    )
    manifest = _json_script(embedded_json) if embedded_json is not None else ""
    summary_section = ""
    artifact_section = ""
    extra_styles = ""
    if summary_metrics:
        extra_styles += """
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 0.75rem;
    }
    .summary-card {
      background: #f6fffd;
      border: 1px solid rgba(15, 118, 110, 0.18);
      border-radius: 14px;
      padding: 0.9rem 1rem;
    }
    .summary-card dt {
      margin: 0;
      font: 600 0.82rem/1.2 "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      color: #476b67;
    }
    .summary-card dd {
      margin: 0.35rem 0 0;
      font: 700 1.2rem/1.25 var(--mono);
      color: var(--ink);
    }
"""
        summary_section = "\n".join(
            [
                '<section><h2>summary</h2><div class="summary-grid">',
                *[
                    (
                        '<dl class="summary-card">'
                        f"<dt>{escape(name)}</dt>"
                        f"<dd>{escape(str(value))}</dd>"
                        "</dl>"
                    )
                    for name, value in summary_metrics
                ],
                "</div></section>",
            ]
        )
    if artifact_links:
        extra_styles += """
    .artifact-list {
      margin: 0;
      padding-left: 1.25rem;
    }
    .artifact-list li + li {
      margin-top: 0.65rem;
    }
    .artifact-list a {
      color: var(--accent);
      font-family: var(--mono);
      text-decoration: none;
    }
    .artifact-list a:hover {
      text-decoration: underline;
    }
    .artifact-note {
      display: block;
      margin-top: 0.2rem;
      color: #51646b;
    }
"""
        artifact_items: list[str] = []
        for label, href, note in artifact_links:
            note_markup = (
                f'<span class="artifact-note">{escape(note)}</span>' if note else ""
            )
            artifact_items.append(
                "<li>"
                f"{escape(label)}: "
                f'<a href="{escape(href, quote=True)}">{escape(href)}</a>'
                f"{note_markup}"
                "</li>"
            )
        artifact_section = "\n".join(
            [
                '<section><h2>artifacts</h2><ul class="artifact-list">',
                *artifact_items,
                "</ul></section>",
            ]
        )
    styles = """
    :root {
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f8fafc;
      --panel: #ffffff;
      --accent: #0f766e;
      --rule: #d6dee8;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }
    body {
      margin: 0;
      padding: 2rem;
      background: linear-gradient(180deg, #eef6f4 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }
    main {
      max-width: 960px;
      margin: 0 auto;
      background: var(--panel);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 2rem;
      box-shadow: 0 24px 80px rgba(15, 118, 110, 0.08);
    }
    h1, h2 {
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }
    h1 {
      margin-top: 0;
      color: var(--accent);
    }
    section + section {
      margin-top: 1.5rem;
      padding-top: 1.5rem;
      border-top: 1px solid var(--rule);
    }
    pre {
      margin: 0;
      white-space: pre-wrap;
      word-break: break-word;
      font-family: var(--mono);
      background: #f6f8fa;
      border-radius: 12px;
      padding: 1rem;
    }"""
    if extra_styles:
        styles = f"{styles}\n{extra_styles.rstrip()}"
    styles = styles.lstrip("\n")
    if not summary_metrics and not artifact_links:
        html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
{styles}
  </style>
</head>
<body>
  <main>
    <h1>{escape(title)}</h1>
    {manifest}
    {body}
  </main>
</body>
</html>
"""
    else:
        main_content = "\n".join(
            block
            for block in (manifest, summary_section, artifact_section, body)
            if block
        )
        html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
{styles}
  </style>
</head>
<body>
  <main>
    <h1>{escape(title)}</h1>
    {main_content}
  </main>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return out_path
