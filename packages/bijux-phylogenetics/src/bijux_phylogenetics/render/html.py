from __future__ import annotations

import json
from html import escape
from pathlib import Path


def _json_script(payload: dict[str, object]) -> str:
    serialized = json.dumps(payload, indent=2, sort_keys=True).replace("</", "<\\/")
    return f'<script id="bijux-report-manifest" type="application/json">{serialized}</script>'


def write_html_report(
    *,
    title: str,
    sections: list[tuple[str, str]],
    out_path: Path,
    embedded_json: dict[str, object] | None = None,
) -> Path:
    """Write a simple standalone HTML report."""
    body = "\n".join(
        f"<section><h2>{escape(name)}</h2><pre>{escape(content)}</pre></section>"
        for name, content in sections
    )
    manifest = _json_script(embedded_json) if embedded_json is not None else ""
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f8fafc;
      --panel: #ffffff;
      --accent: #0f766e;
      --rule: #d6dee8;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
    body {{
      margin: 0;
      padding: 2rem;
      background: linear-gradient(180deg, #eef6f4 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 960px;
      margin: 0 auto;
      background: var(--panel);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 2rem;
      box-shadow: 0 24px 80px rgba(15, 118, 110, 0.08);
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin-top: 0;
      color: var(--accent);
    }}
    section + section {{
      margin-top: 1.5rem;
      padding-top: 1.5rem;
      border-top: 1px solid var(--rule);
    }}
    pre {{
      margin: 0;
      white-space: pre-wrap;
      word-break: break-word;
      font-family: var(--mono);
      background: #f6f8fa;
      border-radius: 12px;
      padding: 1rem;
    }}
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
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return out_path
