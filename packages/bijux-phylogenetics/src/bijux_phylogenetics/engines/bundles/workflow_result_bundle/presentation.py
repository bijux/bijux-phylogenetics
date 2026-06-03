from __future__ import annotations

from html import escape
import json
from pathlib import Path


def write_bundle_report(
    path: Path,
    *,
    workflow: str,
    source_manifest_path: Path,
    input_entries: list[dict[str, str]],
    output_entries: dict[str, str],
    step_manifest_entries: dict[str, str],
    notes: list[str],
    stage_fingerprints: object,
    embedded_manifest: dict[str, object],
) -> Path:
    input_lines = [
        f"{entry['label']}: {entry.get('relative_path', '[checksum only]')}"
        for entry in input_entries
    ]
    sections = [
        ("workflow", f"workflow: {workflow}\nsource_manifest: {source_manifest_path}"),
        ("inputs", "\n".join(input_lines) if input_lines else "none copied"),
        (
            "outputs",
            "\n".join(
                f"{label}: {path_text}" for label, path_text in output_entries.items()
            ),
        ),
        (
            "step-manifests",
            "\n".join(
                f"{label}: {path_text}"
                for label, path_text in step_manifest_entries.items()
            )
            if step_manifest_entries
            else "none",
        ),
    ]
    if stage_fingerprints is not None:
        sections.append(
            (
                "stage-fingerprints",
                json.dumps(stage_fingerprints, indent=2, sort_keys=True),
            )
        )
    if notes:
        sections.append(("notes", "\n".join(f"- {note}" for note in notes)))
    body = "\n".join(
        f"<section><h2>{escape(title)}</h2><pre>{escape(content)}</pre></section>"
        for title, content in sections
    )
    embedded_json = escape(
        json.dumps(embedded_manifest, indent=2, sort_keys=True).replace("</", "<\\/")
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(f"Bijux Workflow Result Bundle: {workflow}")}</title>
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
    <h1>{escape(f"Bijux Workflow Result Bundle: {workflow}")}</h1>
    <script id="bijux-report-manifest" type="application/json">{embedded_json}</script>
    {body}
  </main>
</body>
</html>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path


def render_bundle_readme(
    *,
    workflow: str,
    source_manifest_path: Path,
    includes_input_files: bool,
    missing_input_paths: list[Path],
    output_labels: list[str],
    step_labels: list[str],
) -> str:
    lines = [
        "# Bijux Workflow Result Bundle",
        "",
        f"Workflow: `{workflow}`",
        f"Source manifest: `{source_manifest_path}`",
        "",
        "## Contents",
        "",
        "- `bundle.manifest.json`: bundle inventory, checksums, and required workflow entries.",
        "- `workflow-config.json`: extracted workflow config for reruns and review.",
        "- `workflow-rerun.json`: bundle-local rerun ledger using the copied input files.",
        "- `reports/workflow-report.html`: reviewer-facing summary for the bundled workflow.",
        "- `manifests/`: copied workflow manifest plus any step manifests.",
        "- `inputs/`: copied workflow input files when they were still available at export time.",
        "- `outputs/final/`: copied reviewer-facing workflow outputs.",
        "- `outputs/engine-artifacts/`: copied native engine artifacts declared by the step manifests.",
        "",
        "## Review Notes",
        "",
        f"- copied input files present: `{includes_input_files}`",
        f"- missing source inputs at export time: `{len(missing_input_paths)}`",
        f"- final output labels: `{', '.join(output_labels)}`",
        f"- step manifests: `{', '.join(step_labels) if step_labels else 'none'}`",
        "",
    ]
    if missing_input_paths:
        lines.append("Missing input paths:")
        lines.extend(f"- `{path}`" for path in missing_input_paths)
        lines.append("")
    return "\n".join(lines)
