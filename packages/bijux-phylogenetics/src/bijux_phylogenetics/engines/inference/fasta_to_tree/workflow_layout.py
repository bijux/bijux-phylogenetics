from __future__ import annotations

from pathlib import Path


def _write_tsv(path: Path, *, header: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.read_bytes())
    return destination


def _display_path(path: Path, *, root_dir: Path | None) -> str:
    if root_dir is None:
        return str(path)
    try:
        return str(path.relative_to(root_dir))
    except ValueError:
        return str(path)


def _display_command(
    command: list[str],
    *,
    root_dir: Path | None,
    aliases: dict[str, str] | None = None,
) -> str:
    rendered: list[str] = []
    token_aliases = {} if aliases is None else aliases
    for token in command:
        if token in token_aliases:
            rendered.append(token_aliases[token])
            continue
        if token.startswith("/"):
            rendered.append(_display_path(Path(token), root_dir=root_dir))
        else:
            rendered.append(token)
    return " ".join(rendered)


def _artifact_prefix(out_dir: Path, prefix: str, step_name: str) -> Path:
    return out_dir / "engine-artifacts" / prefix / step_name / step_name


def _final_output_paths(out_dir: Path, prefix: str) -> dict[str, Path]:
    root = out_dir / prefix
    return {
        "alignment": root.with_suffix(".aln"),
        "trimmed_alignment": root.with_suffix(".trimmed.aln"),
        "tree": root.with_suffix(".tree"),
        "log": root.with_suffix(".log"),
        "methods_summary": root.with_suffix(".methods-summary.md"),
        "model_table": root.with_suffix(".model.tsv"),
        "support_table": root.with_suffix(".support.tsv"),
        "manifest": root.with_suffix(".manifest.json"),
        "run_manifest": root.with_suffix(".run.json"),
    }
