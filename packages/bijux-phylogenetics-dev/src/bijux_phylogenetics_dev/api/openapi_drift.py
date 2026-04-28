"""OpenAPI drift checks for repository API contracts."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys
from typing import Any

import yaml

from bijux_phylogenetics_dev.trusted_process import run_text

Schema = dict[str, Any]


def _as_schema(value: object) -> Schema:
    """Coerce to schema."""
    return value if isinstance(value, dict) else {}


def _load_schema(text: str) -> Schema:
    """Load schema."""
    return _as_schema(yaml.safe_load(text))


def _git_show(repo_root: Path, path: str) -> str | None:
    """Handle Git show."""
    git_bin = shutil.which("git")
    if git_bin is None:
        return None
    path_parts = Path(path).parts
    if Path(path).is_absolute() or ".." in path_parts:
        return None
    try:
        completed = run_text(
            [git_bin, "-C", str(repo_root), "show", f"HEAD~1:{path}"],
            check=True,
            capture_output=True,
        )
    except (FileNotFoundError, ValueError):
        return None
    except Exception:
        return None
    return completed.stdout


def _extract_fields(schema: Schema) -> set[str]:
    """Extract fields."""
    fields: set[str] = set()
    components = _as_schema(schema.get("components"))
    schemas = _as_schema(components.get("schemas"))
    for name, payload in schemas.items():
        props = _as_schema(_as_schema(payload).get("properties"))
        for prop in props:
            fields.add(f"{name}.{prop}")
    return fields


def run(repo_root: Path) -> int:
    """Check for breaking OpenAPI field removals without a schema version bump."""
    schema_paths = sorted((repo_root / "apis").glob("*/v1/schema.yaml"))
    if not schema_paths:
        print("No OpenAPI schemas found; skipping.")
        return 0

    failures: list[str] = []
    for schema_path in schema_paths:
        rel = schema_path.relative_to(repo_root).as_posix()
        current_schema = _load_schema(schema_path.read_text(encoding="utf-8"))
        previous_text = _git_show(repo_root, rel)
        if not previous_text:
            continue
        previous_schema = _load_schema(previous_text)
        removed_fields = _extract_fields(previous_schema) - _extract_fields(
            current_schema
        )
        if not removed_fields:
            continue
        current_version = (current_schema.get("info") or {}).get("version")
        previous_version = (previous_schema.get("info") or {}).get("version")
        if current_version == previous_version:
            failures.append(f"{rel}: breaking change without version bump")
            for field in sorted(removed_fields):
                failures.append(f"{rel}: removed {field}")

    if failures:
        print("Breaking OpenAPI changes detected:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the OpenAPI drift check."""
    parser = argparse.ArgumentParser(
        description="Detect breaking OpenAPI changes without version bumps."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root that contains the apis/ contract tree.",
    )
    return parser.parse_args()


def main() -> int:
    """Execute the OpenAPI drift check and return its process exit code."""
    args = parse_args()
    return run(args.repo_root.resolve())


if __name__ == "__main__":
    sys.exit(main())
