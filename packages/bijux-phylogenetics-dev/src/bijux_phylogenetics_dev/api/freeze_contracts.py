"""Repository API freeze checks for schema, pinned JSON, and hash digests."""

from __future__ import annotations

import argparse
from datetime import date, datetime
from enum import Enum
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import yaml


def _load_artifact(path: Path) -> Any:
    """Load artifact."""
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".json":
        return json.loads(text)
    return yaml.safe_load(text)


def _canonicalize(payload: Any) -> Any:
    """Handle canonicalize."""
    if isinstance(payload, dict):
        return {
            str(key): _canonicalize(value)
            for key, value in sorted(payload.items(), key=lambda item: str(item[0]))
        }
    if isinstance(payload, list):
        return [_canonicalize(item) for item in payload]
    if isinstance(payload, tuple):
        return [_canonicalize(item) for item in payload]
    if isinstance(payload, (datetime, date)):
        return payload.isoformat()
    if isinstance(payload, Enum):
        return payload.value
    if isinstance(payload, Path):
        return str(payload)
    return payload


def _extract_hash_value(path: Path) -> str | None:
    """Extract hash value."""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("sha256:"):
            return line.split(":", 1)[1].strip()
    return None


def run(repo_root: Path) -> int:
    """Validate pinned API freeze artifacts for every checked-in schema."""
    failures: list[str] = []
    schema_paths = sorted((repo_root / "apis").glob("*/v1/schema.yaml"))
    if not schema_paths:
        print("No OpenAPI schemas found under apis/*/v1/schema.yaml", file=sys.stderr)
        return 1

    for schema_path in schema_paths:
        package_dir = schema_path.parent.parent.name
        pinned_path = schema_path.with_name("pinned_openapi.json")
        hash_path = schema_path.with_name("schema.hash")
        if not pinned_path.exists():
            failures.append(f"{package_dir}: missing pinned_openapi.json")
            continue
        if not hash_path.exists():
            failures.append(f"{package_dir}: missing schema.hash")
            continue

        expected = _canonicalize(_load_artifact(schema_path))
        actual = _canonicalize(_load_artifact(pinned_path))
        if expected != actual:
            failures.append(
                f"{package_dir}: pinned_openapi.json does not match schema.yaml"
            )

        digest = hashlib.sha256(
            schema_path.read_text(encoding="utf-8").encode("utf-8")
        ).hexdigest()
        schema_digest = _extract_hash_value(hash_path)
        if schema_digest != digest:
            failures.append(f"{package_dir}: schema.hash does not match schema.yaml")

    if failures:
        print("API freeze contract violations detected:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the API freeze validator."""
    parser = argparse.ArgumentParser(
        description="Validate checked-in API freeze contracts for a repository root."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root that contains the apis/ contract tree.",
    )
    return parser.parse_args()


def main() -> int:
    """Execute the API freeze validator and return its process exit code."""
    args = parse_args()
    return run(args.repo_root.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
