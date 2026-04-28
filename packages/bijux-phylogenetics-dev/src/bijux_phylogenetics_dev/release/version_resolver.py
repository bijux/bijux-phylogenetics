"""Version resolver helpers."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys
import tomllib
from typing import Any, cast

from bijux_phylogenetics_dev.trusted_process import run_text


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Resolve a package version from pyproject metadata or package tags."
    )
    parser.add_argument("--pyproject", required=True, help="Path to pyproject.toml")
    parser.add_argument("--package-name", required=True, help="Package slug/tag prefix")
    return parser.parse_args()


def _pyproject_data(pyproject_path: Path) -> dict[str, object]:
    """Load and validate the pyproject table."""
    data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"invalid pyproject table in {pyproject_path}")
    return cast(dict[str, object], data)


def _resolve_hatch_version(pyproject_path: Path) -> str | None:
    """Resolve the version through Hatch when the project is VCS-versioned."""
    result = run_text(
        [sys.executable, "-m", "hatch", "version"],
        capture_output=True,
        check=False,
        cwd=pyproject_path.parent,
    )
    if result.returncode != 0:
        return None
    version = result.stdout.strip()
    return version or None


def _tag_glob(pyproject: dict[str, object], package_name: str) -> str:
    """Infer the matching tag glob from Hatch version settings."""
    tool_table = cast(dict[str, Any], pyproject.get("tool", {}))
    hatch_table = cast(dict[str, Any], tool_table.get("hatch", {}))
    hatch_version = cast(dict[str, Any], hatch_table.get("version", {}))
    tag_pattern = hatch_version.get("tag-pattern")
    if isinstance(tag_pattern, str):
        marker = "(?P<version>"
        if tag_pattern.startswith("^") and marker in tag_pattern:
            prefix = tag_pattern[1:].split(marker, 1)[0]
            if prefix:
                return f"{prefix}*"
    return f"{package_name}/v*"


def _git_executable() -> str:
    """Resolve the git executable path for trusted subprocess calls."""
    resolved = shutil.which("git")
    if resolved is None:
        raise SystemExit("git executable not found")
    return resolved


def resolve_version(pyproject_path: Path, package_name: str) -> str:
    """Resolve the package version from static metadata, Hatch, or tags."""
    pyproject = _pyproject_data(pyproject_path)
    project = cast(dict[str, Any], pyproject.get("project", {}))
    version = project.get("version")
    if isinstance(version, str) and version:
        return version

    hatch_version = _resolve_hatch_version(pyproject_path)
    if hatch_version:
        return hatch_version

    tag_process = run_text(
        [
            _git_executable(),
            "tag",
            "--sort=v:refname",
            "--list",
            _tag_glob(pyproject, package_name),
        ],
        capture_output=True,
        check=False,
    )
    tags = [line.strip() for line in tag_process.stdout.splitlines() if line.strip()]
    if tags:
        tag = tags[-1]
        if "/v" in tag:
            return tag.rsplit("/v", 1)[1]
        return tag
    return "0.0.0"


def main() -> int:
    """Run the command-line entry point."""
    args = parse_args()
    print(resolve_version(Path(args.pyproject), args.package_name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
