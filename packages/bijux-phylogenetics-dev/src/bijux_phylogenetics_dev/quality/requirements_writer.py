"""Render package dependency groups into pip-audit requirements input."""

from __future__ import annotations

import argparse
from pathlib import Path
import tomllib
from typing import Any


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Write package dependency groups as requirements.txt content.",
    )
    parser.add_argument("--pyproject", required=True, help="Path to pyproject.toml.")
    parser.add_argument("--group", required=True, choices=("prod", "dev"))
    parser.add_argument(
        "--optional-group",
        action="append",
        default=[],
        help="Optional dependency group to include for the selected output.",
    )
    parser.add_argument("--output", required=True, help="Output requirements path.")
    return parser.parse_args(argv)


def _coerce_table(value: object) -> dict[str, Any]:
    """Return a TOML table when the value is one."""
    return value if isinstance(value, dict) else {}


def _coerce_list(value: object) -> list[str]:
    """Return a string list when the value is one."""
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, str)]


def load_dependency_lines(
    pyproject_path: Path, group: str, optional_groups: list[str]
) -> list[str]:
    """Load dependency lines from one package pyproject file."""
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    project_table = _coerce_table(pyproject.get("project"))
    requirements = list(_coerce_list(project_table.get("dependencies")))
    if group != "dev":
        return requirements
    optional_dependencies = _coerce_table(project_table.get("optional-dependencies"))
    for optional_group in optional_groups:
        requirements.extend(_coerce_list(optional_dependencies.get(optional_group)))
    return requirements


def companion_dev_pyproject(pyproject_path: Path) -> Path | None:
    """Resolve the sibling repository dev package when present."""
    package_dir = pyproject_path.parent
    if not package_dir.name.startswith("bijux-phylogenetics"):
        return None
    sibling = package_dir.parent / f"{package_dir.name}-dev" / "pyproject.toml"
    if sibling.is_file():
        return sibling
    return None


def unique_requirements(requirements: list[str]) -> list[str]:
    """Deduplicate requirement lines while preserving input order."""
    seen: set[str] = set()
    deduplicated: list[str] = []
    for requirement in requirements:
        if requirement in seen:
            continue
        seen.add(requirement)
        deduplicated.append(requirement)
    return deduplicated


def main(argv: list[str] | None = None) -> int:
    """Write requirements.txt content for the requested dependency groups."""
    args = parse_args(argv)
    pyproject_path = Path(args.pyproject).resolve()
    output_path = Path(args.output).resolve()

    requirements = load_dependency_lines(
        pyproject_path, args.group, list(args.optional_group)
    )
    sibling_dev_pyproject = companion_dev_pyproject(pyproject_path)
    if sibling_dev_pyproject is not None:
        requirements.extend(
            load_dependency_lines(
                sibling_dev_pyproject, args.group, list(args.optional_group)
            )
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "\n".join(unique_requirements(requirements)) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
