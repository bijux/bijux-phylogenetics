"""Run deptry with repository-owned configuration merged into a package scan."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import sys
import tempfile
import tomllib
from typing import Any

from bijux_phylogenetics_dev.trusted_process import run_text

TomlTable = dict[str, Any]


def _as_table(value: object) -> TomlTable:
    """Coerce to table."""
    return value if isinstance(value, dict) else {}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the deptry scan wrapper."""
    parser = argparse.ArgumentParser(
        description="Run deptry with repository-owned config merged into a package pyproject.toml.",
    )
    parser.add_argument(
        "--config", required=True, help="Path to the repo-owned deptry TOML file."
    )
    parser.add_argument(
        "--project-dir",
        required=True,
        help="Path to the package root containing pyproject.toml.",
    )
    parser.add_argument("--deptry-bin", default="", help="Deptry executable to invoke.")
    parser.add_argument(
        "roots",
        nargs="*",
        default=["."],
        help="Roots to scan, relative to the project dir.",
    )
    return parser.parse_args()


def resolve_relative_command(command: list[str], project_dir: Path) -> list[str]:
    """Resolve an executable path relative to the package directory when needed."""
    executable = Path(command[0]).expanduser()
    if executable.is_absolute():
        return [os.fspath(executable), *command[1:]]
    if len(executable.parts) == 1:
        resolved = shutil.which(command[0])
        if resolved is None:
            raise SystemExit(f"deptry executable not found: {command[0]}")
        return [resolved, *command[1:]]

    candidates = (
        project_dir / executable,
        project_dir.parent / executable,
        project_dir.parent.parent / executable,
    )
    for candidate in candidates:
        if candidate.is_file():
            return [os.fspath(candidate.resolve()), *command[1:]]
    raise SystemExit(f"deptry executable not found: {command[0]}")


def resolve_deptry_command(deptry_bin: str, project_dir: Path) -> list[str]:
    """Resolve the deptry command, defaulting to the active Python module entrypoint."""
    deptry_command = shlex.split(deptry_bin)
    if not deptry_command:
        return resolve_relative_command([sys.executable, "-m", "deptry"], project_dir)
    return resolve_relative_command(deptry_command, project_dir)


def merge_deptry_config(
    config_path: Path, package_slug: str, package_pyproject: TomlTable
) -> TomlTable:
    """Merge repository deptry defaults with package-specific overrides."""
    root_config = tomllib.loads(config_path.read_text(encoding="utf-8"))
    tool_table = _as_table(root_config.get("tool"))
    base_config = _as_table(tool_table.get("deptry"))
    package_configs = _as_table(
        _as_table(tool_table.get("repo_deptry")).get("packages")
    )
    package_override = _as_table(package_configs.get(package_slug))

    merged_config: TomlTable = dict(base_config)
    for key, value in package_override.items():
        current = merged_config.get(key)
        if isinstance(current, list) and isinstance(value, list):
            merged_config[key] = list(dict.fromkeys([*current, *value]))
        elif isinstance(current, dict) and isinstance(value, dict):
            merged_config[key] = {**current, **value}
        else:
            merged_config[key] = value

    if "known_first_party" not in merged_config:
        merged_config["known_first_party"] = [package_slug.replace("-", "_")]
    optional_dependencies = _as_table(
        _as_table(package_pyproject.get("project")).get("optional-dependencies")
    )
    available_groups = set(optional_dependencies)
    for key in ("pep621_dev_dependency_groups", "optional_dependencies_dev_groups"):
        dev_groups = merged_config.get(key)
        if not isinstance(dev_groups, list):
            continue
        filtered_groups = [group for group in dev_groups if group in available_groups]
        if filtered_groups:
            merged_config[key] = filtered_groups
        else:
            merged_config.pop(key, None)
    return merged_config


def render_deptry_config(config: dict[str, object]) -> str:
    """Render the merged deptry configuration as TOML text."""
    lines = ["[tool.deptry]"]
    config_copy = dict(config)
    package_module_name_map = config_copy.pop("package_module_name_map", None)
    for key, value in config_copy.items():
        lines.append(f"{key} = {render_toml_value(value)}")
    if isinstance(package_module_name_map, dict) and package_module_name_map:
        lines.append("")
        lines.append("[tool.deptry.package_module_name_map]")
        for key, value in package_module_name_map.items():
            lines.append(f"{json.dumps(str(key))} = {json.dumps(str(value))}")
    return "\n".join(lines)


def render_toml_value(value: object) -> str:
    """Render a Python value into the limited TOML surface used by deptry config."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, list):
        return "[" + ", ".join(render_toml_value(item) for item in value) + "]"
    if isinstance(value, dict):
        items = ", ".join(
            f"{key} = {render_toml_value(item)}" for key, item in value.items()
        )
        return "{ " + items + " }"
    return str(value)


def main() -> int:
    """Run deptry against a temporary merged pyproject configuration."""
    args = parse_args()
    project_dir = Path(args.project_dir).resolve()
    config_path = Path(args.config).resolve()
    pyproject_path = project_dir / "pyproject.toml"

    if not pyproject_path.is_file():
        raise SystemExit(f"pyproject.toml not found in {project_dir}")
    if not config_path.is_file():
        raise SystemExit(f"deptry config not found at {config_path}")

    pyproject_text = pyproject_path.read_text(encoding="utf-8").rstrip()
    package_pyproject = tomllib.loads(pyproject_text)
    deptry_command = resolve_deptry_command(args.deptry_bin, project_dir)
    config_text = render_deptry_config(
        merge_deptry_config(config_path, project_dir.name, package_pyproject)
    )
    merged_text = f"{pyproject_text}\n\n{config_text}\n"

    with tempfile.TemporaryDirectory(prefix="deptry-") as tmpdir:
        merged_pyproject = Path(tmpdir) / "pyproject.toml"
        merged_pyproject.write_text(merged_text, encoding="utf-8")
        completed = run_text(
            [*deptry_command, "--config", os.fspath(merged_pyproject), *args.roots],
            cwd=project_dir,
            check=False,
        )
        return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
