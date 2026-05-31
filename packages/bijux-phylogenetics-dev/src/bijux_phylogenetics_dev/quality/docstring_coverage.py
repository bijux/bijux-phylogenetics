"""Measure docstring coverage without external tooling dependencies."""

from __future__ import annotations

import argparse
import ast
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import Any


@dataclass(frozen=True)
class CoverageConfig:
    """Own the coverage rules applied to Python sources."""

    fail_under: float
    ignore_init_method: bool
    ignore_module: bool
    ignore_private: bool
    ignore_semiprivate: bool


@dataclass(frozen=True)
class FileCoverage:
    """Capture the coverage result for one Python file."""

    path: Path
    documented: int
    total: int

    @property
    def percentage(self) -> float:
        """Return the file coverage percentage."""
        if self.total == 0:
            return 100.0
        return (self.documented / self.total) * 100.0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Measure docstring coverage for Python source files.",
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help="Source files or directories to scan.",
    )
    parser.add_argument(
        "--pyproject",
        default="pyproject.toml",
        help="Path to the package pyproject.toml file.",
    )
    parser.add_argument(
        "--fail-under",
        type=float,
        default=None,
        help="Minimum required coverage percentage.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Accepted for compatibility with the previous tooling surface.",
    )
    return parser.parse_args(argv)


def _coerce_table(value: object) -> dict[str, Any]:
    """Return a TOML table when the value is one."""
    return value if isinstance(value, dict) else {}


def _coerce_bool(value: object, default: bool) -> bool:
    """Return a boolean option with a default fallback."""
    return value if isinstance(value, bool) else default


def _coerce_float(value: object, default: float) -> float:
    """Return a float option with a default fallback."""
    if isinstance(value, (int, float)):
        return float(value)
    return default


def load_config(
    pyproject_path: Path, fail_under_override: float | None
) -> CoverageConfig:
    """Load coverage policy from the package pyproject file."""
    pyproject = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))
    tool_table = _coerce_table(pyproject.get("tool"))
    interrogate_table = _coerce_table(tool_table.get("interrogate"))
    fail_under = _coerce_float(interrogate_table.get("fail-under"), 100.0)
    if fail_under_override is not None:
        fail_under = fail_under_override
    return CoverageConfig(
        fail_under=fail_under,
        ignore_init_method=_coerce_bool(
            interrogate_table.get("ignore-init-method"), False
        ),
        ignore_module=_coerce_bool(interrogate_table.get("ignore-module"), False),
        ignore_private=_coerce_bool(interrogate_table.get("ignore-private"), False),
        ignore_semiprivate=_coerce_bool(
            interrogate_table.get("ignore-semiprivate"), False
        ),
    )


def iter_python_files(paths: Iterable[Path]) -> list[Path]:
    """Collect Python files beneath the requested roots."""
    discovered: set[Path] = set()
    for path in paths:
        if path.is_file():
            if path.suffix == ".py":
                discovered.add(path.resolve())
            continue
        if not path.exists():
            raise SystemExit(f"path not found: {path}")
        for candidate in path.rglob("*.py"):
            if candidate.is_file():
                discovered.add(candidate.resolve())
    return sorted(discovered)


def should_count_name(name: str, config: CoverageConfig) -> bool:
    """Decide whether a definition name participates in coverage."""
    is_private = name.startswith("__") and not name.endswith("__")
    is_semiprivate = (
        name.startswith("_") and not name.startswith("__") and not name.endswith("__")
    )
    if config.ignore_private and is_private:
        return False
    return not (config.ignore_semiprivate and is_semiprivate)


def iter_documentable_nodes(
    body: Sequence[ast.stmt], config: CoverageConfig
) -> Iterable[ast.FunctionDef | ast.ClassDef]:
    """Yield documentable nodes in source order."""
    for statement in body:
        if isinstance(statement, (ast.FunctionDef, ast.ClassDef)):
            yield statement
            if isinstance(statement, ast.ClassDef) and should_count_name(
                statement.name, config
            ):
                yield from iter_documentable_nodes(statement.body, config)


def analyze_file(path: Path, config: CoverageConfig) -> FileCoverage:
    """Measure docstring coverage for one Python file."""
    module = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    documented = 0
    total = 0

    if not config.ignore_module:
        total += 1
        if ast.get_docstring(module) is not None:
            documented += 1

    for node in iter_documentable_nodes(module.body, config):
        if isinstance(node, ast.ClassDef):
            if not should_count_name(node.name, config):
                continue
        elif isinstance(node, ast.FunctionDef):
            if config.ignore_init_method and node.name == "__init__":
                continue
            if not should_count_name(node.name, config):
                continue
        total += 1
        if ast.get_docstring(node) is not None:
            documented += 1

    return FileCoverage(path=path, documented=documented, total=total)


def relative_display_path(path: Path, root: Path) -> str:
    """Render a stable relative path when possible."""
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def emit_report(
    records: Sequence[FileCoverage], config: CoverageConfig, root: Path
) -> int:
    """Write the text report and return the matching process status."""
    documented = sum(record.documented for record in records)
    total = sum(record.total for record in records)
    overall = 100.0 if total == 0 else (documented / total) * 100.0

    for record in records:
        print(
            f"FILE: {relative_display_path(record.path, root)} | {record.percentage:.1f}%"
        )

    passed = overall >= config.fail_under
    verdict = "PASSED" if passed else "FAILED"
    print(f"RESULT: {verdict} ({overall:.1f}%)")
    return 0 if passed else 1


def main(argv: Sequence[str] | None = None) -> int:
    """Run the docstring coverage measurement command."""
    args = parse_args(argv)
    root = Path.cwd().resolve()
    config = load_config(Path(args.pyproject).resolve(), args.fail_under)
    files = iter_python_files(Path(path).resolve() for path in args.paths)
    records = [analyze_file(path, config) for path in files]
    return emit_report(records, config, root)


if __name__ == "__main__":
    raise SystemExit(main())
