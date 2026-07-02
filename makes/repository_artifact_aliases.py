#!/usr/bin/env python3
"""Materialize repository artifact aliases with repository-local root policy."""

from __future__ import annotations

import argparse
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from types import ModuleType
from typing import Callable, cast


AliasMaterializer = Callable[..., None]
PackageDirectoryDiscoverer = Callable[..., list[Path]]


def _shared_script_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "bijux-py"
        / "repository"
        / "artifact_aliases.py"
    )


def _load_shared_alias_module() -> ModuleType:
    script_path = _shared_script_path()
    spec = spec_from_file_location("bijux_shared_artifact_aliases", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"unable to load shared artifact alias script: {script_path}"
        )
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Materialize repository and package artifact alias links."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    root_parser = subparsers.add_parser(
        "root",
        help="materialize repository aliases and package aliases under packages/",
    )
    root_parser.add_argument("--repo-root", required=True, type=Path)
    root_parser.add_argument("--packages-dir", type=Path)

    package_parser = subparsers.add_parser(
        "package",
        help="materialize aliases for one package root",
    )
    package_parser.add_argument("--repo-root", required=True, type=Path)
    package_parser.add_argument("--package-dir", required=True, type=Path)

    return parser.parse_args()


def _materialize_root_aliases_preserving_local_venv(
    *, module: ModuleType, repo_root: Path
) -> None:
    root_alias_layout = cast(dict[str, Path], module.ROOT_ALIAS_LAYOUT)
    materialize_alias = cast(AliasMaterializer, module._materialize_alias)

    for alias_name, target_rel in root_alias_layout.items():
        link_path = repo_root / alias_name
        target_path = repo_root / target_rel
        if alias_name == ".venv" and link_path.exists() and not link_path.is_symlink():
            if not link_path.is_dir():
                raise RuntimeError(
                    f"refusing to preserve non-directory root alias path '{link_path}'"
                )
            target_path.mkdir(parents=True, exist_ok=True)
            continue
        materialize_alias(link_path=link_path, target_path=target_path)


def _materialize_root_and_package_aliases(
    *, module: ModuleType, repo_root: Path, packages_dir: Path
) -> None:
    discover_package_dirs = cast(
        PackageDirectoryDiscoverer,
        module._discover_package_dirs,
    )
    materialize_package_aliases = cast(
        AliasMaterializer,
        module._materialize_package_aliases,
    )

    _materialize_root_aliases_preserving_local_venv(
        module=module,
        repo_root=repo_root,
    )
    for package_dir in discover_package_dirs(packages_dir=packages_dir):
        materialize_package_aliases(repo_root=repo_root, package_dir=package_dir)


def main() -> int:
    args = _parse_args()
    module = _load_shared_alias_module()
    repo_root = args.repo_root.resolve()

    if args.command == "root":
        packages_dir = (
            args.packages_dir.resolve()
            if args.packages_dir is not None
            else repo_root / "packages"
        )
        _materialize_root_and_package_aliases(
            module=module,
            repo_root=repo_root,
            packages_dir=packages_dir,
        )
        return 0

    if args.command == "package":
        materialize_package_aliases = cast(
            AliasMaterializer,
            module._materialize_package_aliases,
        )
        materialize_package_aliases(
            repo_root=repo_root,
            package_dir=args.package_dir.resolve(),
        )
        return 0

    raise AssertionError(f"unsupported command: {args.command}")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
