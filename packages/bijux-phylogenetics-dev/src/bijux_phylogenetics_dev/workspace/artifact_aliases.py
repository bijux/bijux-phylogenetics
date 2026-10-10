"""Materialize repository artifact aliases from the shared governed layout."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from types import ModuleType
from typing import cast

AliasMaterializer = Callable[..., None]
PackageDirectoryDiscoverer = Callable[..., list[Path]]


def _repository_root() -> Path:
    current_path = Path(__file__).resolve()
    for parent in current_path.parents:
        if (
            parent / "makes" / "bijux-py" / "repository" / "artifact_aliases.py"
        ).is_file():
            return parent
    raise RuntimeError(f"unable to locate repository root from {current_path}")


def _shared_script_path() -> Path:
    return (
        _repository_root() / "makes" / "bijux-py" / "repository" / "artifact_aliases.py"
    )


def _load_shared_alias_module() -> ModuleType:
    script_path = _shared_script_path()
    spec = spec_from_file_location("bijux_shared_artifact_aliases", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"unable to load shared artifact alias script: {script_path}"
        )
    module = module_from_spec(spec)
    # Dataclass annotation resolution needs its owning module registered during import.
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(spec.name, None)
        raise
    return module


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
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

    return parser.parse_args(argv)


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
    materialize_root_aliases = cast(
        AliasMaterializer,
        module._materialize_root_aliases,
    )
    materialize_root_aliases(repo_root=repo_root)
    for package_dir in discover_package_dirs(packages_dir=packages_dir):
        materialize_package_aliases(repo_root=repo_root, package_dir=package_dir)


def main(argv: Sequence[str] | None = None) -> int:
    """Materialize governed artifact aliases for the repository or one package."""
    args = _parse_args(argv)
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
