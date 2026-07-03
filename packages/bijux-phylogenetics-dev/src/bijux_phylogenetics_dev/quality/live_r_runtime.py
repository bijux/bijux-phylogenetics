from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

DEFAULT_LIVE_R_TEST_PACKAGES = (
    "ape",
    "jsonlite",
    "phylolm",
    "phytools",
    "geiger",
    "nlme",
)
DEFAULT_CRAN_REPOSITORY = "https://cloud.r-project.org"


def repository_r_library(repo_root: Path) -> Path:
    return repo_root / "artifacts" / "r-lib"


def _quote_r_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _package_vector_expression(package_names: tuple[str, ...]) -> str:
    return "c(" + ", ".join(_quote_r_string(name) for name in package_names) + ")"


def _base_environment(repo_root: Path) -> dict[str, str]:
    environment = dict(os.environ)
    environment.setdefault("R_LIBS_USER", str(repository_r_library(repo_root)))
    return environment


def _resolve_rscript_executable(rscript_executable: str | None) -> str:
    if rscript_executable:
        return rscript_executable
    resolved = shutil.which("Rscript")
    if resolved is None:
        raise RuntimeError(
            "Rscript is required for full-suite live R validation but was not found on PATH."
        )
    return resolved


def _run_rscript(
    repo_root: Path,
    *,
    expression: str,
    rscript_executable: str,
    environment: dict[str, str],
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [rscript_executable, "-e", expression],
        capture_output=True,
        check=False,
        cwd=repo_root,
        env=environment,
        text=True,
    )


def query_live_r_test_packages(
    repo_root: Path,
    *,
    package_names: tuple[str, ...] = DEFAULT_LIVE_R_TEST_PACKAGES,
    rscript_executable: str | None = None,
) -> dict[str, bool]:
    resolved_rscript = _resolve_rscript_executable(rscript_executable)
    r_library = repository_r_library(repo_root)
    r_library.mkdir(parents=True, exist_ok=True)
    environment = _base_environment(repo_root)
    expression = "\n".join(
        (
            f"lib <- {_quote_r_string(str(r_library))}",
            "dir.create(lib, recursive = TRUE, showWarnings = FALSE)",
            ".libPaths(c(lib, .libPaths()))",
            f"packages <- {_package_vector_expression(package_names)}",
            (
                'for (package_name in packages) cat(package_name, "\\t", '
                'requireNamespace(package_name, quietly = TRUE), "\\n", sep = "")'
            ),
        )
    )
    completed = _run_rscript(
        repo_root,
        expression=expression,
        rscript_executable=resolved_rscript,
        environment=environment,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Unable to inspect live R test package availability.\n"
            + completed.stderr.strip()
        )
    availability: dict[str, bool] = {}
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        package_name, _, availability_value = line.partition("\t")
        availability[package_name] = availability_value.strip() == "TRUE"
    missing_from_output = [name for name in package_names if name not in availability]
    if missing_from_output:
        raise RuntimeError(
            "R package availability inspection did not report every requested package: "
            + ", ".join(missing_from_output)
        )
    return availability


def install_live_r_test_packages(
    repo_root: Path,
    *,
    package_names: tuple[str, ...],
    cran_repository: str = DEFAULT_CRAN_REPOSITORY,
    rscript_executable: str | None = None,
) -> None:
    resolved_rscript = _resolve_rscript_executable(rscript_executable)
    r_library = repository_r_library(repo_root)
    r_library.mkdir(parents=True, exist_ok=True)
    environment = _base_environment(repo_root)
    expression = "\n".join(
        (
            f"lib <- {_quote_r_string(str(r_library))}",
            f"repo <- {_quote_r_string(cran_repository)}",
            "dir.create(lib, recursive = TRUE, showWarnings = FALSE)",
            ".libPaths(c(lib, .libPaths()))",
            f"packages <- {_package_vector_expression(package_names)}",
            "options(repos = c(CRAN = repo))",
            (
                "install.packages(packages, lib = lib, quiet = FALSE, "
                "Ncpus = max(1L, parallel::detectCores(logical = FALSE)))"
            ),
        )
    )
    completed = _run_rscript(
        repo_root,
        expression=expression,
        rscript_executable=resolved_rscript,
        environment=environment,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Installing live R test packages failed.\n"
            + completed.stderr.strip()
            + ("\n" + completed.stdout.strip() if completed.stdout.strip() else "")
        )


def ensure_live_r_test_packages(
    repo_root: Path,
    *,
    package_names: tuple[str, ...] = DEFAULT_LIVE_R_TEST_PACKAGES,
    cran_repository: str = DEFAULT_CRAN_REPOSITORY,
    rscript_executable: str | None = None,
) -> dict[str, object]:
    resolved_rscript = _resolve_rscript_executable(rscript_executable)
    initial_availability = query_live_r_test_packages(
        repo_root,
        package_names=package_names,
        rscript_executable=resolved_rscript,
    )
    missing_before = tuple(
        package_name
        for package_name in package_names
        if not initial_availability[package_name]
    )
    installed_now: tuple[str, ...] = ()
    if missing_before:
        install_live_r_test_packages(
            repo_root,
            package_names=missing_before,
            cran_repository=cran_repository,
            rscript_executable=resolved_rscript,
        )
        installed_now = missing_before
    final_availability = query_live_r_test_packages(
        repo_root,
        package_names=package_names,
        rscript_executable=resolved_rscript,
    )
    missing_after = tuple(
        package_name
        for package_name in package_names
        if not final_availability[package_name]
    )
    if missing_after:
        raise RuntimeError(
            "Full-suite live R validation is still missing required packages after installation: "
            + ", ".join(missing_after)
        )
    return {
        "status": "ready",
        "rscript_executable": resolved_rscript,
        "cran_repository": cran_repository,
        "r_library": str(repository_r_library(repo_root)),
        "required_packages": list(package_names),
        "already_available": [
            package_name for package_name in package_names if package_name not in installed_now
        ],
        "installed_now": list(installed_now),
    }


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ensure the governed live R test packages are available for full-suite pytest lanes."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        required=True,
        help="Repository root that owns artifacts/r-lib.",
    )
    parser.add_argument(
        "--cran-repository",
        default=DEFAULT_CRAN_REPOSITORY,
        help="CRAN repository used when packages need installation.",
    )
    parser.add_argument(
        "--rscript-executable",
        default=None,
        help="Override the Rscript executable used for package inspection and installation.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional JSON report path.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_argument_parser()
    args = parser.parse_args(argv)
    report = ensure_live_r_test_packages(
        args.repo_root.resolve(),
        cran_repository=args.cran_repository,
        rscript_executable=args.rscript_executable,
    )
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(f"live R test runtime ready: {', '.join(report['required_packages'])}")
    if report["installed_now"]:
        print("installed now: " + ", ".join(report["installed_now"]))
    else:
        print("installed now: none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
