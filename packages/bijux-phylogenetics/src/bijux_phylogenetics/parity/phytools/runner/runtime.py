from __future__ import annotations

from importlib import metadata
import os
from pathlib import Path
import shutil

# Parity helpers invoke repository-owned reference commands under governed paths.
import subprocess  # nosec B404


def repository_root() -> Path:
    return Path(__file__).resolve().parents[7]


def phytools_runner_path() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "resources"
        / "reference"
        / "phytools_parity_runner.R"
    )


def failure_root() -> Path:
    return repository_root() / "artifacts" / "phytools-parity-failures"


def reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def bijux_commit() -> str | None:
    git_executable = shutil.which("git")
    if git_executable is None:
        return None
    result = subprocess.run(  # nosec B603
        [git_executable, "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None
