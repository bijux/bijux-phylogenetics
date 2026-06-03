from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

PACKAGE_ROOT = Path(__file__).resolve().parents[2]

REQUIRED_DISTRIBUTION_RESOURCES = [
    "bijux_phylogenetics/resources/examples/alignments/example_alignment.fasta",
    "bijux_phylogenetics/resources/examples/metadata/example_traits.tsv",
    "bijux_phylogenetics/resources/examples/trees/example_tree.nwk",
    "bijux_phylogenetics/resources/datasets/mammals/primate_comparative/tree.nwk",
    "bijux_phylogenetics/resources/reference/phytools_parity_runner.R",
]


def _clean_subprocess_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    return environment


def _build_python() -> str:
    candidate = PACKAGE_ROOT / ".venv" / "bin" / "python"
    if candidate.exists():
        return str(candidate)
    return sys.executable


def build_installable_distributions(output_root: Path) -> tuple[Path, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            _build_python(),
            "-m",
            "build",
            "--wheel",
            "--sdist",
            "--outdir",
            str(output_root),
        ],
        cwd=PACKAGE_ROOT,
        env=_clean_subprocess_environment(),
        check=True,
        capture_output=True,
        text=True,
    )
    wheel_paths = sorted(output_root.glob("*.whl"))
    sdist_paths = sorted(output_root.glob("*.tar.gz"))
    assert len(wheel_paths) == 1
    assert len(sdist_paths) == 1
    return wheel_paths[0], sdist_paths[0]


def assert_distribution_contains_packaged_resources(distribution_path: Path) -> None:
    if distribution_path.suffix == ".whl":
        with zipfile.ZipFile(distribution_path) as archive:
            members = archive.namelist()
    else:
        with tarfile.open(distribution_path, "r:gz") as archive:
            members = archive.getnames()
    missing = [
        resource
        for resource in REQUIRED_DISTRIBUTION_RESOURCES
        if not any(member.endswith(resource) for member in members)
    ]
    assert not missing, (
        f"missing packaged resources in {distribution_path.name}: {missing}"
    )


def create_clean_virtualenv(venv_root: Path) -> Path:
    subprocess.run(
        [sys.executable, "-m", "venv", str(venv_root)],
        cwd=PACKAGE_ROOT,
        env=_clean_subprocess_environment(),
        check=True,
        capture_output=True,
        text=True,
    )
    return venv_root / "bin" / "python"


def install_distribution(venv_python: Path, distribution_path: Path) -> None:
    subprocess.run(
        [
            str(venv_python),
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            str(distribution_path),
        ],
        cwd=PACKAGE_ROOT,
        env=_clean_subprocess_environment(),
        check=True,
        capture_output=True,
        text=True,
    )


def run_installed_cli(
    venv_root: Path,
    arguments: list[str],
    *,
    cwd: Path,
) -> subprocess.CompletedProcess[str]:
    executable = venv_root / "bin" / "bijux-phylogenetics"
    return subprocess.run(
        [str(executable), *arguments],
        cwd=cwd,
        env=_clean_subprocess_environment(),
        check=True,
        capture_output=True,
        text=True,
    )


def copy_installed_example_inputs(
    venv_python: Path, destination: Path
) -> dict[str, Path]:
    completed = subprocess.run(
        [
            str(venv_python),
            "-c",
            (
                "from pathlib import Path; "
                "import json; "
                "from bijux_phylogenetics.core import copy_example_inputs; "
                "copied = copy_example_inputs(Path(__import__('sys').argv[1])); "
                "print(json.dumps({name: str(path) for name, path in copied.items()}, sort_keys=True))"
            ),
            str(destination),
        ],
        cwd=destination.parent,
        env=_clean_subprocess_environment(),
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    return {name: Path(path) for name, path in payload.items()}
