from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = next(
    parent
    for parent in Path(__file__).resolve().parents
    if (parent / "packages").is_dir() and (parent / "configs").is_dir()
)
WORKFLOWS_DIR = REPO_ROOT / ".github" / "workflows"
EXPECTED_RELEASE_PACKAGES = {"bijux-phylogenetics", "phylogenetic"}
EXPECTED_VERIFY_PACKAGES = {
    "bijux-phylogenetics",
    "phylogenetic",
    "bijux-phylogenetics-dev",
}


def _workflow(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    assert isinstance(data, dict)
    workflow: dict[str, Any] = {}
    for key, value in data.items():
        normalized_key = "on" if key is True else key
        if isinstance(normalized_key, str):
            workflow[normalized_key] = value
    return workflow


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _matrix_include(job: dict[str, Any]) -> list[dict[str, Any]]:
    strategy = _as_dict(job.get("strategy"))
    matrix = _as_dict(strategy.get("matrix"))
    include = matrix.get("include", [])
    return include if isinstance(include, list) else []


def _release_var(name: str) -> str:
    for line in (
        (REPO_ROOT / ".github" / "release.env").read_text(encoding="utf-8").splitlines()
    ):
        if line.startswith(f"{name}="):
            value = line.split("=", 1)[1].strip()
            return value.strip("'")
    raise AssertionError(f"missing release env variable: {name}")


def test_release_matrices_include_canonical_and_alias_packages() -> None:
    for variable in (
        "BIJUX_RELEASE_BUILD_MATRIX_JSON",
        "BIJUX_PYPI_PACKAGE_MATRIX_JSON",
        "BIJUX_GHCR_RELEASE_PACKAGE_MATRIX_JSON",
    ):
        matrix = json.loads(_release_var(variable))
        slugs = {entry["package_slug"] for entry in matrix}
        assert slugs == EXPECTED_RELEASE_PACKAGES


def test_verify_workflow_checks_canonical_alias_and_dev_packages() -> None:
    workflow = _workflow(WORKFLOWS_DIR / "verify.yml")
    jobs = _as_dict(workflow.get("jobs"))
    package_job = _as_dict(jobs.get("package"))

    assert package_job.get("uses") == "./.github/workflows/ci.yml"

    include = _matrix_include(package_job)
    found = {entry["package_slug"] for entry in include if isinstance(entry, dict)}
    assert found == EXPECTED_VERIFY_PACKAGES

    runtime_entry = next(
        entry for entry in include if entry["package_slug"] == "bijux-phylogenetics"
    )
    assert runtime_entry["check_targets"] == (
        "[\"quality\", \"security\", \"api\", \"openapi-drift\", \"build\", \"sbom\"]"
    )
