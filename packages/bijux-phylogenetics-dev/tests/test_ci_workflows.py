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
        '["quality", "security", "api", "openapi-drift", "build", "sbom"]'
    )


def test_evidence_governance_workflow_separates_contracts_and_cleanroom_jobs() -> None:
    workflow = _workflow(WORKFLOWS_DIR / "evidence-governance.yml")
    jobs = _as_dict(workflow.get("jobs"))

    evidence_contracts = _as_dict(jobs.get("evidence-contracts"))
    evidence_cleanroom = _as_dict(jobs.get("evidence-cleanroom"))

    assert evidence_contracts["runs-on"] == "ubuntu-latest"
    assert evidence_cleanroom["runs-on"] == "ubuntu-latest"
    contract_run_steps = [
        step["run"]
        for step in evidence_contracts.get("steps", [])
        if isinstance(step, dict) and isinstance(step.get("run"), str)
    ]
    cleanroom_run_steps = [
        step["run"]
        for step in evidence_cleanroom.get("steps", [])
        if isinstance(step, dict) and isinstance(step.get("run"), str)
    ]
    contract_text = "\n".join(contract_run_steps)
    cleanroom_text = "\n".join(cleanroom_run_steps)

    assert "make validate-evidence-book" in contract_text
    assert "make check-evidence-artifacts" in contract_text
    assert "make check-artifact-governance" in contract_text
    assert "make rerun-evidence-cleanroom EVIDENCE_STUDY_ID=primate-longevity-signal EVIDENCE_IDS=evidence-002" in cleanroom_text
    assert "make rerun-evidence-cleanroom EVIDENCE_STUDY_ID=primate-pgls-and-signal EVIDENCE_IDS=evidence-002" in cleanroom_text


def test_publish_readiness_workflow_keeps_report_and_gate_jobs_separate() -> None:
    workflow = _workflow(WORKFLOWS_DIR / "publish-readiness.yml")
    jobs = _as_dict(workflow.get("jobs"))

    report_job = _as_dict(jobs.get("publish-readiness-report"))
    gate_job = _as_dict(jobs.get("release-readiness-gate"))

    report_steps = [
        step["run"]
        for step in report_job.get("steps", [])
        if isinstance(step, dict) and isinstance(step.get("run"), str)
    ]
    gate_steps = [
        step["run"]
        for step in gate_job.get("steps", [])
        if isinstance(step, dict) and isinstance(step.get("run"), str)
    ]
    report_text = "\n".join(report_steps)
    gate_text = "\n".join(gate_steps)

    assert "make check-config-ssot" in report_text
    assert "make report-release-readiness" in report_text
    assert "make check-release-readiness" in gate_text
    assert gate_job["if"] == "${{ github.event_name == 'workflow_dispatch' && inputs.enforce_release_gate }}"
