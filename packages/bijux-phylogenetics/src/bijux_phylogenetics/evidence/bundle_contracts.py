from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .study_contracts import load_study_contract

REQUIRED_BUNDLE_LOCAL_ARTIFACTS = (
    "reference.R",
    "analysis.py",
    "checks.json",
    "report.md",
    "provenance.json",
)
ARTIFACT_JSON_FILENAMES = {"checks.json", "provenance.json"}
RESULTS_DIRNAME = "results"
REQUIRED_BUNDLE_RESULT_ARTIFACTS = (
    f"{RESULTS_DIRNAME}/README.md",
    f"{RESULTS_DIRNAME}/manifest.json",
)
RESULT_ARTIFACT_JSON_FILENAMES = {f"{RESULTS_DIRNAME}/manifest.json"}
PRIMARY_OUTPUT_EXCLUDED_FILENAMES = {
    "README.md",
    "manifest.json",
    "claims.json",
    "claim_verdicts.json",
    "reviewer-summary.json",
    "reviewer-summary.md",
    "inputs.manifest.json",
    *REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
}
LOCAL_ARTIFACT_PURPOSES = {
    "reference.R": "r-reference-program",
    "analysis.py": "python-analysis-program",
    "checks.json": "machine-check-contract",
    "report.md": "human-report",
    "provenance.json": "provenance-record",
}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def iter_bundle_roots(repo_root: Path) -> list[Path]:
    studies_root = repo_root / "evidence-book" / "studies"
    return sorted(path for path in studies_root.glob("*/evidence-*") if path.is_dir())


def _relative(path: Path, repo_root: Path) -> str:
    return path.relative_to(repo_root).as_posix()


def _primary_output_paths(bundle_root: Path, repo_root: Path) -> list[str]:
    output_paths: list[str] = []
    for child in sorted(bundle_root.iterdir()):
        if not child.is_file() or child.name in PRIMARY_OUTPUT_EXCLUDED_FILENAMES:
            continue
        output_paths.append(_relative(child, repo_root))
    return output_paths


def _bundle_manifest(bundle_root: Path) -> dict[str, Any]:
    return _load_json(bundle_root / "manifest.json")


def _study_support_scripts(study_root: Path, repo_root: Path) -> dict[str, Any]:
    reference_root = study_root / "reference"
    reference_r_scripts = [
        _relative(path, repo_root)
        for path in sorted(reference_root.glob("*.R"))
        if path.is_file()
    ]
    reference_python_scripts = [
        _relative(path, repo_root)
        for path in sorted(reference_root.glob("*.py"))
        if path.is_file()
    ]
    return {
        "build_script_path": None,
        "reference_r_scripts": reference_r_scripts,
        "reference_python_scripts": reference_python_scripts,
    }


def _source_basis_locators(bundle_manifest: dict[str, Any]) -> list[str]:
    locators: list[str] = []
    for entry in bundle_manifest.get("source_basis", []):
        if isinstance(entry, dict):
            locator = entry.get("locator")
            if isinstance(locator, str) and locator:
                locators.append(locator)
    return locators


def build_bundle_artifact_contract(repo_root: Path, bundle_root: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    bundle_root = bundle_root.resolve()
    study_root = bundle_root.parent
    study_manifest = load_study_contract(study_root)
    bundle_manifest = _bundle_manifest(bundle_root)
    support_scripts = _study_support_scripts(study_root, repo_root)
    provenance_path = study_root / "provenance"
    provenance_locators = [
        _relative(path, repo_root)
        for path in sorted(provenance_path.glob("*.json"))
        if path.is_file()
    ]
    results_root = bundle_root / RESULTS_DIRNAME
    return {
        "study_id": str(bundle_manifest["study_id"]),
        "evidence_id": str(bundle_manifest["evidence_id"]),
        "study_title": str(study_manifest["study_title"]),
        "evidence_title": str(bundle_manifest["evidence_title"]),
        "summary": str(bundle_manifest["summary"]),
        "comparison_mode": str(bundle_manifest["comparison_mode"]),
        "claim_ids": list(bundle_manifest.get("claim_ids", [])),
        "claim_tags": list(bundle_manifest.get("claim_tags", [])),
        "verdict": dict(bundle_manifest.get("verdict", {})),
        "limitations": list(bundle_manifest.get("limitations", [])),
        "source_intake_policy": str(study_manifest.get("source_intake_policy", "")),
        "study_provenance_locator": str(
            study_manifest.get("provenance_descriptor_locator", "")
        ),
        "study_dataset_registry_locator": str(
            study_manifest.get("dataset_registry_locator", "")
        ),
        "source_basis_locators": _source_basis_locators(bundle_manifest),
        "primary_output_paths": _primary_output_paths(bundle_root, repo_root),
        "bundle_relative_path": _relative(bundle_root, repo_root),
        "results_relative_path": _relative(results_root, repo_root),
        "build_script_path": support_scripts["build_script_path"],
        "reference_r_scripts": support_scripts["reference_r_scripts"],
        "reference_python_scripts": support_scripts["reference_python_scripts"],
        "bundle_provenance_paths": provenance_locators,
        "required_local_artifacts": list(REQUIRED_BUNDLE_LOCAL_ARTIFACTS),
        "required_result_artifacts": list(REQUIRED_BUNDLE_RESULT_ARTIFACTS),
        "local_artifact_purposes": dict(LOCAL_ARTIFACT_PURPOSES),
    }
