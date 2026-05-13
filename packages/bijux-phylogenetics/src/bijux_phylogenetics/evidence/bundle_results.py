from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .bundle_contracts import RESULTS_DIRNAME

RESULTS_README_RELATIVE_PATH = f"{RESULTS_DIRNAME}/README.md"
RESULTS_MANIFEST_RELATIVE_PATH = f"{RESULTS_DIRNAME}/manifest.json"


def render_results_manifest(contract: dict[str, Any]) -> dict[str, Any]:
    results_relative_path = str(contract["results_relative_path"])
    primary_output_paths = list(contract["primary_output_paths"])
    return {
        "schema_version": 1,
        "study_id": contract["study_id"],
        "evidence_id": contract["evidence_id"],
        "results_directory": results_relative_path,
        "governed_primary_output_count": len(primary_output_paths),
        "governed_primary_outputs": primary_output_paths,
        "local_execution_products": [
            f"{results_relative_path}/reference-contract.json",
            f"{results_relative_path}/analysis-run.json",
        ],
        "bundle_local_artifacts": [
            f"{contract['bundle_relative_path']}/{filename}"
            for filename in contract["required_local_artifacts"]
        ],
        "pass_when": [
            "bundle-local authored artifacts remain separate from machine outputs",
            "governed primary outputs stay inventoried explicitly for reviewers",
            "local evidence scripts can write ephemeral execution products inside the evidence-owned results directory",
        ],
    }


def render_results_readme(contract: dict[str, Any]) -> str:
    primary_output_paths = list(contract["primary_output_paths"])
    lines = [
        f"# Results for {contract['evidence_id']}",
        "",
        "This directory is the evidence-owned execution workspace for local reruns.",
        "The governed machine outputs tracked for review remain enumerated below even",
        "when they live at the evidence root for compatibility with existing study builders.",
        "",
        "## Governed Primary Outputs",
        "",
    ]
    for path in primary_output_paths:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "## Local Execution Products",
            "",
            "- `reference-contract.json`",
            "- `analysis-run.json`",
            "",
            "These files are written by the bundle-local `reference.R` and `analysis.py`",
            "programs when contributors rerun one evidence unit directly.",
            "",
        ]
    )
    return "\n".join(lines)


def build_bundle_result_artifacts(
    repo_root: Path, bundle_root: Path, contract: dict[str, Any] | None = None
) -> dict[str, str | dict[str, Any]]:
    if contract is None:
        from .bundle_contracts import build_bundle_artifact_contract

        contract = build_bundle_artifact_contract(repo_root, bundle_root)
    return {
        RESULTS_README_RELATIVE_PATH: render_results_readme(contract),
        RESULTS_MANIFEST_RELATIVE_PATH: render_results_manifest(contract),
    }


def write_result_artifact(path: Path, payload: str | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, dict):
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return
    path.write_text(payload, encoding="utf-8")
