"""Govern per-evidence input manifests for the repository evidence-book."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from bijux_phylogenetics.evidence.bundle_contracts import (
    REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
    RESULTS_DIRNAME,
)

INPUT_MANIFEST_FILENAME = "inputs.manifest.json"
SKIP_BUNDLE_FILENAMES = {
    "manifest.json",
    "README.md",
    "reviewer-summary.json",
    "reviewer-summary.md",
    "claims.json",
    "claim_verdicts.json",
    INPUT_MANIFEST_FILENAME,
    *REQUIRED_BUNDLE_LOCAL_ARTIFACTS,
}
DATA_LIKE_EXTENSIONS = {
    ".csv",
    ".tsv",
    ".nwk",
    ".nex",
    ".json",
    ".txt",
    ".xlsx",
    ".rdata",
}
LOCAL_INPUT_EXTENSIONS = {".csv", ".tsv", ".nwk", ".nex", ".xlsx", ".rdata"}
EXPLICIT_LOCAL_INPUT_FILENAMES = {"expected-failure-cases.json"}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_inputs(
    manifest: dict[str, Any], bundle_root: Path, repo_root: Path
) -> list[dict[str, Any]]:
    source_inputs: list[dict[str, Any]] = []
    source_basis = manifest.get("source_basis", [])
    if not isinstance(source_basis, list):
        return source_inputs
    for index, entry in enumerate(source_basis, start=1):
        if not isinstance(entry, dict):
            continue
        locator = entry.get("locator")
        if not isinstance(locator, str) or not locator:
            continue
        source_payload: dict[str, Any] = {
            "input_id": f"source-input-{index:03d}",
            "kind": entry.get("kind", "unknown"),
            "label": entry.get("label", locator),
            "locator": locator,
            "read_only": not locator.startswith(
                bundle_root.relative_to(repo_root).as_posix()
            ),
        }
        if not locator.startswith("external:"):
            path = repo_root / locator
            if path.is_file():
                source_payload["sha256"] = _sha256(path)
            elif path.exists():
                source_payload["entry_count"] = sum(
                    1 for child in path.rglob("*") if child.is_file()
                )
        source_inputs.append(source_payload)
    return source_inputs


def _governed_local_artifacts(
    manifest: dict[str, Any], bundle_root: Path, repo_root: Path
) -> list[dict[str, Any]]:
    manifest_locators = {
        entry.get("locator")
        for entry in manifest.get("source_basis", [])
        if isinstance(entry, dict) and isinstance(entry.get("locator"), str)
    }
    local_inputs: list[dict[str, Any]] = []
    local_index = 0
    for child in sorted(bundle_root.iterdir()):
        if not child.is_file() or child.name in SKIP_BUNDLE_FILENAMES:
            continue
        relative_path = child.relative_to(repo_root).as_posix()
        is_manifest_declared_input = relative_path in manifest_locators
        is_data_like = child.suffix.lower() in DATA_LIKE_EXTENSIONS
        if not is_manifest_declared_input and not is_data_like:
            continue
        local_index += 1
        local_inputs.append(
            {
                "input_id": f"local-input-{local_index:03d}",
                "kind": (
                    "copied-reference-fragment"
                    if child.name.startswith("reference_")
                    else "governed-local-input"
                ),
                "label": child.name,
                "locator": relative_path,
                "read_only": False,
                "sha256": _sha256(child),
            }
        )
    results_root = bundle_root / RESULTS_DIRNAME
    if results_root.is_dir():
        for child in sorted(results_root.iterdir()):
            if not child.is_file():
                continue
            local_index += 1
            local_inputs.append(
                {
                    "input_id": f"local-input-{local_index:03d}",
                    "kind": "governed-result-surface",
                    "label": f"{RESULTS_DIRNAME}/{child.name}",
                    "locator": child.relative_to(repo_root).as_posix(),
                    "read_only": False,
                    "sha256": _sha256(child),
                }
            )
    return local_inputs


def _classify_local_input(entry: dict[str, Any]) -> str | None:
    label = str(entry.get("label", ""))
    suffix = Path(label).suffix.lower()
    if label.startswith("reference_"):
        return "copied-reference-fragment"
    if label.startswith("random_"):
        return "derived-local-input"
    if label in EXPLICIT_LOCAL_INPUT_FILENAMES:
        return "derived-local-input"
    if suffix in LOCAL_INPUT_EXTENSIONS:
        return "derived-local-input"
    return None


def _local_inputs(local_artifacts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    local_inputs: list[dict[str, Any]] = []
    for artifact in local_artifacts:
        input_class = _classify_local_input(artifact)
        if input_class is None:
            continue
        local_inputs.append(
            {
                **artifact,
                "input_class": input_class,
            }
        )
    return local_inputs


def build_inputs_manifest(repo_root: Path, bundle_root: Path) -> dict[str, Any]:
    manifest = _load_json(bundle_root / "manifest.json")
    source_inputs = _source_inputs(manifest, bundle_root, repo_root)
    local_artifacts = _governed_local_artifacts(manifest, bundle_root, repo_root)
    local_inputs = _local_inputs(local_artifacts)
    return {
        "schema_version": 2,
        "study_id": manifest["study_id"],
        "evidence_id": manifest["evidence_id"],
        "source_input_count": len(source_inputs),
        "governed_local_artifact_count": len(local_artifacts),
        "governed_local_artifacts": local_artifacts,
        "local_input_count": len(local_inputs),
        "source_inputs": source_inputs,
        "local_inputs": local_inputs,
    }


def iter_bundle_roots(repo_root: Path) -> list[Path]:
    studies_root = repo_root / "evidence-book" / "studies"
    return sorted(path for path in studies_root.glob("*/evidence-*") if path.is_dir())


def _selected_bundle_roots(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[Path]:
    roots = iter_bundle_roots(repo_root)
    if study_id is not None:
        roots = [root for root in roots if root.parent.name == study_id]
    if evidence_id is not None:
        roots = [root for root in roots if root.name == evidence_id]
    return roots


def sync_inputs_manifests(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[Path]:
    written: list[Path] = []
    for bundle_root in _selected_bundle_roots(
        repo_root,
        study_id=study_id,
        evidence_id=evidence_id,
    ):
        payload = build_inputs_manifest(repo_root, bundle_root)
        target = bundle_root / INPUT_MANIFEST_FILENAME
        _write_json(target, payload)
        written.append(target)
    return written


def check_inputs_manifests(
    repo_root: Path,
    *,
    study_id: str | None = None,
    evidence_id: str | None = None,
) -> list[str]:
    mismatches: list[str] = []
    for bundle_root in _selected_bundle_roots(
        repo_root,
        study_id=study_id,
        evidence_id=evidence_id,
    ):
        target = bundle_root / INPUT_MANIFEST_FILENAME
        expected = build_inputs_manifest(repo_root, bundle_root)
        if not target.is_file():
            mismatches.append(
                f"{target.relative_to(repo_root)}: missing {INPUT_MANIFEST_FILENAME}"
            )
            continue
        actual = _load_json(target)
        if actual != expected:
            mismatches.append(f"{target.relative_to(repo_root)}: stale input manifest")
    return mismatches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync or validate per-evidence input manifests for the evidence-book."
    )
    parser.add_argument("command", choices=("sync", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--study-id", default="")
    parser.add_argument("--evidence-id", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    study_id = args.study_id or None
    evidence_id = args.evidence_id or None
    if args.command == "sync":
        written = sync_inputs_manifests(
            repo_root,
            study_id=study_id,
            evidence_id=evidence_id,
        )
        print(
            json.dumps(
                {
                    "written": [
                        path.relative_to(repo_root).as_posix() for path in written
                    ]
                },
                indent=2,
            )
        )
        return 0
    mismatches = check_inputs_manifests(
        repo_root,
        study_id=study_id,
        evidence_id=evidence_id,
    )
    if mismatches:
        raise SystemExit(
            "evidence input manifest check failed:\n" + "\n".join(mismatches)
        )
    print(
        json.dumps(
            {
                "status": "ok",
                "bundle_count": len(
                    _selected_bundle_roots(
                        repo_root,
                        study_id=study_id,
                        evidence_id=evidence_id,
                    )
                ),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
