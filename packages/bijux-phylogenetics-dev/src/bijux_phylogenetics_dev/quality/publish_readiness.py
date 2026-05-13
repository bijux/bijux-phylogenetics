"""Build publish-readiness reports from repository-owned package and evidence metadata."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import shutil
import tomllib
from typing import Any

from bijux_phylogenetics.evidence.study_contracts import load_study_contract

from ..trusted_process import run_text
from .config_ssot import build_config_ssot_report
from .evidence_inputs import INPUT_MANIFEST_FILENAME, check_inputs_manifests
from .execution_surfaces import build_execution_surfaces_report
from .package_boundaries import build_package_boundary_report
from .package_bundles import (
    build_dependency_policy_report,
    load_package_bundle_policies,
    load_publication_readiness_settings,
    load_target_package_bundle_policies,
)
from .policies import PUBLICATION_READINESS_POLICY_PATH

TomlTable = dict[str, Any]
JsonObject = dict[str, Any]
GOVERNED_JUNK_NAMES = {".DS_Store", "__pycache__"}


@dataclass(frozen=True)
class ReadinessIssue:
    """Describe one publish-readiness problem."""

    code: str
    path: str
    message: str


@dataclass(frozen=True)
class Scorecard:
    """Summarize one publish-readiness dimension."""

    name: str
    status: str
    score: int
    max_score: int
    summary: str
    issue_count: int


def _load_json(path: Path) -> JsonObject:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _load_toml(path: Path) -> TomlTable:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _as_dict(value: object) -> JsonObject:
    return value if isinstance(value, dict) else {}


def _as_list(value: object) -> list[object]:
    return value if isinstance(value, list) else []


def _is_repo_relative_locator(locator: str) -> bool:
    return not locator.startswith("external:")


def _source_intake_policy_matches(intake_policy: str, sources: list[object]) -> bool:
    locators = [
        str(locator)
        for source in sources
        if isinstance(source, dict)
        for locator in [source.get("locator")]
        if isinstance(locator, str) and locator
    ]
    if intake_policy == "repository-owned-source":
        return bool(locators) and all(
            _is_repo_relative_locator(locator) for locator in locators
        )
    if intake_policy == "read-only-external-source":
        return any(not _is_repo_relative_locator(locator) for locator in locators)
    return False


def _output_classification(path: Path) -> str:
    if path.suffix == ".md":
        return "human-report"
    if path.suffix in {".json", ".csv", ".tsv", ".nwk", ".txt"}:
        return "machine-or-reference-output"
    return "other-governed-output"


def _as_str_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [entry for entry in value if isinstance(entry, str)]


def _issue_metadata(code: str) -> dict[str, str]:
    if code.startswith(("config-", "missing-root-make-target")):
        return {
            "domain": "standards",
            "owner_surface": "repository-governance",
            "severity": "error",
        }
    if code in {
        "runtime-example-fixture-coupling",
        "runtime-example-resource-missing",
        "runtime-evidence-helper-coupling",
    }:
        return {
            "domain": "runtime-boundary",
            "owner_surface": "bijux-phylogenetics",
            "severity": "error",
        }
    if code in {
        "missing-target-shape-package",
        "missing-target-package-policy",
        "runtime-owns-forbidden-subpackage",
    }:
        return {
            "domain": "repository-shape",
            "owner_surface": "package-layout",
            "severity": "error",
        }
    if code.startswith(("missing-evidence", "stale-evidence")):
        return {
            "domain": "evidence-program",
            "owner_surface": "evidence-book",
            "severity": "error",
        }
    if code.startswith("dependency-") or "package" in code or "runtime" in code:
        return {
            "domain": "package-policy",
            "owner_surface": "packaging",
            "severity": "error",
        }
    return {
        "domain": "governance",
        "owner_surface": "repository",
        "severity": "error",
    }


def scan_runtime_owner_drift(
    repo_root: Path, *, publication_settings: JsonObject
) -> list[ReadinessIssue]:
    """Report runtime surfaces that still depend on maintainer-owned material."""
    issues: list[ReadinessIssue] = []
    target_shape = _as_dict(publication_settings.get("target_repository_shape"))
    runtime_example_root = str(target_shape.get("runtime_example_input_root", ""))
    if runtime_example_root and not (repo_root / runtime_example_root).is_dir():
        issues.append(
            ReadinessIssue(
                code="runtime-example-resource-missing",
                path=runtime_example_root,
                message="runtime example workflows must ship packaged example inputs from the governed repository path",
            )
        )
    demo_path = (
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "core"
        / "demo.py"
    )
    if demo_path.is_file():
        demo_text = demo_path.read_text(encoding="utf-8")
        if "tests/fixtures" in demo_text:
            issues.append(
                ReadinessIssue(
                    code="runtime-example-fixture-coupling",
                    path=demo_path.relative_to(repo_root).as_posix(),
                    message="runtime example workflow still reads maintainer test fixtures instead of packaged example resources",
                )
            )
    for relative_path in (
        "packages/bijux-phylogenetics/src/bijux_phylogenetics/core/demo.py",
        "packages/bijux-phylogenetics/src/bijux_phylogenetics/bayesian/evidence.py",
        "packages/bijux-phylogenetics/src/bijux_phylogenetics/engines/evidence.py",
    ):
        path = repo_root / relative_path
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        if "bijux_phylogenetics.evidence.bundles" in text:
            issues.append(
                ReadinessIssue(
                    code="runtime-evidence-helper-coupling",
                    path=relative_path,
                    message="runtime-owned workflow still imports evidence bundle helpers instead of the governed provenance bundle contract",
                )
            )
    return issues


def iter_study_roots(repo_root: Path) -> list[Path]:
    """Return every governed evidence study root."""
    studies_root = repo_root / "evidence-book" / "studies"
    return sorted(path for path in studies_root.iterdir() if path.is_dir())


def _git_available(repo_root: Path) -> bool:
    return (repo_root / ".git").exists()


def _is_ignored_junk(repo_root: Path, path: Path) -> bool:
    if not _git_available(repo_root):
        return False
    git_bin = shutil.which("git")
    if git_bin is None:
        return False
    completed = run_text(
        [
            git_bin,
            "-C",
            str(repo_root),
            "check-ignore",
            "-q",
            str(path.relative_to(repo_root)),
        ],
        capture_output=True,
        check=False,
    )
    return completed.returncode == 0


def _study_payload(study_root: Path) -> JsonObject:
    return load_study_contract(study_root)


def _package_pyprojects(repo_root: Path) -> dict[str, JsonObject]:
    packages: dict[str, JsonObject] = {}
    for path in sorted((repo_root / "packages").glob("*/pyproject.toml")):
        payload = _load_toml(path)
        project = payload.get("project", {})
        if isinstance(project, dict):
            name = project.get("name")
            if isinstance(name, str):
                packages[name] = payload
    return packages


def scan_governed_junk(repo_root: Path) -> list[ReadinessIssue]:
    """Find junk files inside governed repository surfaces."""
    issues: list[ReadinessIssue] = []
    for root_name in ("evidence-book", "packages", "docs", "configs", "makes"):
        root = repo_root / root_name
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if _is_ignored_junk(repo_root, path):
                continue
            if path.name == ".DS_Store":
                issues.append(
                    ReadinessIssue(
                        code="governed-junk",
                        path=path.relative_to(repo_root).as_posix(),
                        message="governed repository surfaces must not contain .DS_Store files",
                    )
                )
            elif path.is_dir() and path.name == "__pycache__":
                issues.append(
                    ReadinessIssue(
                        code="governed-junk",
                        path=path.relative_to(repo_root).as_posix(),
                        message="governed repository surfaces must not contain __pycache__ directories",
                    )
                )
    return issues


def build_evidence_reproducibility_inventory(
    repo_root: Path,
    *,
    publication_settings: JsonObject | None = None,
) -> JsonObject:
    """Summarize study metadata, dataset registries, and checksum coverage."""
    publication_settings = publication_settings or {}
    required_input_manifest = str(
        publication_settings.get(
            "required_evidence_input_manifest",
            INPUT_MANIFEST_FILENAME,
        )
    )
    required_bundle_artifacts = tuple(
        _as_str_list(publication_settings.get("required_evidence_bundle_artifacts"))
    )
    studies: list[JsonObject] = []
    issues: list[ReadinessIssue] = []
    repo_dataset_checksum_count = 0
    evidence_manifest_count = 0
    evidence_input_manifest_count = 0
    evidence_output_checksum_count = 0
    missing_bundle_artifact_count = 0
    missing_input_manifest_count = 0
    output_classification_counts = {
        "human-report": 0,
        "machine-or-reference-output": 0,
        "other-governed-output": 0,
    }
    stale_input_manifest_issues = check_inputs_manifests(repo_root)

    for study_root in iter_study_roots(repo_root):
        study = _study_payload(study_root)
        study_id = str(study["study_id"])
        provenance_locator = study.get("provenance_descriptor_locator", "")
        dataset_registry_locator = study.get("dataset_registry_locator", "")
        study_scope = _as_dict(study.get("study_scope"))
        source_intake_policy = study.get("source_intake_policy", "")

        if not isinstance(provenance_locator, str) or not provenance_locator:
            issues.append(
                ReadinessIssue(
                    code="missing-study-provenance",
                    path=f"evidence-book/studies/{study_id}/README.md",
                    message="study metadata must declare a provenance descriptor locator",
                )
            )
            provenance_payload: JsonObject = {}
        else:
            provenance_path = repo_root / provenance_locator
            if not provenance_path.is_file():
                issues.append(
                    ReadinessIssue(
                        code="missing-study-provenance-file",
                        path=provenance_locator,
                        message="declared provenance descriptor file does not exist",
                    )
                )
                provenance_payload = {}
            else:
                provenance_payload = _load_json(provenance_path)
                if provenance_payload.get("intake_policy") != source_intake_policy:
                    issues.append(
                        ReadinessIssue(
                            code="provenance-intake-policy-mismatch",
                            path=provenance_locator,
                            message="study provenance descriptor must repeat the study source_intake_policy exactly",
                        )
                    )
                sources = _as_list(provenance_payload.get("sources"))
                if not _source_intake_policy_matches(source_intake_policy, sources):
                    issues.append(
                        ReadinessIssue(
                            code="provenance-intake-policy-mismatch",
                            path=provenance_locator,
                            message="study provenance intake_policy must match the actual source locator ownership model",
                        )
                    )
                declared_source_count = provenance_payload.get("source_count")
                if declared_source_count != len(sources):
                    issues.append(
                        ReadinessIssue(
                            code="provenance-source-count-mismatch",
                            path=provenance_locator,
                            message="provenance descriptor source_count must match the number of sources",
                        )
                    )
                for entry in sources:
                    source = _as_dict(entry)
                    required_fields = {
                        "source_id",
                        "kind",
                        "label",
                        "locator",
                        "read_only",
                    }
                    missing_fields = sorted(
                        field for field in required_fields if field not in source
                    )
                    if missing_fields:
                        issues.append(
                            ReadinessIssue(
                                code="incomplete-provenance-source",
                                path=provenance_locator,
                                message="provenance source is missing required fields: "
                                + ", ".join(missing_fields),
                            )
                        )

        if (
            not isinstance(dataset_registry_locator, str)
            or not dataset_registry_locator
        ):
            issues.append(
                ReadinessIssue(
                    code="missing-dataset-registry",
                    path=f"evidence-book/studies/{study_id}/README.md",
                    message="study metadata must declare a dataset registry locator",
                )
            )
            dataset_registry: JsonObject = {}
        else:
            dataset_registry_path = repo_root / dataset_registry_locator
            if not dataset_registry_path.is_file():
                issues.append(
                    ReadinessIssue(
                        code="missing-dataset-registry-file",
                        path=dataset_registry_locator,
                        message="declared dataset registry file does not exist",
                    )
                )
                dataset_registry = {}
            else:
                dataset_registry = _load_json(dataset_registry_path)
                datasets = _as_list(dataset_registry.get("datasets"))
                declared_dataset_count = dataset_registry.get("dataset_count")
                if declared_dataset_count != len(datasets):
                    issues.append(
                        ReadinessIssue(
                            code="dataset-count-mismatch",
                            path=dataset_registry_locator,
                            message="dataset registry dataset_count must match the number of datasets",
                        )
                    )
                for entry in datasets:
                    dataset = _as_dict(entry)
                    required_fields = {
                        "dataset_id",
                        "kind",
                        "label",
                        "locator",
                        "schema_summary",
                    }
                    missing_fields = sorted(
                        field for field in required_fields if field not in dataset
                    )
                    if missing_fields:
                        issues.append(
                            ReadinessIssue(
                                code="incomplete-dataset-entry",
                                path=dataset_registry_locator,
                                message="dataset entry is missing required fields: "
                                + ", ".join(missing_fields),
                            )
                        )

        if not study_scope:
            issues.append(
                ReadinessIssue(
                    code="missing-study-scope",
                    path=f"evidence-book/studies/{study_id}/README.md",
                    message="study metadata must declare a structured study_scope section",
                )
            )
        else:
            coverage_focus = _as_list(study_scope.get("coverage_focus"))
            untouched_source_locators = _as_list(
                study_scope.get("untouched_source_locators")
            )
            if not coverage_focus:
                issues.append(
                    ReadinessIssue(
                        code="missing-study-coverage-focus",
                        path=f"evidence-book/studies/{study_id}/README.md",
                        message="study_scope must describe the study coverage focus",
                    )
                )
            if not untouched_source_locators:
                issues.append(
                    ReadinessIssue(
                        code="missing-untouched-source-locators",
                        path=f"evidence-book/studies/{study_id}/README.md",
                        message="study_scope must declare the untouched source locators used by the study",
                    )
                )
        if source_intake_policy not in {
            "read-only-external-source",
            "repository-owned-source",
        }:
            issues.append(
                ReadinessIssue(
                    code="invalid-source-intake-policy",
                    path=f"evidence-book/studies/{study_id}/README.md",
                    message="study metadata must declare a supported source_intake_policy",
                )
            )

        datasets = _as_list(dataset_registry.get("datasets"))
        dataset_entries: list[JsonObject] = []
        for dataset_payload in datasets:
            entry = _as_dict(dataset_payload)
            locator = entry.get("locator")
            if isinstance(locator, str) and _is_repo_relative_locator(locator):
                target = repo_root / locator
                if target.exists() and target.is_file():
                    entry["sha256"] = _sha256(target)
                    repo_dataset_checksum_count += 1
                elif target.exists():
                    entry["entry_count"] = sum(
                        1 for _ in target.rglob("*") if _.is_file()
                    )
                else:
                    issues.append(
                        ReadinessIssue(
                            code="missing-dataset-locator",
                            path=locator,
                            message="repo-relative dataset locator does not exist",
                        )
                    )
            dataset_entries.append(entry)

        evidence_entries: list[JsonObject] = []
        for evidence_root in sorted(study_root.glob("evidence-*")):
            manifest_path = evidence_root / "manifest.json"
            if not manifest_path.is_file():
                issues.append(
                    ReadinessIssue(
                        code="missing-evidence-manifest",
                        path=evidence_root.relative_to(repo_root).as_posix(),
                        message="evidence directory is missing manifest.json",
                    )
                )
                continue
            evidence_manifest_count += 1
            manifest = _load_json(manifest_path)
            input_manifest_path = evidence_root / required_input_manifest
            if input_manifest_path.is_file():
                evidence_input_manifest_count += 1
                input_manifest = _load_json(input_manifest_path)
                if input_manifest.get("schema_version", 0) < 2:
                    issues.append(
                        ReadinessIssue(
                            code="stale-evidence-input-manifest",
                            path=input_manifest_path.relative_to(repo_root).as_posix(),
                            message="governed input manifest must use the schema that distinguishes local inputs from governed local output artifacts",
                        )
                    )
                if "governed_local_artifact_count" not in input_manifest:
                    issues.append(
                        ReadinessIssue(
                            code="missing-evidence-input-artifact-inventory",
                            path=input_manifest_path.relative_to(repo_root).as_posix(),
                            message="governed input manifest must enumerate the full local artifact inventory alongside true local inputs",
                        )
                    )
            else:
                missing_input_manifest_count += 1
                issues.append(
                    ReadinessIssue(
                        code="missing-evidence-input-manifest",
                        path=input_manifest_path.relative_to(repo_root).as_posix(),
                        message="evidence bundle is missing the governed input manifest companion file",
                    )
                )
            bundle_local_artifacts: list[JsonObject] = []
            for artifact_name in required_bundle_artifacts:
                artifact_path = evidence_root / artifact_name
                bundle_local_artifacts.append(
                    {
                        "path": artifact_path.relative_to(repo_root).as_posix(),
                        "present": artifact_path.is_file(),
                    }
                )
                if not artifact_path.is_file():
                    missing_bundle_artifact_count += 1
                    issues.append(
                        ReadinessIssue(
                            code="missing-evidence-bundle-artifact",
                            path=artifact_path.relative_to(repo_root).as_posix(),
                            message="evidence bundle is missing one of the required governed local artifacts",
                        )
                    )
            local_output_checksums: dict[str, str] = {}
            local_output_classification_counts = {
                "human-report": 0,
                "machine-or-reference-output": 0,
                "other-governed-output": 0,
            }
            for child in sorted(evidence_root.iterdir()):
                if child.name in {"README.md", "manifest.json"}:
                    continue
                if child.is_file():
                    classification = _output_classification(child)
                    local_output_checksums[child.relative_to(repo_root).as_posix()] = (
                        _sha256(child)
                    )
                    evidence_output_checksum_count += 1
                    output_classification_counts[classification] += 1
                    local_output_classification_counts[classification] += 1
            evidence_entries.append(
                {
                    "evidence_id": manifest.get("evidence_id", evidence_root.name),
                    "bundle_local_artifacts": bundle_local_artifacts,
                    "has_input_manifest": input_manifest_path.is_file(),
                    "input_manifest_locator": input_manifest_path.relative_to(
                        repo_root
                    ).as_posix(),
                    "governed_local_artifact_count": (
                        int(input_manifest.get("governed_local_artifact_count", 0))
                        if input_manifest_path.is_file()
                        else 0
                    ),
                    "local_input_count": (
                        int(input_manifest.get("local_input_count", 0))
                        if input_manifest_path.is_file()
                        else 0
                    ),
                    "output_classification_counts": local_output_classification_counts,
                    "output_checksum_count": len(local_output_checksums),
                    "output_checksums": local_output_checksums,
                }
            )

        studies.append(
            {
                "study_id": study_id,
                "source_intake_policy": source_intake_policy,
                "provenance_descriptor_locator": provenance_locator,
                "dataset_registry_locator": dataset_registry_locator,
                "study_scope": study_scope,
                "provenance": provenance_payload,
                "datasets": dataset_entries,
                "evidences": evidence_entries,
            }
        )

    junk_issues = scan_governed_junk(repo_root)
    issues.extend(junk_issues)
    for mismatch in stale_input_manifest_issues:
        path_text, _, detail = mismatch.partition(": ")
        issues.append(
            ReadinessIssue(
                code="stale-evidence-input-manifest",
                path=path_text,
                message=detail or "governed input manifest is stale",
            )
        )

    return {
        "evidence_manifest_count": evidence_manifest_count,
        "evidence_input_manifest_count": evidence_input_manifest_count,
        "evidence_output_checksum_count": evidence_output_checksum_count,
        "governed_junk_issue_count": len(junk_issues),
        "issues": [asdict(issue) for issue in issues],
        "missing_bundle_artifact_count": missing_bundle_artifact_count,
        "missing_input_manifest_count": missing_input_manifest_count,
        "output_classification_counts": output_classification_counts,
        "repo_dataset_checksum_count": repo_dataset_checksum_count,
        "schema_version": 1,
        "stale_input_manifest_count": len(stale_input_manifest_issues),
        "study_count": len(studies),
        "studies": studies,
    }


def check_publish_readiness(
    repo_root: Path,
    *,
    json_out: Path | None = None,
) -> JsonObject:
    """Build the publish-readiness report and fail if the repository is not ready."""
    payload = build_publish_readiness_report(repo_root.resolve())
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["summary"]["overall_status"] != "ready":
        raise SystemExit(
            "publish-readiness report failed: repository still has unresolved readiness issues"
        )
    return payload


def build_publish_readiness_report(repo_root: Path) -> JsonObject:
    """Build the repository publish-readiness report and scorecards."""
    repo_root = repo_root.resolve()
    publication_settings = load_publication_readiness_settings(repo_root)
    package_payloads = _package_pyprojects(repo_root)
    publishable_bundle_policies = load_package_bundle_policies(repo_root)
    target_bundle_policies = load_target_package_bundle_policies(repo_root)
    dependency_policy = build_dependency_policy_report(repo_root)
    config_report = build_config_ssot_report(repo_root)
    execution_surfaces = build_execution_surfaces_report(repo_root)
    package_boundaries = build_package_boundary_report(repo_root)
    evidence_inventory = build_evidence_reproducibility_inventory(
        repo_root,
        publication_settings=publication_settings,
    )
    runtime_owner_issues = scan_runtime_owner_drift(
        repo_root,
        publication_settings=publication_settings,
    )

    package_issues: list[ReadinessIssue] = []
    standards_issues: list[ReadinessIssue] = []
    shape_issues: list[ReadinessIssue] = []
    runtime = _as_dict(package_payloads.get("bijux-phylogenetics"))
    alias = _as_dict(package_payloads.get("phylogenetic"))
    dev = _as_dict(package_payloads.get("bijux-phylogenetics-dev"))

    runtime_deps = set(_as_dict(runtime.get("project")).get("dependencies", []))
    alias_deps = set(_as_dict(alias.get("project")).get("dependencies", []))
    dev_deps = set(_as_dict(dev.get("project")).get("dependencies", []))

    if any(dep.startswith("phylogenetic") for dep in runtime_deps):
        package_issues.append(
            ReadinessIssue(
                code="runtime-depends-on-alias",
                path="packages/bijux-phylogenetics/pyproject.toml",
                message="runtime package must not depend on the compatibility alias package",
            )
        )
    if alias_deps != {"bijux-phylogenetics>=0.1.0,<1.0"}:
        package_issues.append(
            ReadinessIssue(
                code="alias-dependency-drift",
                path="packages/phylogenetic/pyproject.toml",
                message="compatibility alias package should depend directly on the canonical runtime package only",
            )
        )
    if any(dep.startswith("bijux-phylogenetics") for dep in dev_deps):
        package_issues.append(
            ReadinessIssue(
                code="dev-package-runtime-dependency",
                path="packages/bijux-phylogenetics-dev/pyproject.toml",
                message="maintainer package should not depend on runtime packages as install-time requirements",
            )
        )

    for package_name, payload in package_payloads.items():
        hatch = _as_dict(_as_dict(payload.get("tool")).get("hatch"))
        build = _as_dict(hatch.get("build"))
        targets = _as_dict(build.get("targets"))
        if "sdist" not in targets or "wheel" not in targets:
            package_issues.append(
                ReadinessIssue(
                    code="missing-build-targets",
                    path=f"packages/{package_name}/pyproject.toml",
                    message="publishable package must declare both sdist and wheel build targets",
                )
            )

    package_issues.extend(
        ReadinessIssue(**issue) for issue in dependency_policy["issues"]
    )
    package_issues.extend(
        ReadinessIssue(**issue) for issue in package_boundaries["issues"]
    )
    standards_issues.extend(
        ReadinessIssue(**issue) for issue in execution_surfaces["issues"]
    )

    expected_publishable_packages = sorted(
        _as_str_list(publication_settings.get("expected_publishable_packages"))
    )
    actual_package_names = sorted(package_payloads)
    if sorted(actual_package_names) != expected_publishable_packages:
        package_issues.append(
            ReadinessIssue(
                code="publishable-package-set-drift",
                path="packages",
                message="observed repository packages do not match the governed publishable package set",
            )
        )

    target_shape_packages = sorted(
        _as_str_list(publication_settings.get("target_shape_packages"))
    )
    missing_target_shape_packages = sorted(
        set(target_shape_packages).difference(actual_package_names)
    )
    for package_name in missing_target_shape_packages:
        shape_issues.append(
            ReadinessIssue(
                code="missing-target-shape-package",
                path=f"packages/{package_name}",
                message="target repository shape requires an owned package that does not exist yet",
            )
        )
    for package_name in target_shape_packages:
        if (
            package_name in publishable_bundle_policies
            or package_name in target_bundle_policies
        ):
            continue
        shape_issues.append(
            ReadinessIssue(
                code="missing-target-package-policy",
                path=PUBLICATION_READINESS_POLICY_PATH.as_posix(),
                message="target repository shape package must have either a publishable or target-only bundle policy",
            )
        )

    package_issues.extend(runtime_owner_issues)

    root_make_path = repo_root / "makes" / "root.mk"
    root_make_text = (
        root_make_path.read_text(encoding="utf-8") if root_make_path.is_file() else ""
    )
    for target in _as_str_list(publication_settings.get("required_root_make_targets")):
        if target not in root_make_text:
            standards_issues.append(
                ReadinessIssue(
                    code="missing-root-make-target",
                    path="makes/root.mk",
                    message=f"repository root make surface must declare {target}",
                )
            )

    runtime_src_root = (
        repo_root / "packages" / "bijux-phylogenetics" / "src" / "bijux_phylogenetics"
    )
    observed_runtime_subpackages: list[str] = []
    if runtime_src_root.is_dir():
        observed_runtime_subpackages = sorted(
            path.name for path in runtime_src_root.iterdir() if path.is_dir()
        )
        for package_name in _as_str_list(
            publication_settings.get("forbidden_runtime_subpackages")
        ):
            if (runtime_src_root / package_name).is_dir():
                shape_issues.append(
                    ReadinessIssue(
                        code="runtime-owns-forbidden-subpackage",
                        path=f"packages/bijux-phylogenetics/src/bijux_phylogenetics/{package_name}",
                        message="runtime package still owns a subpackage that belongs in a separate consumer layer",
                    )
                )

    evidence_issues = [
        issue
        for issue in evidence_inventory["issues"]
        if issue["code"] not in {"governed-junk"}
    ]
    junk_issues = [
        issue
        for issue in evidence_inventory["issues"]
        if issue["code"] == "governed-junk"
    ]

    package_scorecard = Scorecard(
        name="package-boundaries",
        status="ready" if not package_issues and not shape_issues else "blocked",
        score=max(0, 4 - len(package_issues) - len(shape_issues)),
        max_score=4,
        summary="Package metadata, publishable package boundaries, and runtime ownership must align with the governed repository shape.",
        issue_count=len(package_issues),
    )
    evidence_scorecard = Scorecard(
        name="evidence-program",
        status="ready" if not evidence_issues else "blocked",
        score=max(0, 5 - len(evidence_issues)),
        max_score=5,
        summary="Evidence bundles must expose provenance, governed inputs, and side-by-side code surfaces before publication claims are allowed.",
        issue_count=len(evidence_issues),
    )
    config_scorecard = Scorecard(
        name="config-and-standards",
        status="ready"
        if config_report.issue_count == 0 and not standards_issues
        else "blocked",
        score=2 if config_report.issue_count == 0 and not standards_issues else 0,
        max_score=2,
        summary="Repository-owned config and make surfaces must satisfy the governed standards contract.",
        issue_count=config_report.issue_count + len(standards_issues),
    )
    provenance_scorecard = Scorecard(
        name="reproducibility-and-provenance",
        status="ready"
        if not junk_issues
        and evidence_inventory["repo_dataset_checksum_count"]
        and evidence_inventory["evidence_input_manifest_count"]
        and not evidence_inventory["stale_input_manifest_count"]
        else "needs-work",
        score=3
        if not junk_issues and evidence_inventory["repo_dataset_checksum_count"]
        else 1,
        max_score=3,
        summary="Study datasets, governed local inputs, and evidence outputs must remain traceable and free from junk or stale manifest drift.",
        issue_count=len(junk_issues)
        + int(evidence_inventory["stale_input_manifest_count"])
        + int(evidence_inventory["missing_input_manifest_count"]),
    )
    blocker_register = []
    for issue in (
        package_issues
        + standards_issues
        + shape_issues
        + [ReadinessIssue(**issue) for issue in evidence_issues]
    ):
        blocker_register.append({**asdict(issue), **_issue_metadata(issue.code)})
    runtime_closure_blockers = [
        issue["code"]
        for issue in blocker_register
        if issue["code"]
        in {
            "runtime-depends-on-alias",
            "alias-dependency-drift",
            "dev-package-runtime-dependency",
            "dependency-role-drift",
            "forbidden-runtime-top-level-export",
            "alias-local-surface-drift",
            "repo-import-boundary-drift",
            "invalid-compatibility-locator",
            "unsupported-compatibility-module",
            "missing-compatibility-module",
            "missing-compatibility-export",
            "missing-build-targets",
            "dependency-policy-drift",
            "publishable-package-set-drift",
            "runtime-example-fixture-coupling",
            "runtime-example-resource-missing",
            "runtime-evidence-helper-coupling",
            "runtime-owns-forbidden-subpackage",
        }
    ]
    evidence_closure_blockers = [
        issue["code"]
        for issue in blocker_register
        if issue["code"]
        in {
            "missing-study-provenance",
            "missing-study-provenance-file",
            "missing-dataset-registry",
            "missing-dataset-registry-file",
            "missing-evidence-manifest",
            "missing-evidence-input-manifest",
            "missing-evidence-input-artifact-inventory",
            "stale-evidence-input-manifest",
            "missing-evidence-bundle-artifact",
            "missing-target-package-policy",
        }
    ]
    standards_closure_blockers = [
        issue["code"]
        for issue in blocker_register
        if issue["code"]
        in {
            "missing-root-make-target",
            "missing-governed-root-target",
            "missing-governed-tox-env",
            "tox-command-drift",
        }
        or issue["code"].startswith("config-")
    ]
    closure_criteria = {
        "runtime_publishable": {
            "status": "ready" if not runtime_closure_blockers else "blocked",
            "pass_when": [
                "canonical runtime, alias, and maintainer package boundaries match governed dependency policy",
                "runtime top-level exports do not leak evidence-only helper surfaces",
                "compatibility alias package remains a thin alias surface with no local scientific ownership",
                "runtime-to-evidence compatibility is declared against supported public runtime locators",
                "runtime-owned modules remain limited to public comparative APIs, provenance contracts, and packaged example resources",
                "runtime example workflows consume packaged example resources instead of maintainer test fixtures",
                "runtime-owned workflow bundling uses the governed provenance bundle contract instead of evidence-only helper modules",
                "all governed publishable packages declare wheel and sdist targets",
            ],
            "blocker_codes": runtime_closure_blockers,
        },
        "evidence_program_publishable": {
            "status": "ready" if not evidence_closure_blockers else "blocked",
            "pass_when": [
                "every evidence bundle has a governed input manifest",
                "every evidence input manifest distinguishes actual local inputs from governed local output artifacts",
                "every evidence bundle carries side-by-side reference and python code surfaces",
                "study provenance and dataset registries resolve cleanly",
                "the evidence-book stays publishable as a governed repository surface with human-readable study and evidence ownership",
            ],
            "blocker_codes": evidence_closure_blockers,
        },
        "repository_standards_aligned": {
            "status": "ready"
            if not standards_closure_blockers and config_report.issue_count == 0
            else "blocked",
            "pass_when": [
                "repository-owned config SSOT audit is clean",
                "repository make, tox, and workflow execution surfaces are governed explicitly and stay separated by owned responsibilities",
                "root make surface exposes the governed publish-readiness commands",
                "the target repository shape is encoded explicitly in repository-owned policy instead of being implied by backlog prose",
            ],
            "blocker_codes": standards_closure_blockers,
        },
    }
    release_gate_status = "ready" if not blocker_register else "blocked"
    release_gate: JsonObject = {
        "status": release_gate_status,
        "publish_allowed": not blocker_register,
        "superficial_completion_refused": bool(blocker_register),
        "blocked_by": [issue["code"] for issue in blocker_register],
        "summary": (
            "publication claims are refused until package boundaries, evidence decomposition, and standards closure are all real"
            if blocker_register
            else "publication claims are supported by governed repository evidence"
        ),
    }
    publication_scorecard = Scorecard(
        name="publication-closure",
        status=release_gate_status,
        score=1 if release_gate_status == "ready" else 0,
        max_score=1,
        summary="The repository must refuse superficial publication claims until the governed closure criteria are met.",
        issue_count=len(blocker_register),
    )
    overall_status = (
        "blocked"
        if blocker_register
        else "ready"
        if all(
            score.status == "ready"
            for score in (
                package_scorecard,
                evidence_scorecard,
                config_scorecard,
                provenance_scorecard,
                publication_scorecard,
            )
        )
        else "needs-work"
    )

    return {
        "package_count": len(package_payloads),
        "package_names": sorted(package_payloads),
        "publication_settings": publication_settings,
        "target_repository_shape": _as_dict(
            publication_settings.get("target_repository_shape")
        ),
        "package_bundle_policy": {
            "publishable_policy_names": sorted(publishable_bundle_policies),
            "target_policy_names": sorted(target_bundle_policies),
        },
        "package_issues": [asdict(issue) for issue in package_issues],
        "dependency_policy": dependency_policy,
        "package_boundaries": package_boundaries,
        "execution_surfaces": execution_surfaces,
        "config_ssot": config_report.to_dict(),
        "evidence_inventory": evidence_inventory,
        "repository_shape": {
            "actual_package_names": actual_package_names,
            "expected_publishable_packages": expected_publishable_packages,
            "target_shape_packages": target_shape_packages,
            "missing_target_shape_packages": missing_target_shape_packages,
            "observed_runtime_subpackages": observed_runtime_subpackages,
        },
        "closure_criteria": closure_criteria,
        "blocker_register": {
            "issue_count": len(blocker_register),
            "issues": blocker_register,
        },
        "release_gate": release_gate,
        "scorecards": {
            "package_boundaries": asdict(package_scorecard),
            "evidence_program": asdict(evidence_scorecard),
            "config_and_standards": asdict(config_scorecard),
            "reproducibility_and_provenance": asdict(provenance_scorecard),
            "publication_closure": asdict(publication_scorecard),
        },
        "summary": {
            "overall_status": overall_status,
            "blocker_count": len(blocker_register),
            "study_count": evidence_inventory["study_count"],
            "evidence_manifest_count": evidence_inventory["evidence_manifest_count"],
            "evidence_input_manifest_count": evidence_inventory[
                "evidence_input_manifest_count"
            ],
        },
    }


def _write_json(path: Path, payload: JsonObject) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for publish-readiness reporting."""
    parser = argparse.ArgumentParser(
        description="Build repository publish-readiness reports from evidence and package metadata."
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root to inspect. Defaults to the current directory.",
    )
    parser.add_argument(
        "--json-out",
        default="",
        help="Optional path to write the report as JSON.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit non-zero when the repository is not publish-ready.",
    )
    return parser.parse_args()


def main() -> int:
    """Run the publish-readiness report CLI."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    json_out = Path(args.json_out).resolve() if args.json_out else None
    if args.check:
        payload = check_publish_readiness(repo_root, json_out=json_out)
    else:
        payload = build_publish_readiness_report(repo_root)
        if json_out is not None:
            _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
