"""Build publish-readiness reports from repository-owned package and evidence metadata."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import tomllib
from typing import Any

from .config_ssot import build_config_ssot_report

TomlTable = dict[str, Any]
JsonObject = dict[str, object]
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
    return json.loads(path.read_text(encoding="utf-8"))


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


def _output_classification(path: Path) -> str:
    if path.suffix == ".md":
        return "human-report"
    if path.suffix in {".json", ".csv", ".tsv", ".nwk", ".txt"}:
        return "machine-or-reference-output"
    return "other-governed-output"


def iter_study_roots(repo_root: Path) -> list[Path]:
    studies_root = repo_root / "evidence-book" / "studies"
    return sorted(path for path in studies_root.iterdir() if path.is_dir())


def _git_available(repo_root: Path) -> bool:
    return (repo_root / ".git").exists()


def _is_ignored_junk(repo_root: Path, path: Path) -> bool:
    if not _git_available(repo_root):
        return False
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(repo_root),
            "check-ignore",
            "-q",
            str(path.relative_to(repo_root)),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.returncode == 0


def _study_payload(study_root: Path) -> JsonObject:
    return _load_json(study_root / "study.json")


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


def build_evidence_reproducibility_inventory(repo_root: Path) -> JsonObject:
    """Summarize study metadata, dataset registries, and checksum coverage."""
    studies: list[JsonObject] = []
    issues: list[ReadinessIssue] = []
    repo_dataset_checksum_count = 0
    evidence_manifest_count = 0
    evidence_output_checksum_count = 0
    output_classification_counts = {
        "human-report": 0,
        "machine-or-reference-output": 0,
        "other-governed-output": 0,
    }

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
                    path=f"evidence-book/studies/{study_id}/study.json",
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
                    path=f"evidence-book/studies/{study_id}/study.json",
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
                    path=f"evidence-book/studies/{study_id}/study.json",
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
                        path=f"evidence-book/studies/{study_id}/study.json",
                        message="study_scope must describe the study coverage focus",
                    )
                )
            if not untouched_source_locators:
                issues.append(
                    ReadinessIssue(
                        code="missing-untouched-source-locators",
                        path=f"evidence-book/studies/{study_id}/study.json",
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
                    path=f"evidence-book/studies/{study_id}/study.json",
                    message="study metadata must declare a supported source_intake_policy",
                )
            )

        datasets = _as_list(dataset_registry.get("datasets"))
        dataset_entries: list[JsonObject] = []
        for dataset in datasets:
            entry = _as_dict(dataset)
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

    return {
        "evidence_manifest_count": evidence_manifest_count,
        "evidence_output_checksum_count": evidence_output_checksum_count,
        "governed_junk_issue_count": len(junk_issues),
        "issues": [asdict(issue) for issue in issues],
        "output_classification_counts": output_classification_counts,
        "repo_dataset_checksum_count": repo_dataset_checksum_count,
        "schema_version": 1,
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
    if payload["summary"]["overall_status"] != "ready":  # type: ignore[index]
        raise SystemExit(
            "publish-readiness report failed: repository still has unresolved readiness issues"
        )
    return payload


def build_publish_readiness_report(repo_root: Path) -> JsonObject:
    """Build the repository publish-readiness report and scorecards."""
    repo_root = repo_root.resolve()
    package_payloads = _package_pyprojects(repo_root)
    config_report = build_config_ssot_report(repo_root)
    evidence_inventory = build_evidence_reproducibility_inventory(repo_root)

    package_issues: list[ReadinessIssue] = []
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

    evidence_issues = [
        issue
        for issue in evidence_inventory["issues"]  # type: ignore[index]
        if issue["code"] not in {"governed-junk"}
    ]
    junk_issues = [
        issue
        for issue in evidence_inventory["issues"]  # type: ignore[index]
        if issue["code"] == "governed-junk"
    ]

    package_scorecard = Scorecard(
        name="package-boundaries",
        status="ready" if not package_issues else "blocked",
        score=max(0, 3 - len(package_issues)),
        max_score=3,
        summary="Package metadata keeps runtime, alias, and maintainer responsibilities separated.",
        issue_count=len(package_issues),
    )
    evidence_scorecard = Scorecard(
        name="evidence-program",
        status="ready" if not evidence_issues else "needs-work",
        score=max(0, 4 - len(evidence_issues)),
        max_score=4,
        summary="Evidence studies expose provenance descriptors, dataset registries, and evidence manifests for publish-time review.",
        issue_count=len(evidence_issues),
    )
    config_scorecard = Scorecard(
        name="config-and-standards",
        status="ready" if config_report.issue_count == 0 else "blocked",
        score=1 if config_report.issue_count == 0 else 0,
        max_score=1,
        summary="Repository-owned config SSOT audit is clean.",
        issue_count=config_report.issue_count,
    )
    provenance_scorecard = Scorecard(
        name="reproducibility-and-provenance",
        status="ready"
        if not junk_issues and evidence_inventory["repo_dataset_checksum_count"]  # type: ignore[index]
        else "needs-work",
        score=2 if not junk_issues else 1,
        max_score=2,
        summary="Study datasets and evidence outputs are traceable and governed surfaces are free from junk files.",
        issue_count=len(junk_issues),
    )

    return {
        "package_count": len(package_payloads),
        "package_names": sorted(package_payloads),
        "package_issues": [asdict(issue) for issue in package_issues],
        "config_ssot": config_report.to_dict(),
        "evidence_inventory": evidence_inventory,
        "scorecards": {
            "package_boundaries": asdict(package_scorecard),
            "evidence_program": asdict(evidence_scorecard),
            "config_and_standards": asdict(config_scorecard),
            "reproducibility_and_provenance": asdict(provenance_scorecard),
        },
        "summary": {
            "overall_status": "ready"
            if all(
                score.status == "ready"
                for score in (
                    package_scorecard,
                    evidence_scorecard,
                    config_scorecard,
                    provenance_scorecard,
                )
            )
            else "needs-work",
            "study_count": evidence_inventory["study_count"],
            "evidence_manifest_count": evidence_inventory["evidence_manifest_count"],
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
