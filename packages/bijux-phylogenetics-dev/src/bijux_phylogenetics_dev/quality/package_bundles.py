"""Audit package dependency policy and built publish artifacts."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tomllib
from typing import Any
from uuid import uuid4
import zipfile

from ..trusted_process import run_text
from .policies import PUBLICATION_READINESS_POLICY_PATH

DEFAULT_ARTIFACTS_ROOT = Path("artifacts/root/package-bundles")


@dataclass(frozen=True)
class BundleIssue:
    """Describe one package bundle policy or artifact violation."""

    code: str
    path: str
    message: str


@dataclass(frozen=True)
class PackageBundlePolicy:
    """Define the governed publish artifact policy for one package."""

    package_name: str
    package_dir: str
    wheel_module_root: str
    allowed_dependencies: tuple[str, ...]
    required_sdist_entries: tuple[str, ...]
    required_wheel_entries: tuple[str, ...]
    forbidden_archive_prefixes: tuple[str, ...]


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_str_tuple(values: object) -> tuple[str, ...]:
    if not isinstance(values, list):
        return ()
    return tuple(value for value in values if isinstance(value, str))


def load_publication_readiness_settings(repo_root: Path) -> dict[str, Any]:
    """Load the repository publication-readiness settings."""
    payload = _load_toml(repo_root / PUBLICATION_READINESS_POLICY_PATH)
    return _as_dict(
        _as_dict(_as_dict(payload.get("tool")).get("bijux_phylogenetics")).get(
            "publication_readiness"
        )
    )


def load_package_bundle_policies(repo_root: Path) -> dict[str, PackageBundlePolicy]:
    """Load publishable package bundle policies from the governed settings."""
    readiness = load_publication_readiness_settings(repo_root)
    package_policy = _as_dict(readiness.get("package_policy"))
    return _load_policy_entries(package_policy)


def load_target_package_bundle_policies(
    repo_root: Path,
) -> dict[str, PackageBundlePolicy]:
    """Load governed target-package bundle policies from the repository settings."""
    readiness = load_publication_readiness_settings(repo_root)
    package_policy = _as_dict(readiness.get("target_package_policy"))
    return _load_policy_entries(package_policy)


def _load_policy_entries(
    package_policy: dict[str, Any],
) -> dict[str, PackageBundlePolicy]:
    policies: dict[str, PackageBundlePolicy] = {}
    for package_name, entry in package_policy.items():
        values = _as_dict(entry)
        policies[package_name] = PackageBundlePolicy(
            package_name=package_name,
            package_dir=str(values["package_dir"]),
            wheel_module_root=str(values["wheel_module_root"]),
            allowed_dependencies=_as_str_tuple(values.get("allowed_dependencies")),
            required_sdist_entries=_as_str_tuple(values.get("required_sdist_entries")),
            required_wheel_entries=_as_str_tuple(values.get("required_wheel_entries")),
            forbidden_archive_prefixes=_as_str_tuple(
                values.get("forbidden_archive_prefixes")
            ),
        )
    return policies


def build_dependency_policy_report(repo_root: Path) -> dict[str, Any]:
    """Build the dependency-policy section of the package bundle report."""
    policies = load_package_bundle_policies(repo_root)
    issues: list[BundleIssue] = []
    packages: list[dict[str, Any]] = []
    for package_name, policy in sorted(policies.items()):
        pyproject = _load_toml(repo_root / policy.package_dir / "pyproject.toml")
        project = _as_dict(pyproject.get("project"))
        dependencies = tuple(
            dependency
            for dependency in project.get("dependencies", [])
            if isinstance(dependency, str)
        )
        if dependencies != policy.allowed_dependencies:
            issues.append(
                BundleIssue(
                    code="dependency-policy-drift",
                    path=f"{policy.package_dir}/pyproject.toml",
                    message=(
                        f"{package_name} dependencies {dependencies!r} do not match "
                        f"the governed policy {policy.allowed_dependencies!r}"
                    ),
                )
            )
        packages.append(
            {
                "package_name": package_name,
                "package_dir": policy.package_dir,
                "dependencies": list(dependencies),
                "allowed_dependencies": list(policy.allowed_dependencies),
            }
        )
    return {
        "schema_version": 1,
        "package_count": len(packages),
        "packages": packages,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def _archive_members(path: Path) -> set[str]:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            return set(archive.namelist())
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as archive:
            members = set()
            for member in archive.getmembers():
                if (
                    member.name
                    and not member.isdir()
                    and not member.issym()
                    and not member.islnk()
                ):
                    _, _, relative = member.name.partition("/")
                    members.add(relative)
            return members
    raise ValueError(f"unsupported archive type: {path}")


def _archive_link_issues(path: Path, package_name: str) -> list[BundleIssue]:
    if not path.name.endswith(".tar.gz"):
        return []
    issues: list[BundleIssue] = []
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.name or not (member.issym() or member.islnk()):
                continue
            link_kind = "symbolic link" if member.issym() else "hard link"
            issues.append(
                BundleIssue(
                    code="linked-archive-entry",
                    path=path.name,
                    message=(
                        f"{package_name} sdist contains {link_kind} "
                        f"{member.name} -> {member.linkname}"
                    ),
                )
            )
    return issues


def _find_single(dist_dir: Path, pattern: str) -> Path:
    matches = sorted(dist_dir.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one artifact for pattern {pattern} in {dist_dir}"
        )
    return matches[0]


def audit_package_bundle_directory(
    policy: PackageBundlePolicy,
    dist_dir: Path,
) -> dict[str, Any]:
    """Audit one built wheel and sdist directory against the governed policy."""
    issues: list[BundleIssue] = []
    wheel_path = _find_single(dist_dir, "*.whl")
    sdist_path = _find_single(dist_dir, "*.tar.gz")
    wheel_members = _archive_members(wheel_path)
    sdist_members = _archive_members(sdist_path)
    issues.extend(_archive_link_issues(sdist_path, policy.package_name))

    for entry in policy.required_wheel_entries:
        if entry not in wheel_members:
            issues.append(
                BundleIssue(
                    code="missing-wheel-entry",
                    path=wheel_path.name,
                    message=f"{policy.package_name} wheel is missing {entry}",
                )
            )
    for entry in policy.required_sdist_entries:
        if entry not in sdist_members:
            issues.append(
                BundleIssue(
                    code="missing-sdist-entry",
                    path=sdist_path.name,
                    message=f"{policy.package_name} sdist is missing {entry}",
                )
            )
    for prefix in policy.forbidden_archive_prefixes:
        for archive_name, members in (
            (wheel_path.name, wheel_members),
            (sdist_path.name, sdist_members),
        ):
            if any(member.startswith(prefix) for member in members):
                issues.append(
                    BundleIssue(
                        code="forbidden-archive-prefix",
                        path=archive_name,
                        message=f"{policy.package_name} artifact includes forbidden prefix {prefix}",
                    )
                )

    return {
        "package_name": policy.package_name,
        "wheel_path": wheel_path.as_posix(),
        "sdist_path": sdist_path.as_posix(),
        "wheel_member_count": len(wheel_members),
        "sdist_member_count": len(sdist_members),
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def build_package_bundle_report(
    repo_root: Path,
    *,
    artifacts_root: Path,
    build_artifacts: bool,
    publish_artifacts: bool = True,
) -> dict[str, Any]:
    """Build the full publish artifact bundle report for governed packages."""
    policies = load_package_bundle_policies(repo_root)
    target_policies = load_target_package_bundle_policies(repo_root)
    settings = load_publication_readiness_settings(repo_root)
    package_reports: list[dict[str, Any]] = []
    issues: list[BundleIssue] = []
    expected_publishable_packages = tuple(
        entry
        for entry in settings.get("expected_publishable_packages", [])
        if isinstance(entry, str)
    )
    if tuple(sorted(policies)) != tuple(sorted(expected_publishable_packages)):
        issues.append(
            BundleIssue(
                code="publishable-package-policy-drift",
                path=PUBLICATION_READINESS_POLICY_PATH.as_posix(),
                message="package bundle policy entries must match the governed publishable package set exactly",
            )
        )
    target_shape_packages = tuple(
        entry
        for entry in settings.get("target_shape_packages", [])
        if isinstance(entry, str)
    )
    target_policy_reports: list[dict[str, Any]] = []
    for package_name in target_shape_packages:
        target_policy = target_policies.get(package_name)
        package_dir = (
            target_policy.package_dir
            if target_policy is not None
            else f"packages/{package_name}"
        )
        if package_name not in policies and target_policy is None:
            issues.append(
                BundleIssue(
                    code="missing-target-package-policy",
                    path=PUBLICATION_READINESS_POLICY_PATH.as_posix(),
                    message=f"target repository shape package {package_name} is missing a target package bundle policy",
                )
            )
        target_policy_reports.append(
            {
                "package_name": package_name,
                "has_publishable_policy": package_name in policies,
                "has_target_policy": target_policy is not None,
                "package_dir": package_dir,
                "package_exists": (repo_root / package_dir).is_dir(),
            }
        )
    artifacts_root.mkdir(parents=True, exist_ok=True)
    staging_root = (
        artifacts_root / f".package-bundles-{uuid4().hex}" if build_artifacts else None
    )
    if staging_root is not None:
        shutil.rmtree(staging_root, ignore_errors=True)
        staging_root.mkdir(parents=True, exist_ok=True)

    try:
        for package_name, policy in sorted(policies.items()):
            dist_dir = artifacts_root / package_name
            audit_dir = dist_dir
            if build_artifacts:
                if staging_root is None:
                    raise ValueError("staging root is required when building artifacts")
                audit_dir = staging_root / package_name
                shutil.rmtree(audit_dir, ignore_errors=True)
                audit_dir.mkdir(parents=True, exist_ok=True)
                run_text(
                    [
                        sys.executable,
                        "-m",
                        "build",
                        "--wheel",
                        "--sdist",
                        "--outdir",
                        str(audit_dir),
                        str(repo_root / policy.package_dir),
                    ],
                    check=True,
                    capture_output=True,
                )
            report = audit_package_bundle_directory(policy, audit_dir)
            if build_artifacts and publish_artifacts:
                shutil.rmtree(dist_dir, ignore_errors=True)
                shutil.copytree(audit_dir, dist_dir)
                report["wheel_path"] = str(dist_dir / Path(report["wheel_path"]).name)
                report["sdist_path"] = str(dist_dir / Path(report["sdist_path"]).name)
            package_reports.append(report)
            for issue in report["issues"]:
                issues.append(BundleIssue(**issue))
    finally:
        if staging_root is not None:
            shutil.rmtree(staging_root, ignore_errors=True)

    return {
        "schema_version": 1,
        "package_count": len(package_reports),
        "packages": package_reports,
        "expected_publishable_packages": list(expected_publishable_packages),
        "policy_package_names": sorted(policies),
        "target_package_policies": target_policy_reports,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def check_package_bundles(
    repo_root: Path,
    *,
    artifacts_root: Path | None = None,
    json_out: Path | None = None,
) -> dict[str, Any]:
    """Raise when the package bundle report contains any issues."""
    payload = build_package_bundle_report(
        repo_root.resolve(),
        artifacts_root=(artifacts_root or repo_root / DEFAULT_ARTIFACTS_ROOT).resolve(),
        build_artifacts=True,
        publish_artifacts=False,
    )
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["issue_count"]:
        raise SystemExit("package bundle audit failed")
    return payload


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the package bundle audit."""
    parser = argparse.ArgumentParser(
        description="Audit package dependency policy and publishable build bundles."
    )
    parser.add_argument("command", choices=("dependencies", "report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--artifacts-root", default=str(DEFAULT_ARTIFACTS_ROOT))
    parser.add_argument("--json-out", default="")
    return parser.parse_args()


def main() -> int:
    """Run the package bundle CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    artifacts_root = Path(args.artifacts_root).resolve()
    json_out = Path(args.json_out).resolve() if args.json_out else None
    if args.command == "dependencies":
        payload = build_dependency_policy_report(repo_root)
    else:
        if args.command == "check":
            payload = check_package_bundles(
                repo_root,
                artifacts_root=artifacts_root,
                json_out=json_out,
            )
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        payload = build_package_bundle_report(
            repo_root,
            artifacts_root=artifacts_root,
            build_artifacts=True,
        )
    if json_out is not None:
        _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
