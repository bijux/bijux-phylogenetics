"""Audit repository-owned config SSOT rules for bijux-phylogenetics."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
import tomllib
from typing import Any

from .policies import CONFIG_SSOT_POLICY_PATH

TomlTable = dict[str, Any]
LOCAL_MYPY_REFERENCE = re.compile(
    r"(?P<ref>\$\(PROJECT_DIR\)/mypy\.ini|packages/[^/\s\"']+/mypy\.ini)"
)


@dataclass(frozen=True)
class ConfigSsotIssue:
    """Describe one repository config SSOT violation."""

    code: str
    path: str
    message: str


@dataclass(frozen=True)
class ConfigSsotReport:
    """Structured repository config SSOT audit result."""

    repo_root: str
    policy_path: str
    issue_count: int
    issues: list[ConfigSsotIssue]
    required_root_files: list[str]
    package_config_allowlist: list[str]
    package_local_config_files: list[str]
    audited_paths: list[str]

    def to_dict(self) -> dict[str, object]:
        """Render the report as a JSON-serializable mapping."""
        payload = asdict(self)
        payload["issues"] = [asdict(issue) for issue in self.issues]
        return payload


@dataclass(frozen=True)
class ConfigSsotPolicy:
    """Repository policy for config source-of-truth enforcement."""

    required_root_files: tuple[str, ...]
    forbidden_package_config_filenames: tuple[str, ...]
    allowed_package_config_paths: tuple[str, ...]
    audit_paths: tuple[str, ...]
    expected_root_config_dir: str
    expected_root_make_dir: str
    expected_mypy_config_path: str


def _as_table(value: object) -> TomlTable:
    return value if isinstance(value, dict) else {}


def _load_toml(path: Path) -> TomlTable:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def load_config_ssot_policy(repo_root: Path) -> ConfigSsotPolicy:
    """Load the repository config SSOT policy."""
    policy_path = repo_root / CONFIG_SSOT_POLICY_PATH
    payload = _load_toml(policy_path)
    tool = _as_table(payload.get("tool"))
    workspace = _as_table(tool.get("bijux_phylogenetics"))
    config_ssot = _as_table(workspace.get("config_ssot"))

    def _tuple_from_list(name: str) -> tuple[str, ...]:
        values = config_ssot.get(name, [])
        if not isinstance(values, list) or not all(isinstance(item, str) for item in values):
            raise ValueError(f"{policy_path}: {name} must be a list of strings")
        return tuple(values)

    expected_root_config_dir = config_ssot.get("expected_root_config_dir")
    expected_root_make_dir = config_ssot.get("expected_root_make_dir")
    expected_mypy_config_path = config_ssot.get("expected_mypy_config_path")
    if not all(
        isinstance(value, str)
        for value in (
            expected_root_config_dir,
            expected_root_make_dir,
            expected_mypy_config_path,
        )
    ):
        raise ValueError(f"{policy_path}: expected root config fields must be strings")

    return ConfigSsotPolicy(
        required_root_files=_tuple_from_list("required_root_files"),
        forbidden_package_config_filenames=_tuple_from_list(
            "forbidden_package_config_filenames"
        ),
        allowed_package_config_paths=_tuple_from_list("allowed_package_config_paths"),
        audit_paths=_tuple_from_list("audit_paths"),
        expected_root_config_dir=expected_root_config_dir,
        expected_root_make_dir=expected_root_make_dir,
        expected_mypy_config_path=expected_mypy_config_path,
    )


def _workspace_metadata(repo_root: Path) -> TomlTable:
    payload = _load_toml(repo_root / "pyproject.toml")
    tool = _as_table(payload.get("tool"))
    return _as_table(tool.get("bijux_phylogenetics"))


def build_config_ssot_report(repo_root: Path) -> ConfigSsotReport:
    """Build a structured config SSOT audit report for the repository."""
    repo_root = repo_root.resolve()
    policy = load_config_ssot_policy(repo_root)
    workspace = _workspace_metadata(repo_root)
    issues: list[ConfigSsotIssue] = []

    for relative_path in policy.required_root_files:
        if not (repo_root / relative_path).is_file():
            issues.append(
                ConfigSsotIssue(
                    code="missing-root-config",
                    path=relative_path,
                    message="required repository-owned config file is missing",
                )
            )

    if workspace.get("config_dir") != policy.expected_root_config_dir:
        issues.append(
            ConfigSsotIssue(
                code="workspace-config-dir-drift",
                path="pyproject.toml",
                message=(
                    "tool.bijux_phylogenetics.config_dir must point at the repository "
                    f"config SSOT directory {policy.expected_root_config_dir!r}"
                ),
            )
        )
    if workspace.get("make_dir") != policy.expected_root_make_dir:
        issues.append(
            ConfigSsotIssue(
                code="workspace-make-dir-drift",
                path="pyproject.toml",
                message=(
                    "tool.bijux_phylogenetics.make_dir must point at the repository "
                    f"make SSOT directory {policy.expected_root_make_dir!r}"
                ),
            )
        )

    package_local_configs: list[str] = []
    allowlist = set(policy.allowed_package_config_paths)
    for package_root in sorted((repo_root / "packages").glob("*")):
        if not package_root.is_dir():
            continue
        for filename in policy.forbidden_package_config_filenames:
            candidate = package_root / filename
            if not candidate.exists():
                continue
            relative_path = candidate.relative_to(repo_root).as_posix()
            package_local_configs.append(relative_path)
            if relative_path not in allowlist:
                issues.append(
                    ConfigSsotIssue(
                        code="forbidden-package-config",
                        path=relative_path,
                        message=(
                            "package-local config drift is not allowed; move the rule to "
                            "the repository-owned configs/ surface or register an explicit "
                            "allowlist exception"
                        ),
                    )
                )

    for relative_path in policy.audit_paths:
        path = repo_root / relative_path
        if not path.is_file():
            issues.append(
                ConfigSsotIssue(
                    code="missing-audit-path",
                    path=relative_path,
                    message="configured config SSOT audit path does not exist",
                )
            )
            continue

        text = path.read_text(encoding="utf-8")
        for match in LOCAL_MYPY_REFERENCE.finditer(text):
            issues.append(
                ConfigSsotIssue(
                    code="local-mypy-reference",
                    path=relative_path,
                    message=(
                        "audit path still references package-local mypy config "
                        f"{match.group('ref')!r}; use {policy.expected_mypy_config_path!r}"
                    ),
                )
            )

        if "MYPY_CONFIG" in text and policy.expected_mypy_config_path not in text:
            issues.append(
                ConfigSsotIssue(
                    code="mypy-config-path-drift",
                    path=relative_path,
                    message=(
                        "audit path defines MYPY_CONFIG without pointing at the "
                        f"repository-owned config {policy.expected_mypy_config_path!r}"
                    ),
                )
            )
        if "QUALITY_MYPY_CONFIG" in text and policy.expected_mypy_config_path not in text:
            issues.append(
                ConfigSsotIssue(
                    code="quality-mypy-config-path-drift",
                    path=relative_path,
                    message=(
                        "audit path defines QUALITY_MYPY_CONFIG without pointing at the "
                        f"repository-owned config {policy.expected_mypy_config_path!r}"
                    ),
                )
            )

    return ConfigSsotReport(
        repo_root=repo_root.as_posix(),
        policy_path=CONFIG_SSOT_POLICY_PATH.as_posix(),
        issue_count=len(issues),
        issues=issues,
        required_root_files=list(policy.required_root_files),
        package_config_allowlist=list(policy.allowed_package_config_paths),
        package_local_config_files=sorted(package_local_configs),
        audited_paths=list(policy.audit_paths),
    )


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def check_config_ssot(repo_root: Path, *, json_out: Path | None = None) -> ConfigSsotReport:
    """Raise on config SSOT violations and optionally write the audit report."""
    report = build_config_ssot_report(repo_root)
    if json_out is not None:
        _write_json(json_out, report.to_dict())
    if report.issues:
        details = "\n".join(
            f"- [{issue.code}] {issue.path}: {issue.message}" for issue in report.issues
        )
        raise SystemExit(f"config SSOT audit failed with {report.issue_count} issue(s):\n{details}")
    return report


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for config SSOT audit commands."""
    parser = argparse.ArgumentParser(
        description="Audit repository-owned config source-of-truth rules."
    )
    parser.add_argument(
        "command",
        choices=("audit", "check"),
        help="Emit a JSON audit report or fail on any config SSOT drift.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root to audit. Defaults to the current directory.",
    )
    parser.add_argument(
        "--json-out",
        default="",
        help="Optional path to write the JSON audit report.",
    )
    return parser.parse_args()


def main() -> int:
    """Run the config SSOT audit CLI."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    json_out = Path(args.json_out).resolve() if args.json_out else None

    if args.command == "audit":
        report = build_config_ssot_report(repo_root)
        payload = report.to_dict()
        if json_out is not None:
            _write_json(json_out, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    check_config_ssot(repo_root, json_out=json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
