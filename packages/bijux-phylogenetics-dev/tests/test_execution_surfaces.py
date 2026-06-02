from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics_dev.quality.execution_surfaces import (
    build_execution_surfaces_report,
)
from bijux_phylogenetics_dev.quality.policies import EXECUTION_SURFACES_POLICY_PATH

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / EXECUTION_SURFACES_POLICY_PATH,
        """
[tool.bijux_phylogenetics.execution_surfaces]
required_root_make_targets = [
  "check-evidence-governance:",
  "report-evidence-governance:",
  "check-evidence-completeness:",
  "report-evidence-completeness:",
  "check-execution-surfaces:",
  "report-execution-surfaces:",
]
required_tox_envs = [
  "repository-contracts",
  "config-ssot",
  "evidence-governance",
  "evidence-completeness",
  "publish-readiness",
  "release-readiness-gate",
]

[tool.bijux_phylogenetics.execution_surfaces.tox_commands]
repository-contracts = ["make check-shared-bijux-py check-config-layout check-make-layout help"]
config-ssot = ["make check-config-ssot"]
evidence-governance = ["make check-evidence-governance"]
evidence-completeness = ["make check-evidence-completeness"]
publish-readiness = ["make report-release-readiness"]
release-readiness-gate = ["make check-release-readiness"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "makes" / "root.mk",
        """
check-evidence-governance:
\t@echo ok
report-evidence-governance:
\t@echo ok
check-evidence-completeness:
\t@echo ok
report-evidence-completeness:
\t@echo ok
check-execution-surfaces:
\t@echo ok
report-execution-surfaces:
\t@echo ok
""".strip()
        + "\n",
    )
    _write(
        repo_root / "tox.ini",
        """
[tox]
envlist =
    repository-contracts
    config-ssot
    evidence-governance
    evidence-completeness
    publish-readiness
    release-readiness-gate

[testenv:repository-contracts]
commands =
    make check-shared-bijux-py check-config-layout check-make-layout help

[testenv:config-ssot]
commands =
    make check-config-ssot

[testenv:evidence-governance]
commands =
    make check-evidence-governance

[testenv:evidence-completeness]
commands =
    make check-evidence-completeness

[testenv:publish-readiness]
commands =
    make report-release-readiness

[testenv:release-readiness-gate]
commands =
    make check-release-readiness
""".strip()
        + "\n",
    )
    return repo_root


def test_build_execution_surfaces_report_accepts_governed_repo(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = build_execution_surfaces_report(repo_root)

    assert report["issue_count"] == 0


def test_build_execution_surfaces_report_flags_missing_tox_env(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "tox.ini", "[tox]\nenvlist = repository-contracts\n")

    report = build_execution_surfaces_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "missing-governed-tox-env" in issue_codes


def test_repository_execution_surfaces_report_is_clean() -> None:
    report = build_execution_surfaces_report(REPO_ROOT)

    assert report["issue_count"] == 0
