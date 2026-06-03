from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics_dev.quality.artifact_governance import (
    build_artifact_governance_report,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / "tox.ini",
        """
[tox]
toxworkdir = {tox_root}/artifacts/root/tox
""".strip()
        + "\n",
    )
    _write(
        repo_root / "makes" / "root.mk",
        """
list-evidence-studies:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-studies.json"
build-evidence-book:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json"
build-evidence-study:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-book-build.json"
validate-evidence-book:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-book-validation.json"
report-evidence-completeness:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json"
check-evidence-completeness:
\t@echo "$(ROOT_ARTIFACTS_DIR)/evidence-completeness.json"
sync-evidence-artifacts:
\t@echo sync
check-evidence-artifacts:
\t@echo check
report-evidence-governance:
\t@$(MAKE) report-artifact-governance
check-evidence-governance:
\t@$(MAKE) check-artifact-governance
report-artifact-governance:
\t@echo "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json"
check-artifact-governance:
\t@echo "$(ROOT_ARTIFACTS_DIR)/artifact-governance.json"
report-execution-surfaces:
\t@echo "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json"
check-execution-surfaces:
\t@echo "$(ROOT_ARTIFACTS_DIR)/execution-surfaces.json"
report-package-boundaries:
\t@echo "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json"
check-package-boundaries:
\t@echo "$(ROOT_ARTIFACTS_DIR)/package-boundaries.json"
report-package-bundles:
\t@echo report-bundles
check-package-bundles:
\t@echo check-bundles
report-publish-readiness:
\t@echo report-publish
check-publish-readiness:
\t@echo check-publish
report-release-readiness:
\t@$(MAKE) report-publish-readiness
check-release-readiness:
\t@$(MAKE) check-publish-readiness
""".strip()
        + "\n",
    )
    return repo_root


def test_build_artifact_governance_report_accepts_governed_repo(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = build_artifact_governance_report(repo_root)

    assert report["issue_count"] == 0


def test_build_artifact_governance_report_flags_toxworkdir_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "tox.ini", "[tox]\ntoxworkdir = .tox\n")

    report = build_artifact_governance_report(repo_root)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "tox-artifact-path-drift" in issue_codes


def test_repository_artifact_governance_report_is_clean() -> None:
    report = build_artifact_governance_report(REPO_ROOT)

    assert report["issue_count"] == 0
