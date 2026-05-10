from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics_dev.quality.evidence_inputs import build_inputs_manifest
from bijux_phylogenetics_dev.quality.publish_readiness import (
    build_publish_readiness_report,
    check_publish_readiness,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / "pyproject.toml",
        """
[tool.bijux_phylogenetics]
config_dir = "configs"
make_dir = "makes"
""".strip()
        + "\n",
    )
    _write(
        repo_root / "configs" / "config_ssot.toml",
        """
[tool.bijux_phylogenetics.config_ssot]
required_root_files = ["configs/mypy.ini", "configs/pytest.ini"]
forbidden_package_config_filenames = ["mypy.ini", "pytest.ini"]
allowed_package_config_paths = []
audit_paths = ["tox.ini", "makes/packages/runtime.mk"]
expected_root_config_dir = "configs"
expected_root_make_dir = "makes"
expected_mypy_config_path = "configs/mypy.ini"
""".strip()
        + "\n",
    )
    _write(
        repo_root / "configs" / "publication_readiness.toml",
        """
[tool.bijux_phylogenetics.publication_readiness]
required_evidence_input_manifest = "inputs.manifest.json"
required_evidence_bundle_code_files = ["reference.R", "analysis.py"]
expected_publishable_packages = [
    "bijux-phylogenetics",
    "phylogenetic",
    "bijux-phylogenetics-dev",
]
target_shape_packages = [
    "bijux-phylogenetics",
    "phylogenetic",
    "bijux-phylogenetics-dev",
]
forbidden_runtime_subpackages = ["evidence"]
required_root_make_targets = [
    "report-package-bundles:",
    "check-package-bundles:",
    "report-publish-readiness:",
    "check-publish-readiness:",
]
""".strip()
        + "\n",
    )
    _write(repo_root / "configs" / "mypy.ini", "[mypy]\nstrict = true\n")
    _write(repo_root / "configs" / "pytest.ini", "[pytest]\naddopts = -ra\n")
    _write(repo_root / "tox.ini", "[tox]\nminversion = 4.11\n")
    _write(
        repo_root / "makes" / "root.mk",
        "\n".join(
            [
                "report-package-bundles:",
                "check-package-bundles:",
                "report-publish-readiness:",
                "check-publish-readiness:",
            ]
        )
        + "\n",
    )
    _write(
        repo_root / "makes" / "packages" / "runtime.mk",
        "MYPY_CONFIG = $(MONOREPO_ROOT)/configs/mypy.ini\n",
    )
    _write(
        repo_root / "packages" / "bijux-phylogenetics" / "pyproject.toml",
        """
[project]
name = "bijux-phylogenetics"
dependencies = ["biopython>=1.0"]

[tool.hatch.build.targets.sdist]
include = ["src/bijux_phylogenetics/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/bijux_phylogenetics"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    _write(
        repo_root / "packages" / "phylogenetic" / "pyproject.toml",
        """
[project]
name = "phylogenetic"
dependencies = ["bijux-phylogenetics>=0.1.0,<1.0"]

[tool.hatch.build.targets.sdist]
include = ["src/phylogenetic/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/phylogenetic"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "phylogenetic" / "src" / "phylogenetic" / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    _write(
        repo_root / "packages" / "bijux-phylogenetics-dev" / "pyproject.toml",
        """
[project]
name = "bijux-phylogenetics-dev"
dependencies = ["PyYAML>=6.0"]

[tool.hatch.build.targets.sdist]
include = ["src/bijux_phylogenetics_dev/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/bijux_phylogenetics_dev"]
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics-dev"
        / "src"
        / "bijux_phylogenetics_dev"
        / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "study.json",
        """
{
  "dataset_registry_locator": "evidence-book/studies/demo-study/datasets/registry.json",
  "owner_package": "bijux-phylogenetics",
  "provenance_descriptor_locator": "evidence-book/studies/demo-study/provenance/sources.json",
  "source_intake_policy": "repository-owned-source",
  "study_categories": ["scientific-validation"],
  "study_id": "demo-study",
  "study_scope": {
    "coverage_focus": ["demo-analysis"],
    "untouched_source_locators": ["packages/bijux-phylogenetics/tests/fixtures/demo.tsv"]
  },
  "study_title": "Demo Study",
  "summary": "Minimal governed study."
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "provenance"
        / "sources.json",
        """
{
  "intake_policy": "repository-owned-source",
  "schema_version": 1,
  "source_count": 1,
  "sources": [
    {
      "source_id": "demo-fixture",
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "read_only": true
    }
  ],
  "study_id": "demo-study"
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "datasets"
        / "registry.json",
        """
{
  "dataset_count": 1,
  "datasets": [
    {
      "dataset_id": "dataset-001",
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "schema_summary": "Small governed fixture."
    }
  ],
  "schema_version": 1,
  "study_id": "demo-study"
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "manifest.json",
        """
{
  "evidence_id": "evidence-001",
  "study_id": "demo-study",
  "freshness": {
    "governed_code_paths": ["packages/bijux-phylogenetics/src/bijux_phylogenetics"],
    "last_generated_on": "2026-05-10",
    "source_basis_locators": [
      "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
      "evidence-book/studies/demo-study/evidence-001/parity.json"
    ]
  },
  "ownership": {
    "analytical_surfaces": ["demo-analysis"],
    "owner_package": "bijux-phylogenetics"
  },
  "source_basis": [
    {
      "kind": "repository-fixture",
      "label": "Demo fixture",
      "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv"
    },
    {
      "kind": "repository-reference",
      "label": "Local parity payload",
      "locator": "evidence-book/studies/demo-study/evidence-001/parity.json"
    }
  ]
}
""".strip()
        + "\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "analysis.py",
        "def run() -> dict[str, str]:\n    return {'status': 'ok'}\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "reference.R",
        "status <- 'ok'\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "parity.json",
        '{"status": "ok"}\n',
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "report.md",
        "# Demo report\n",
    )
    _write(
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "results.json",
        '{"status": "ok"}\n',
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "fixtures"
        / "demo.tsv",
        "species\tvalue\nA\t1\n",
    )
    bundle_root = repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001"
    inputs_manifest = build_inputs_manifest(repo_root, bundle_root)
    _write(
        bundle_root / "inputs.manifest.json",
        json.dumps(inputs_manifest, indent=2, sort_keys=True) + "\n",
    )
    return repo_root


def test_build_publish_readiness_report_exposes_repository_blockers() -> None:
    report = build_publish_readiness_report(REPO_ROOT)

    assert report["package_count"] == 3
    assert report["summary"]["overall_status"] == "blocked"
    assert report["summary"]["blocker_count"] > 0
    assert report["summary"]["study_count"] == 4
    assert report["summary"]["evidence_manifest_count"] == 19
    assert report["summary"]["evidence_input_manifest_count"] == 19
    assert report["config_ssot"]["issue_count"] == 0
    assert report["evidence_inventory"]["governed_junk_issue_count"] == 0
    assert report["evidence_inventory"]["repo_dataset_checksum_count"] >= 4
    assert report["evidence_inventory"]["evidence_output_checksum_count"] >= 10
    assert report["scorecards"]["package_boundaries"]["status"] == "blocked"
    assert report["scorecards"]["evidence_program"]["status"] == "blocked"
    blocker_codes = {issue["code"] for issue in report["blocker_register"]["issues"]}
    assert "missing-target-shape-package" in blocker_codes
    assert "runtime-owns-forbidden-subpackage" in blocker_codes
    assert "missing-evidence-bundle-code-file" in blocker_codes
    assert report["release_gate"]["publish_allowed"] is False
    assert report["release_gate"]["superficial_completion_refused"] is True


def test_build_publish_readiness_report_flags_governed_junk(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "docs" / ".DS_Store", "")

    report = build_publish_readiness_report(repo_root)

    assert report["summary"]["overall_status"] == "needs-work"
    assert report["evidence_inventory"]["governed_junk_issue_count"] == 1
    issue_codes = {issue["code"] for issue in report["evidence_inventory"]["issues"]}
    assert "governed-junk" in issue_codes
    assert (
        report["scorecards"]["reproducibility_and_provenance"]["status"] == "needs-work"
    )


def test_build_publish_readiness_report_flags_provenance_policy_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    provenance_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "provenance"
        / "sources.json"
    )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["intake_policy"] = "read-only-external-source"
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    report = build_publish_readiness_report(repo_root)

    issue_codes = {issue["code"] for issue in report["evidence_inventory"]["issues"]}
    assert "provenance-intake-policy-mismatch" in issue_codes
    assert report["scorecards"]["evidence_program"]["status"] == "blocked"


def test_build_publish_readiness_report_accepts_ready_minimal_repository(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = build_publish_readiness_report(repo_root)

    assert report["summary"]["overall_status"] == "ready"
    assert report["summary"]["blocker_count"] == 0
    assert report["scorecards"]["publication_closure"]["status"] == "ready"
    assert report["release_gate"]["publish_allowed"] is True


def test_check_publish_readiness_writes_json_and_fails_when_repo_is_not_ready(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(repo_root / "packages" / "bijux-phylogenetics" / "mypy.ini", "[mypy]\n")
    report_path = repo_root / "artifacts" / "root" / "publish-readiness.json"

    with pytest.raises(SystemExit, match="publish-readiness report failed"):
        check_publish_readiness(repo_root, json_out=report_path)

    assert report_path.is_file()
