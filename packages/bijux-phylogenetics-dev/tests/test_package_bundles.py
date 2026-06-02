from __future__ import annotations

from pathlib import Path
import tarfile
import zipfile

import pytest

from bijux_phylogenetics_dev.quality.package_bundles import (
    audit_package_bundle_directory,
    build_dependency_policy_report,
    build_package_bundle_report,
    check_package_bundles,
    load_package_bundle_policies,
    load_publication_readiness_settings,
)
from bijux_phylogenetics_dev.quality.policies import PUBLICATION_READINESS_POLICY_PATH

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_wheel(path: Path, members: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        for member, payload in members.items():
            archive.writestr(member, payload)


def _write_sdist(path: Path, members: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as archive:
        for member, payload in members.items():
            temp_path = path.parent / member.replace("/", "_")
            temp_path.write_text(payload, encoding="utf-8")
            archive.add(temp_path, arcname=f"demo-0.1.0/{member}")
            temp_path.unlink()


def _write_sdist_with_symbolic_link(
    path: Path,
    *,
    members: dict[str, str],
    linked_member: str,
    link_target: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as archive:
        for member, payload in members.items():
            temp_path = path.parent / member.replace("/", "_")
            temp_path.write_text(payload, encoding="utf-8")
            archive.add(temp_path, arcname=f"demo-0.1.0/{member}")
            temp_path.unlink()
        tarinfo = tarfile.TarInfo(name=f"demo-0.1.0/{linked_member}")
        tarinfo.type = tarfile.SYMTYPE
        tarinfo.linkname = link_target
        archive.addfile(tarinfo)


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / PUBLICATION_READINESS_POLICY_PATH,
        """
[tool.bijux_phylogenetics.publication_readiness]
required_evidence_input_manifest = "inputs.manifest.json"
required_evidence_bundle_code_files = ["reference.R", "analysis.py"]
required_evidence_bundle_artifacts = ["reference.R", "analysis.py", "checks.json", "report.md", "provenance.json"]
expected_publishable_packages = ["demo-runtime", "demo-dev"]
target_shape_packages = ["demo-runtime", "demo-dev", "demo-evidence"]
forbidden_runtime_subpackages = ["evidence"]
required_root_make_targets = ["validate-evidence-book:", "report-evidence-completeness:", "check-evidence-completeness:", "check-evidence-governance:", "sync-evidence-artifacts:", "check-evidence-artifacts:", "check-artifact-governance:", "check-package-bundles:", "check-release-readiness:"]

[tool.bijux_phylogenetics.publication_readiness.package_policy."demo-runtime"]
package_dir = "packages/demo-runtime"
wheel_module_root = "demo_runtime"
allowed_dependencies = ["numpy>=1.0"]
required_sdist_entries = ["README.md", "src/demo_runtime/__init__.py"]
required_wheel_entries = ["demo_runtime/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]

[tool.bijux_phylogenetics.publication_readiness.package_policy."demo-dev"]
package_dir = "packages/demo-dev"
wheel_module_root = "demo_dev"
allowed_dependencies = ["PyYAML>=6.0"]
required_sdist_entries = ["README.md", "src/demo_dev/__init__.py"]
required_wheel_entries = ["demo_dev/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]

[tool.bijux_phylogenetics.publication_readiness.target_package_policy."demo-evidence"]
package_dir = "packages/demo-evidence"
wheel_module_root = "demo_evidence"
allowed_dependencies = ["demo-runtime>=0.1.0,<1.0"]
required_sdist_entries = ["README.md", "src/demo_evidence/__init__.py"]
required_wheel_entries = ["demo_evidence/__init__.py"]
forbidden_archive_prefixes = ["tests/", "docs/"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-runtime" / "pyproject.toml",
        """
[project]
name = "demo-runtime"
version = "0.1.0"
dependencies = ["numpy>=1.0"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.sdist]
include = ["README.md", "src/demo_runtime/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/demo_runtime"]
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-dev" / "pyproject.toml",
        """
[project]
name = "demo-dev"
version = "0.1.0"
dependencies = ["PyYAML>=6.0"]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.sdist]
include = ["README.md", "src/demo_dev/**"]

[tool.hatch.build.targets.wheel]
packages = ["src/demo_dev"]
""".strip()
        + "\n",
    )
    _write(repo_root / "packages" / "demo-runtime" / "README.md", "# Demo runtime\n")
    _write(repo_root / "packages" / "demo-dev" / "README.md", "# Demo dev\n")
    _write(
        repo_root
        / "packages"
        / "demo-runtime"
        / "src"
        / "demo_runtime"
        / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    _write(
        repo_root / "packages" / "demo-dev" / "src" / "demo_dev" / "__init__.py",
        "__version__ = '0.1.0'\n",
    )
    return repo_root


def test_load_publication_readiness_settings_reads_repo_owned_policy() -> None:
    settings = load_publication_readiness_settings(REPO_ROOT)

    assert settings["required_evidence_input_manifest"] == "inputs.manifest.json"
    assert settings["target_shape_packages"] == [
        "bijux-phylogenetics",
        "phylogenetic",
        "bijux-phylogenetics-dev",
    ]


def test_load_package_bundle_policies_reads_repo_owned_policy() -> None:
    policies = load_package_bundle_policies(REPO_ROOT)

    assert set(policies) == {
        "bijux-phylogenetics",
        "phylogenetic",
        "bijux-phylogenetics-dev",
    }
    assert policies["bijux-phylogenetics"].wheel_module_root == "bijux_phylogenetics"
    assert (
        "bijux_phylogenetics_dev/quality/package_install_smoke.py"
        in policies["bijux-phylogenetics-dev"].required_wheel_entries
    )


def test_runtime_publication_policy_allows_secured_xml_dependency() -> None:
    runtime_policy = load_package_bundle_policies(REPO_ROOT)["bijux-phylogenetics"]

    assert runtime_policy.allowed_dependencies == (
        "biopython>=1.87,<2.0",
        "cairosvg>=2.9.0,<3.0",
        "defusedxml>=0.7.1,<1.0",
        "PyYAML>=6.0,<7.0",
    )
    assert (
        "src/bijux_phylogenetics/resources/examples/alignments/example_alignment.fasta"
        in runtime_policy.required_sdist_entries
    )
    assert (
        "src/bijux_phylogenetics/resources/examples/trees/example_tree.nwk"
        in runtime_policy.required_sdist_entries
    )
    assert (
        "src/bijux_phylogenetics/resources/datasets/mammals/primate_comparative/tree.nwk"
        in runtime_policy.required_sdist_entries
    )
    assert (
        "src/bijux_phylogenetics/resources/datasets/mammals/primate_comparative/traits.csv"
        in runtime_policy.required_sdist_entries
    )
    assert (
        "src/bijux_phylogenetics/resources/datasets/pathogens/rabies_cross_host_geography_panel/workflow-config.json"
        in runtime_policy.required_sdist_entries
    )
    assert (
        "bijux_phylogenetics/resources/examples/alignments/example_alignment.fasta"
        in runtime_policy.required_wheel_entries
    )
    assert (
        "bijux_phylogenetics/resources/examples/trees/example_tree.nwk"
        in runtime_policy.required_wheel_entries
    )
    assert (
        "bijux_phylogenetics/resources/datasets/mammals/primate_comparative/tree.nwk"
        in runtime_policy.required_wheel_entries
    )
    assert (
        "bijux_phylogenetics/resources/datasets/mammals/primate_comparative/traits.csv"
        in runtime_policy.required_wheel_entries
    )
    assert (
        "bijux_phylogenetics/resources/datasets/pathogens/rabies_cross_host_geography_panel/workflow-config.json"
        in runtime_policy.required_wheel_entries
    )


@pytest.mark.slow
def test_check_package_bundles_reports_target_package_policy_coverage(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = check_package_bundles(repo_root)

    assert report["issue_count"] == 0
    assert report["target_package_policies"][0]["package_name"] == "demo-runtime"
    assert report["target_package_policies"][-1]["package_name"] == "demo-evidence"
    assert report["target_package_policies"][-1]["has_target_policy"] is True


def test_build_dependency_policy_report_is_clean_for_repository() -> None:
    report = build_dependency_policy_report(REPO_ROOT)

    assert report["issue_count"] == 0
    assert report["package_count"] == 3


def test_build_dependency_policy_report_flags_dependency_drift(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "packages" / "demo-runtime" / "pyproject.toml",
        "[project]\nname='demo-runtime'\ndependencies=['pandas>=2.0']\n",
    )

    report = build_dependency_policy_report(repo_root)

    assert report["issue_count"] == 1
    assert report["issues"][0]["code"] == "dependency-policy-drift"


def test_audit_package_bundle_directory_flags_missing_and_forbidden_entries(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    policy = load_package_bundle_policies(repo_root)["demo-runtime"]
    dist_dir = tmp_path / "dist"
    _write_wheel(
        dist_dir / "demo_runtime-0.1.0-py3-none-any.whl",
        {
            "demo_runtime/__init__.py": "__version__='0.1.0'\n",
            "tests/test_runtime.py": "bad\n",
        },
    )
    _write_sdist(
        dist_dir / "demo_runtime-0.1.0.tar.gz",
        {
            "README.md": "# Demo\n",
            "tests/test_runtime.py": "bad\n",
        },
    )

    report = audit_package_bundle_directory(policy, dist_dir)

    issue_codes = {issue["code"] for issue in report["issues"]}
    assert issue_codes == {"missing-sdist-entry", "forbidden-archive-prefix"}


def test_audit_package_bundle_directory_rejects_linked_sdist_entries(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    policy = load_package_bundle_policies(repo_root)["demo-runtime"]
    dist_dir = tmp_path / "dist"
    _write_wheel(
        dist_dir / "demo_runtime-0.1.0-py3-none-any.whl",
        {
            "demo_runtime/__init__.py": "__version__='0.1.0'\n",
        },
    )
    _write_sdist_with_symbolic_link(
        dist_dir / "demo_runtime-0.1.0.tar.gz",
        members={
            "README.md": "# Demo\n",
            "src/demo_runtime/__init__.py": "__version__='0.1.0'\n",
        },
        linked_member="LICENSE",
        link_target="../../LICENSE",
    )

    report = audit_package_bundle_directory(policy, dist_dir)

    assert report["issue_count"] == 1
    assert report["issues"][0]["code"] == "linked-archive-entry"
    assert "../../LICENSE" in report["issues"][0]["message"]


def test_audit_package_bundle_directory_accepts_compliant_archives(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    policy = load_package_bundle_policies(repo_root)["demo-runtime"]
    dist_dir = tmp_path / "dist"
    _write_wheel(
        dist_dir / "demo_runtime-0.1.0-py3-none-any.whl",
        {"demo_runtime/__init__.py": "__version__='0.1.0'\n"},
    )
    _write_sdist(
        dist_dir / "demo_runtime-0.1.0.tar.gz",
        {
            "README.md": "# Demo\n",
            "src/demo_runtime/__init__.py": "__version__='0.1.0'\n",
        },
    )

    report = audit_package_bundle_directory(policy, dist_dir)

    assert report["issue_count"] == 0


@pytest.mark.slow
def test_check_package_bundles_builds_and_audits_publishable_packages(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    report = check_package_bundles(repo_root)

    assert report["issue_count"] == 0
    assert report["package_count"] == 2


@pytest.mark.slow
def test_check_package_bundles_rebuilds_cleanly_into_staged_output_directories(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)

    first_report = check_package_bundles(repo_root)
    second_report = check_package_bundles(repo_root)

    assert first_report["issue_count"] == 0
    assert second_report["issue_count"] == 0
    assert second_report["package_count"] == 2


@pytest.mark.slow
def test_check_package_bundles_keeps_published_bundle_outputs_unchanged(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    artifacts_root = repo_root / "artifacts" / "root" / "package-bundles"

    published_report = build_package_bundle_report(
        repo_root,
        artifacts_root=artifacts_root,
        build_artifacts=True,
    )
    sentinel = artifacts_root / "demo-runtime" / "sentinel.txt"
    sentinel.write_text("keep\n", encoding="utf-8")

    checked_report = check_package_bundles(repo_root, artifacts_root=artifacts_root)

    assert published_report["issue_count"] == 0
    assert checked_report["issue_count"] == 0
    assert sentinel.read_text(encoding="utf-8") == "keep\n"
