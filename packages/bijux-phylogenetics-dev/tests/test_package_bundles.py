from __future__ import annotations

from pathlib import Path
import tarfile
import zipfile

from bijux_phylogenetics_dev.quality.package_bundles import (
    audit_package_bundle_directory,
    build_dependency_policy_report,
    load_package_bundle_policies,
)

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


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    _write(
        repo_root / "configs" / "publication_readiness.toml",
        """
[tool.bijux_phylogenetics.publication_readiness]
required_evidence_input_manifest = "inputs.manifest.json"
required_evidence_bundle_code_files = ["reference.R", "analysis.py"]
expected_publishable_packages = ["demo-runtime", "demo-dev"]
target_shape_packages = ["demo-runtime", "demo-dev", "demo-evidence"]
forbidden_runtime_subpackages = ["evidence"]
required_root_make_targets = ["check-package-bundles:"]

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
""".strip()
        + "\n",
    )
    _write(
        repo_root / "packages" / "demo-runtime" / "pyproject.toml",
        "[project]\nname='demo-runtime'\ndependencies=['numpy>=1.0']\n",
    )
    _write(
        repo_root / "packages" / "demo-dev" / "pyproject.toml",
        "[project]\nname='demo-dev'\ndependencies=['PyYAML>=6.0']\n",
    )
    return repo_root


def test_load_package_bundle_policies_reads_repo_owned_policy() -> None:
    policies = load_package_bundle_policies(REPO_ROOT)

    assert set(policies) == {
        "bijux-phylogenetics",
        "phylogenetic",
        "bijux-phylogenetics-dev",
    }
    assert policies["bijux-phylogenetics"].wheel_module_root == "bijux_phylogenetics"


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
