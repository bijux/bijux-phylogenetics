from __future__ import annotations

from pathlib import Path
import tomllib


REPO_ROOT = Path(__file__).resolve().parents[3]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_toml(path: Path) -> dict[str, object]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def test_root_readme_describes_native_inference_without_old_false_boundary() -> None:
    readme = _read(REPO_ROOT / "README.md")

    assert "native maximum-likelihood tree inference results" in readme
    assert "native Bayesian public inference entry points" in readme
    assert (
        "a native maximum-likelihood or Bayesian inference engine implemented inside this repository"
        not in readme
    )


def test_public_runtime_foundation_names_native_runtime_scope_honestly() -> None:
    foundation_index = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "foundation"
        / "index.md"
    )

    assert "native likelihood and inference" in foundation_index
    assert "supported native Bayesian inference surfaces" in foundation_index
    assert "it does not claim to be a native maximum-likelihood or Bayesian inference" not in foundation_index


def test_public_native_inference_page_names_benchmark_and_bayesian_contracts() -> (
    None
):
    native_surface_page = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "interfaces"
        / "native-inference-and-benchmarks.md"
    )

    assert "benchmark_native_maximum_likelihood_suite" in native_surface_page
    assert "run_bayesian_inference" in native_surface_page
    assert "expected-assumption-difference" in native_surface_page
    assert "unsupported-case" in native_surface_page


def test_runtime_and_alias_metadata_publish_beta_native_runtime_story() -> None:
    runtime_pyproject = _load_toml(
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "pyproject.toml"
    )
    alias_pyproject = _load_toml(
        REPO_ROOT / "packages" / "phylogenetic" / "pyproject.toml"
    )

    runtime_project = runtime_pyproject["project"]
    alias_project = alias_pyproject["project"]

    assert (
        runtime_project["description"]
        == "Native phylogenetic runtime, comparative workflows, and evidence-linked reporting."
    )
    assert (
        alias_project["description"]
        == "Compatibility alias package for the canonical native phylogenetics runtime."
    )
    assert "Development Status :: 4 - Beta" in runtime_project["classifiers"]
    assert "Development Status :: 4 - Beta" in alias_project["classifiers"]


def test_package_readmes_link_package_roles_to_native_runtime_and_docs() -> None:
    runtime_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "README.md"
    )
    alias_readme = _read(REPO_ROOT / "packages" / "phylogenetic" / "README.md")
    dev_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "README.md"
    )

    assert "Native Inference And Benchmark Surfaces" in runtime_readme
    assert "benchmark, workflow, and reporting surfaces" in alias_readme
    assert "Documentation ownership" in dev_readme
    assert "Release surface alignment" in dev_readme
