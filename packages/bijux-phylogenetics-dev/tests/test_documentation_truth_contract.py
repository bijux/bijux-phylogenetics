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


def test_public_operations_pages_publish_native_flagship_workflows() -> None:
    operations_index = _read(
        REPO_ROOT / "docs" / "public" / "phylogenetics" / "operations" / "index.md"
    )
    common_workflows = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "operations"
        / "common-workflows.md"
    )
    native_ml_workflows = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "operations"
        / "native-maximum-likelihood-workflows.md"
    )
    native_bayesian_workflows = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "operations"
        / "native-bayesian-workflows.md"
    )
    native_benchmark_review = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "operations"
        / "native-benchmark-review.md"
    )

    assert "native maximum-likelihood workflows" in operations_index
    assert "native Bayesian workflows" in operations_index
    assert "native benchmark review" in operations_index
    assert "native maximum-likelihood workflows" in common_workflows
    assert "native Bayesian workflows" in common_workflows
    assert "native benchmark review" in common_workflows
    assert "infer_nucleotide_maximum_likelihood_result_from_alignment" in native_ml_workflows
    assert "run_bayesian_inference" in native_bayesian_workflows
    assert "benchmark native-maximum-likelihood-suite --json" in native_benchmark_review


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
    assert "maximum-likelihood" in runtime_project["keywords"]
    assert "bayesian-inference" in runtime_project["keywords"]
    assert runtime_project["urls"]["Native Workflows"].endswith("/operations/")
    assert runtime_project["urls"]["Evidence Book"].endswith(
        "/public/phylogenetics-evidence-book/"
    )
    assert runtime_project["urls"]["Benchmarks"].endswith(
        "/operations/native-benchmark-review/"
    )
    assert "maximum-likelihood" in alias_project["keywords"]
    assert "bayesian-inference" in alias_project["keywords"]
    assert alias_project["urls"]["Native Workflows"].endswith("/operations/")
    assert alias_project["urls"]["Evidence Book"].endswith(
        "/public/phylogenetics-evidence-book/"
    )
    assert alias_project["urls"]["Benchmarks"].endswith(
        "/operations/native-benchmark-review/"
    )


def test_package_readmes_link_package_roles_to_native_runtime_and_docs() -> None:
    runtime_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "README.md"
    )
    alias_readme = _read(REPO_ROOT / "packages" / "phylogenetic" / "README.md")
    dev_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "README.md"
    )

    assert "Native Inference And Benchmark Surfaces" in runtime_readme
    assert "Native maximum-likelihood workflows" in runtime_readme
    assert "Native Bayesian workflows" in runtime_readme
    assert "Native benchmark review" in runtime_readme
    assert "benchmark, workflow, and reporting surfaces" in alias_readme
    assert "Native maximum-likelihood workflows" in alias_readme
    assert "Native Bayesian workflows" in alias_readme
    assert "Native benchmark review" in alias_readme
    assert "Documentation ownership" in dev_readme
    assert "Native surface documentation" in dev_readme
    assert "Release surface alignment" in dev_readme


def test_internal_docs_name_native_surface_documentation_rule() -> None:
    maintainer_index = _read(
        REPO_ROOT / "docs" / "internal" / "maintain" / "index.md"
    )
    ownership = _read(
        REPO_ROOT
        / "docs"
        / "internal"
        / "maintain"
        / "documentation-ownership.md"
    )
    native_surface_documentation = _read(
        REPO_ROOT
        / "docs"
        / "internal"
        / "maintain"
        / "native-surface-documentation.md"
    )

    assert "native surface documentation" in maintainer_index
    assert (
        "flagship native surfaces should be visible across reference, workflow, and"
        in ownership
    )
    assert "Minimum Public Coverage For One Flagship Native Surface" in (
        native_surface_documentation
    )
