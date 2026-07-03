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


def test_public_interface_guides_explain_surface_selection_and_artifact_use() -> None:
    interfaces_index = _read(
        REPO_ROOT / "docs" / "public" / "phylogenetics" / "interfaces" / "index.md"
    )
    surface_selection = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "interfaces"
        / "surface-selection.md"
    )
    artifact_contracts = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "interfaces"
        / "artifact-contracts.md"
    )
    artifact_consumption = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "interfaces"
        / "artifact-consumption-guide.md"
    )

    assert "surface selection guide" in interfaces_index
    assert "artifact consumption guide" in interfaces_index
    assert "I want one governed end-user workflow or report artifact" in (
        surface_selection
    )
    assert "I want to judge trust, parity, freshness, or evidence closure" in (
        surface_selection
    )
    assert "Stable Families Worth Depending On" in artifact_contracts
    assert "artifact consumption guide" in artifact_contracts
    assert "do not scrape terminal text when the workflow already emits JSON or TSV" in (
        artifact_consumption
    )


def test_public_architecture_pages_describe_package_roles_and_runtime_families() -> (
    None
):
    architecture_index = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "architecture"
        / "index.md"
    )
    package_map = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "architecture"
        / "package-map.md"
    )
    runtime_system_map = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "architecture"
        / "runtime-system-map.md"
    )
    evidence_boundary = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "architecture"
        / "evidence-consumer-boundary.md"
    )

    assert "which package should I install or depend on" in architecture_index
    assert "`bijux-phylogenetics` is the public center of gravity" in package_map
    assert "`phylogenetic` is the shorter compatibility install surface" in package_map
    assert "`bijux-phylogenetics-dev` is repository-owned maintainer tooling" in package_map
    assert "likelihood and inference" in runtime_system_map
    assert "native maximum-likelihood and native Bayesian DNA workflows" in (
        runtime_system_map
    )
    assert "runtime question" in evidence_boundary
    assert "evidence question" in evidence_boundary


def test_public_evidence_and_release_guides_avoid_fragile_summary_claims() -> None:
    evidence_index = _read(
        REPO_ROOT / "docs" / "public" / "phylogenetics-evidence-book" / "index.md"
    )
    evidence_reading_guide = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "reading-guide.md"
    )
    release_review_workflow = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "quality"
        / "release-review-workflow.md"
    )
    release_gates = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "quality"
        / "release-and-publication-gates.md"
    )

    assert "evidence book reading guide" in evidence_index
    assert "Treat those checked-in reports as the current source of truth" in (
        evidence_index
    )
    assert "teaching studies: `0`" not in evidence_index
    assert "migration studies: `0`" not in evidence_index
    assert "Core Reading Order" in evidence_reading_guide
    assert "make report-publish-readiness" in release_review_workflow
    assert "make check-release-readiness" in release_review_workflow
    assert "bijux-phylogenetics report release-truth" in release_review_workflow
    assert "Governed Release Surfaces" in release_gates
    assert "docs say more than metadata" in release_gates


def test_public_quality_and_foundation_guides_explain_trust_boundaries() -> None:
    quality_index = _read(
        REPO_ROOT / "docs" / "public" / "phylogenetics" / "quality" / "index.md"
    )
    validation_lanes = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "quality"
        / "validation-lanes.md"
    )
    test_strategy = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "quality"
        / "test-strategy.md"
    )
    foundation_index = _read(
        REPO_ROOT / "docs" / "public" / "phylogenetics" / "foundation" / "index.md"
    )
    claim_reading_guide = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "foundation"
        / "claim-reading-guide.md"
    )
    evidence_boundary = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics"
        / "foundation"
        / "evidence-boundary.md"
    )

    assert "validation lanes" in quality_index
    assert "runtime behavior lanes" in validation_lanes
    assert "scientific-validation lanes" in validation_lanes
    assert "a runtime unit or integration test does not prove broader scientific closure" in (
        validation_lanes
    )
    assert "validation lanes" in test_strategy
    assert "claim reading guide" in foundation_index
    assert "native ownership claim" in claim_reading_guide
    assert "wrapper claim" in claim_reading_guide
    assert "a documented runtime surface is not automatically evidence-closed" in (
        claim_reading_guide
    )
    assert "claim reading guide" in evidence_boundary


def test_public_pcm_study_guides_name_scope_order_and_boundaries() -> None:
    pcm1_index = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm1-primate-longevity-signal"
        / "index.md"
    )
    pcm1_preparation = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm1-primate-longevity-signal"
        / "data-and-tree-preparation.md"
    )
    pcm1_results = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm1-primate-longevity-signal"
        / "signal-and-ancestral-results.md"
    )
    pcm2_index = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm2-primate-pgls-and-signal"
        / "index.md"
    )
    pcm2_models = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm2-primate-pgls-and-signal"
        / "model-comparisons-and-pgls.md"
    )
    pcm2_boundaries = _read(
        REPO_ROOT
        / "docs"
        / "public"
        / "phylogenetics-evidence-book"
        / "pcm2-primate-pgls-and-signal"
        / "signal-diagnostics-and-open-boundaries.md"
    )

    assert "Read PCM1 In Order" in pcm1_index
    assert "What PCM1 Does Not Prove" in pcm1_index
    assert "What This Part Answers" in pcm1_preparation
    assert "What To Keep In Mind" in pcm1_results
    assert "Read PCM2 In Order" in pcm2_index
    assert "What PCM2 Does Not Prove" in pcm2_index
    assert "What This Part Answers" in pcm2_models
    assert "Reading Rule" in pcm2_boundaries
    assert "open boundaries are part of the study's value" in pcm2_boundaries


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
    assert runtime_project["urls"]["Interface Guide"].endswith(
        "/interfaces/surface-selection/"
    )
    assert runtime_project["urls"]["Artifact Guide"].endswith(
        "/interfaces/artifact-consumption-guide/"
    )
    assert runtime_project["urls"]["Release Guide"].endswith(
        "/quality/release-review-workflow/"
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
    assert alias_project["urls"]["Interface Guide"].endswith(
        "/interfaces/surface-selection/"
    )
    assert alias_project["urls"]["Artifact Guide"].endswith(
        "/interfaces/artifact-consumption-guide/"
    )
    assert alias_project["urls"]["Release Guide"].endswith(
        "/quality/release-review-workflow/"
    )


def test_package_readmes_link_package_roles_to_native_runtime_and_docs() -> None:
    runtime_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics" / "README.md"
    )
    alias_readme = _read(REPO_ROOT / "packages" / "phylogenetic" / "README.md")
    dev_readme = _read(
        REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "README.md"
    )
    root_readme = _read(REPO_ROOT / "README.md")

    assert "Surface selection guide" in root_readme
    assert "Artifact consumption guide" in root_readme
    assert "Validation lanes" in root_readme
    assert "Claim reading guide" in root_readme
    assert "Surface selection guide" in runtime_readme
    assert "Artifact consumption guide" in runtime_readme
    assert "Native Inference And Benchmark Surfaces" in runtime_readme
    assert "Native maximum-likelihood workflows" in runtime_readme
    assert "Native Bayesian workflows" in runtime_readme
    assert "Native benchmark review" in runtime_readme
    assert "Runtime architecture" in runtime_readme
    assert "Release review workflow" in runtime_readme
    assert "Evidence-book guide" in runtime_readme
    assert "benchmark, workflow, and reporting surfaces" in alias_readme
    assert "Surface selection guide" in alias_readme
    assert "Artifact consumption guide" in alias_readme
    assert "Native maximum-likelihood workflows" in alias_readme
    assert "Native Bayesian workflows" in alias_readme
    assert "Native benchmark review" in alias_readme
    assert "Runtime architecture" in alias_readme
    assert "Release review workflow" in alias_readme
    assert "Evidence-book guide" in alias_readme
    assert "Documentation ownership" in dev_readme
    assert "Native surface documentation" in dev_readme
    assert "Public trust surface maintenance" in dev_readme
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


def test_internal_docs_name_public_trust_surface_maintenance_rule() -> None:
    internal_index = _read(REPO_ROOT / "docs" / "internal" / "index.md")
    maintainer_index = _read(
        REPO_ROOT / "docs" / "internal" / "maintain" / "index.md"
    )
    trust_surface_maintenance = _read(
        REPO_ROOT
        / "docs"
        / "internal"
        / "maintain"
        / "public-trust-surface-maintenance.md"
    )

    assert "public trust surface maintenance" in internal_index
    assert "public trust surface maintenance" in maintainer_index
    assert "runtime architecture pages" in trust_surface_maintenance
    assert "The common failure is not one missing page. It is partial alignment" in (
        trust_surface_maintenance
    )


def test_internal_docs_name_public_interface_and_artifact_maintenance_rule() -> (
    None
):
    internal_index = _read(REPO_ROOT / "docs" / "internal" / "index.md")
    maintainer_index = _read(
        REPO_ROOT / "docs" / "internal" / "maintain" / "index.md"
    )
    interface_maintenance = _read(
        REPO_ROOT
        / "docs"
        / "internal"
        / "maintain"
        / "public-interface-and-artifact-maintenance.md"
    )

    assert "public interface and artifact maintenance" in internal_index
    assert "public interface and artifact maintenance" in maintainer_index
    assert "surface-selection guidance" in interface_maintenance
    assert "The usual failure is partial interface drift" in interface_maintenance


def test_internal_docs_name_trust_lane_and_study_guide_maintenance_rule() -> None:
    internal_index = _read(REPO_ROOT / "docs" / "internal" / "index.md")
    maintainer_index = _read(
        REPO_ROOT / "docs" / "internal" / "maintain" / "index.md"
    )
    trust_lane_maintenance = _read(
        REPO_ROOT
        / "docs"
        / "internal"
        / "maintain"
        / "trust-lane-and-study-guide-maintenance.md"
    )

    assert "trust lane and study guide maintenance" in internal_index
    assert "trust lane and study guide maintenance" in maintainer_index
    assert "quality pages that explain validation lanes" in trust_lane_maintenance
    assert "The failure here is trust compression" in trust_lane_maintenance
