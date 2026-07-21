from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import tomllib
from typing import cast

REPO_ROOT = Path(__file__).resolve().parents[3]
DOCS_ROOT = REPO_ROOT / "docs"
PRODUCT_ROOT = DOCS_ROOT / "01-bijux-phylogenetics"
DOMAINS_ROOT = DOCS_ROOT / "02-bijux-phylogenetics-domains"
EVIDENCE_GUIDE_ROOT = DOCS_ROOT / "03-bijux-phylogenetics-evidence-book"
PARITIES_ROOT = DOCS_ROOT / "04-bijux-phylogenetics-parities"
MAINTAINER_ROOT = DOCS_ROOT / "05-bijux-phylogenetics-maintain"
EVIDENCE_RECORD_ROOT = REPO_ROOT / "evidence-book" / "studies"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_toml(path: Path) -> dict[str, object]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _verdict_counts(study_id: str) -> Counter[str]:
    statuses: list[str] = []
    for path in sorted(
        (EVIDENCE_RECORD_ROOT / study_id).glob("evidence-*/manifest.json")
    ):
        payload = json.loads(_read(path))
        statuses.append(cast(dict[str, str], payload["verdict"])["status"])
    return Counter(statuses)


def test_documentation_hierarchy_has_canonical_domain_and_redirect_paths() -> None:
    assert PRODUCT_ROOT.is_dir()
    assert DOMAINS_ROOT.is_dir()
    assert EVIDENCE_GUIDE_ROOT.is_dir()
    assert PARITIES_ROOT.is_dir()
    assert MAINTAINER_ROOT.is_dir()
    assert not (DOCS_ROOT / "02-bijux-phylogenetics-evidence-book").exists()
    assert not (DOCS_ROOT / "03-bijux-phylogenetics-maintain").exists()
    assert not (DOCS_ROOT / "04-bijux-phylogenetics-maintain").exists()

    config = _read(REPO_ROOT / "mkdocs.yml")
    assert (
        "02-bijux-phylogenetics-evidence-book/index.md: "
        "03-bijux-phylogenetics-evidence-book/index.md"
    ) in config
    assert (
        "03-bijux-phylogenetics-maintain/index.md: "
        "05-bijux-phylogenetics-maintain/index.md"
    ) in config
    assert (
        "04-bijux-phylogenetics-maintain/index.md: "
        "05-bijux-phylogenetics-maintain/index.md"
    ) in config


def test_root_readme_publishes_native_external_and_evidence_boundaries() -> None:
    readme = _read(REPO_ROOT / "README.md")

    assert "infer_nucleotide_maximum_likelihood_result" in readme
    assert "run_bayesian_inference" in readme
    assert "MAFFT" in readme
    assert "external engines" in readme
    assert "Scientific domains" in readme
    assert "Evidence Book" in readme
    assert "evidence-backed" in readme
    assert "bijux-phylogenetics demo run" in readme
    assert "artifacts/capability-demo/manifest.json" in readme
    assert "Read A Result At The Right Layer" in readme
    assert "Evidence Book bundle" in readme


def test_site_landing_connects_first_run_to_scientific_review() -> None:
    landing = _read(DOCS_ROOT / "index.md")

    assert "Run The Installed Runtime" in landing
    assert "Read the\nmanifest before the reports" in landing
    assert "Follow A Scientific Question" in landing
    assert "it is not a numerical inference\nbenchmark" in landing
    assert "Locate The Authoritative Record" in landing
    assert (
        "Runtime parity and Evidence Book status use different vocabularies" in landing
    )


def test_foundation_and_architecture_define_owned_runtime_contracts() -> None:
    foundation = _read(PRODUCT_ROOT / "foundation" / "index.md")
    scope = _read(PRODUCT_ROOT / "foundation" / "repository-scope-and-limits.md")
    claims = _read(PRODUCT_ROOT / "foundation" / "claim-reading-guide.md")
    evidence_boundary = _read(PRODUCT_ROOT / "foundation" / "evidence-boundary.md")
    ownership = _read(PRODUCT_ROOT / "foundation" / "product-surface-and-ownership.md")
    architecture = _read(PRODUCT_ROOT / "architecture" / "index.md")
    system_map = _read(PRODUCT_ROOT / "architecture" / "runtime-system-map.md")
    evidence_consumer = _read(
        PRODUCT_ROOT / "architecture" / "evidence-consumer-boundary.md"
    )

    assert "finite-state likelihood" in foundation
    assert "supported Bayesian DNA inference" in foundation
    assert "Contract Propagation" in foundation
    assert "Refusal Is Part Of The Contract" in foundation
    assert "Native And External Computation" in scope
    assert "Responsibility At Failure Boundaries" in scope
    assert "Capability And Assurance Are Independent Axes" in scope
    assert "Qualifiers That Change Meaning" in claims
    assert "Combining Proofs" in claims
    assert "Promotion Requirements" in evidence_boundary
    assert "Valid Boundary States" in evidence_boundary
    assert "Ownership Matrix" in ownership
    assert "Dependency Direction" in architecture
    assert "Layer Contracts" in architecture
    assert "Scientific State Transformations" in architecture
    assert "Failure Attribution" in system_map
    assert "State Retained Across Boundaries" in system_map
    assert "Partial And Refused Runs" in system_map
    assert "Dependency Rule" in evidence_consumer
    assert "Consumption Contract" in evidence_consumer
    assert "Material Changes And Freshness" in evidence_consumer


def test_product_handbook_routes_analysis_and_preserves_run_identity() -> None:
    product = _read(PRODUCT_ROOT / "index.md")
    operations = _read(PRODUCT_ROOT / "operations" / "index.md")

    assert "Select An Analysis Route" in product
    assert "A Complete Result Is More Than An Estimate" in product
    assert "Locate The Owning Contract" in product
    assert "parity registry, observation rows, and failure artifacts" in product
    assert "claim identifier, source, inputs, checks, freshness, verdict" in product
    assert "Establish The Run Identity First" in operations
    assert "Minimum Handoff Record" in operations
    assert "Changing presentation or omitting the failed check" in operations


def test_programmatic_interfaces_name_curated_and_frozen_contracts() -> None:
    api_surface = _read(PRODUCT_ROOT / "interfaces" / "api-surface.md")
    native_surface = _read(
        PRODUCT_ROOT / "interfaces" / "native-inference-and-benchmarks.md"
    )
    artifacts = _read(PRODUCT_ROOT / "interfaces" / "artifact-consumption-guide.md")

    for name in (
        "run_sequence_to_tree_workflow",
        "run_comparative_model_workflow",
        "run_ancestral_reconstruction_workflow",
        "run_configured_phylo_workflow",
    ):
        assert name in api_surface
    assert "freeze-only schema promise" in api_surface
    assert "NATIVE_MAXIMUM_LIKELIHOOD_ONLY_ENV_VAR" in api_surface
    assert "benchmark_native_maximum_likelihood_suite" in native_surface
    assert "expected-assumption-difference" in native_surface
    assert "Consume a workflow bundle from its manifest outward" in artifacts


def test_interface_guides_preserve_failure_compatibility_and_completeness() -> None:
    interface_index = _read(PRODUCT_ROOT / "interfaces" / "index.md")
    selection = _read(PRODUCT_ROOT / "interfaces" / "surface-selection.md")
    contracts = _read(PRODUCT_ROOT / "interfaces" / "artifact-contracts.md")
    consumption = _read(PRODUCT_ROOT / "interfaces" / "artifact-consumption-guide.md")

    assert "Contract Comparison" in interface_index
    assert "Compatibility Questions" in interface_index
    assert "Choose By Failure And Recovery" in selection
    assert "Migration Test" in selection
    assert "Completeness And Status" in contracts
    assert "Integrity Boundaries" in contracts
    assert "Establish Provenance Before Values" in consumption
    assert "Handle Incomplete Bundles" in consumption


def test_native_inference_guides_preserve_search_and_chain_diagnostics() -> None:
    maximum_likelihood = _read(
        PRODUCT_ROOT / "operations" / "native-maximum-likelihood-workflows.md"
    )
    bayesian = _read(PRODUCT_ROOT / "operations" / "native-bayesian-workflows.md")
    boundaries = _read(PRODUCT_ROOT / "operations" / "operational-boundaries.md")

    assert "start_tree_count" in maximum_likelihood
    assert "accepted search trace" in maximum_likelihood
    assert "bootstrap_replicate_count" in maximum_likelihood
    assert "run_bayesian_inference" in bayesian
    assert "effective sample size" in bayesian
    assert "distinct_topology_count" in bayesian
    assert "Failure Handling" in boundaries
    assert "Publication Threshold" in boundaries


def test_domain_map_covers_scientific_ownership_and_invariants() -> None:
    domain_index = _read(DOMAINS_ROOT / "index.md")
    data_index = _read(DOMAINS_ROOT / "data-and-structure" / "index.md")
    inference_index = _read(DOMAINS_ROOT / "inference" / "index.md")
    evolution_index = _read(DOMAINS_ROOT / "evolutionary-analysis" / "index.md")
    assurance_index = _read(DOMAINS_ROOT / "scientific-assurance" / "index.md")
    tree_domain = _read(
        DOMAINS_ROOT / "data-and-structure" / "trees-alignments-and-taxa.md"
    )
    inference_domain = _read(DOMAINS_ROOT / "inference" / "likelihood-and-inference.md")
    comparative_domain = _read(
        DOMAINS_ROOT / "evolutionary-analysis" / "comparative-and-ancestral-methods.md"
    )
    parsimony_domain = _read(DOMAINS_ROOT / "inference" / "parsimony.md")
    simulation_domain = _read(
        DOMAINS_ROOT / "scientific-assurance" / "simulation-and-datasets.md"
    )
    assurance_domain = _read(
        DOMAINS_ROOT / "scientific-assurance" / "validation-benchmarks-and-reporting.md"
    )

    for domain in (
        "Trees, alignments, and taxa",
        "Likelihood and inference",
        "Comparative and ancestral methods",
        "Parsimony",
        "Simulation and datasets",
        "Validation, benchmarks, and reporting",
    ):
        assert domain in domain_index
    assert "Taxon Reconciliation" in tree_domain
    assert "Tree Sets And Uncertainty" in tree_domain
    assert "Finite-State Likelihood" in inference_domain
    assert "Bayesian Inference" in inference_domain
    assert "PGLS And Trait Covariance" in comparative_domain
    assert "Ancestral Reconstruction" in comparative_domain
    assert "Searching Tree Space" in parsimony_domain
    assert "Support And Consensus" in parsimony_domain
    assert "Native Simulation Families" in simulation_domain
    assert "Preventing Simulation Circularity" in simulation_domain
    assert "Benchmark Design" in assurance_domain
    assert "From Benchmark To Evidence" in assurance_domain
    assert "Domain Handoffs" in domain_index
    assert "Cross-Domain Handoff Record" in domain_index
    assert "Attribute Failure To The Deciding Domain" in domain_index
    assert "Intake Contract" in data_index
    assert "Refuse Ambiguous Inputs" in data_index
    assert "Select The Inferential Contract" in inference_index
    assert "Review The Search, Not Only The Winner" in inference_index
    assert "Interpret The Run State" in inference_index
    assert "Completion is not a scientific verdict" in inference_index
    assert "Match The Question To The Model Family" in evolution_index
    assert "Sensitivity Axes" in evolution_index
    assert "Choose Evidence That Matches The Claim" in assurance_index
    assert "Preserve Negative Evidence" in assurance_index
    assert "Status Namespaces" in assurance_domain
    assert "live parity case" in assurance_domain
    assert "resolution workflow" in assurance_domain
    assert "matched, mismatched, not comparable, open" not in assurance_domain


def test_parity_handbook_separates_live_reference_and_adapter_boundaries() -> None:
    parity_index = _read(PARITIES_ROOT / "index.md")
    r_ecosystem = _read(PARITIES_ROOT / "r-ecosystem" / "index.md")
    engines = _read(PARITIES_ROOT / "external-engines" / "index.md")
    correspondence = _read(PARITIES_ROOT / "correspondence" / "index.md")

    for tool in ("APE", "GEIGER", "phytools", "MAFFT", "IQ-TREE2", "BEAST"):
        assert tool in parity_index
    assert "RUN_LIVE_PARITY=1" in r_ecosystem
    assert "computational owners" in engines
    assert "Runtime Case Status" in correspondence
    assert "Promotion To An Evidence Verdict" in correspondence
    assert "`passed`:" in correspondence
    assert "`failed`:" in correspondence
    assert "`skipped`:" in correspondence
    assert "--ape-shard load_tree" in correspondence
    assert "Three Records, Three Meanings" in parity_index
    assert "Matched · mismatched" not in parity_index


def test_r_parity_guides_publish_registry_scope_and_non_claims() -> None:
    r_root = PARITIES_ROOT / "r-ecosystem"
    ape = _read(r_root / "ape.md")
    geiger = _read(r_root / "geiger.md")
    phytools = _read(r_root / "phytools.md")

    assert "180 cases across 35" in ape
    assert "ape::dist.dna" not in ape
    assert "`dist.dna` registry" in ape
    assert "48 cases across 17" in geiger
    assert "optimizer triage" in geiger
    assert "65 cases across 21" in phytools
    assert "exact histories are not expected" in phytools


def test_engine_correspondence_guides_preserve_computational_ownership() -> None:
    engine_root = PARITIES_ROOT / "external-engines"
    alignment = _read(engine_root / "alignment-and-trimming.md")
    inference = _read(engine_root / "tree-inference.md")
    bayesian = _read(engine_root / "bayesian-inference.md")

    assert "MAFFT owns" in alignment
    assert "numerical parity claim" in alignment
    assert "incomplete bundles" in inference
    assert "remain different statistical" in inference
    assert "quantities and are labeled" in inference
    assert "prepared-only analysis" in bayesian
    assert "does not dispatch through MrBayes or BEAST" in bayesian


def test_correspondence_matrix_distinguishes_verified_and_open_surfaces() -> None:
    study_references = _read(PARITIES_ROOT / "r-ecosystem" / "study-references.md")
    matrix = _read(PARITIES_ROOT / "correspondence" / "verification-matrix.md")
    reproduction = _read(
        PARITIES_ROOT / "correspondence" / "artifacts-and-reproduction.md"
    )

    assert "no standalone live `nlme` registry" in study_references
    assert "tip count of 75" in study_references
    assert "plot-only" in study_references
    assert "no posterior parity claim" in matrix
    assert "External adapter" in matrix
    assert "Mode, Outcome, And Verdict" in matrix
    assert "not_comparable` verdicts and open posterior follow-through" in matrix
    assert "--manifest artifacts/parity/checked-manifest.json" in reproduction
    assert "Reconstruct The Denominator" in reproduction
    assert (
        "Preserve mismatches, skips, warnings, and incomplete outputs" in reproduction
    )


def test_quality_guides_separate_tests_benchmarks_and_release_evidence() -> None:
    quality_index = _read(PRODUCT_ROOT / "quality" / "index.md")
    lanes = _read(PRODUCT_ROOT / "quality" / "validation-lanes.md")
    honesty = _read(PRODUCT_ROOT / "quality" / "evidence-honesty-and-limits.md")
    test_strategy = _read(PRODUCT_ROOT / "quality" / "test-strategy.md")
    benchmark = _read(PRODUCT_ROOT / "quality" / "native-benchmark-evidence.md")
    release = _read(PRODUCT_ROOT / "quality" / "release-and-publication-gates.md")
    review = _read(PRODUCT_ROOT / "quality" / "release-review-workflow.md")

    assert "Quality Decision Record" in quality_index
    assert "Stronger Context Can Reveal New Failure" in quality_index
    assert "Select A Lane From The Claim" in lanes
    assert "mismatch_unexplained" in lanes
    assert "Governed Verdict Meanings" in honesty
    assert "Observation Status And Bundle Verdict" in honesty
    assert "engine_contract" in test_strategy
    assert "engine_real" in test_strategy
    assert "scientific_validation" in test_strategy
    assert "make package-verify" in test_strategy
    assert "native-bug" in benchmark
    assert "Missing And Failed Cases" in benchmark
    assert "make check-release-readiness" in release
    assert "Open Scientific Verdicts" in release
    assert "Required Evidence Set" in review
    assert "Review Procedure" in review
    assert "Release Blockers" in review


def test_maintainer_docs_define_ownership_and_native_publication() -> None:
    maintain_root = REPO_ROOT / "docs" / "05-bijux-phylogenetics-maintain"
    ownership = _read(
        maintain_root / "documentation-system" / "documentation-ownership.md"
    )
    native = _read(
        maintain_root / "documentation-system" / "native-surface-documentation.md"
    )

    assert "Resolve A Statement To Its Owner" in ownership
    assert "Duplication And Projection" in ownership
    assert "Native Surface Publication Contract" in native
    assert "Keep Capability And Assurance Separate" in native
    assert "Publication Blockers" in native
    assert "release story" not in native


def test_primate_guides_publish_values_and_open_boundaries() -> None:
    pcm1 = EVIDENCE_GUIDE_ROOT / "pcm1-primate-longevity-signal"
    pcm2 = EVIDENCE_GUIDE_ROOT / "pcm2-primate-pgls-and-signal"
    preparation = _read(pcm1 / "data-and-tree-preparation.md")
    signal = _read(pcm1 / "signal-and-ancestral-results.md")
    models = _read(pcm2 / "model-comparisons-and-pgls.md")
    boundaries = _read(pcm2 / "signal-diagnostics-and-open-boundaries.md")

    assert "tree–trait join" in preparation
    assert "evidence-009" in preparation
    assert "matched_with_tolerance" in signal
    assert "not_comparable" in signal
    assert "Pagel-λ PGLS" in models
    assert "likelihood-ratio" in models
    assert "intercept-mode likelihood sweep" in boundaries
    assert "boundary register" in boundaries


def test_evidence_guides_publish_verdict_distribution_and_model_semantics() -> None:
    evidence_index = _read(EVIDENCE_GUIDE_ROOT / "index.md")
    pcm3 = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm3-primate-continuous-and-categorical-regression"
        / "index.md"
    )
    pcm4 = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models"
        / "index.md"
    )
    pcm5 = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm5-cichlid-phylogenetic-generalized-linear-mixed-models"
        / "index.md"
    )

    assert "Current Verdict Ledger" in evidence_index
    assert "Evidence Strength Is Multidimensional" in evidence_index
    assert "Read Coefficients From The Design Matrix" in pcm3
    assert "All eight bundle manifests record `not_comparable`" in pcm3
    assert "Model Components" in pcm4
    assert "All seven bundle manifests record `not_comparable`" in pcm4
    assert "Preserve The Scale Of Every Quantity" in pcm5
    assert "All nine bundle manifests record `not_comparable`" in pcm5
    assert _verdict_counts("primate-longevity-signal") == Counter(
        {"matched": 8, "matched_with_tolerance": 1}
    )
    assert _verdict_counts("primate-pgls-and-signal") == Counter(
        {"matched": 3, "matched_with_tolerance": 6, "not_comparable": 1}
    )
    assert _verdict_counts("primate-continuous-and-categorical-regression") == Counter(
        {"not_comparable": 8}
    )
    assert _verdict_counts(
        "cichlid-bayesian-linear-and-phylogenetic-mixed-models"
    ) == Counter({"not_comparable": 7})
    assert _verdict_counts(
        "cichlid-phylogenetic-generalized-linear-mixed-models"
    ) == Counter({"not_comparable": 9})


def test_package_metadata_points_to_canonical_documentation_paths() -> None:
    runtime = cast(
        dict[str, object],
        _load_toml(REPO_ROOT / "packages" / "bijux-phylogenetics" / "pyproject.toml")[
            "project"
        ],
    )
    alias = cast(
        dict[str, object],
        _load_toml(REPO_ROOT / "packages" / "phylogenetic" / "pyproject.toml")[
            "project"
        ],
    )
    maintainer = cast(
        dict[str, object],
        _load_toml(
            REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "pyproject.toml"
        )["project"],
    )

    runtime_urls = cast(dict[str, str], runtime["urls"])
    alias_urls = cast(dict[str, str], alias["urls"])
    maintainer_urls = cast(dict[str, str], maintainer["urls"])
    assert runtime_urls["Evidence Book"].endswith(
        "/03-bijux-phylogenetics-evidence-book/"
    )
    assert alias_urls["Evidence Book"].endswith(
        "/03-bijux-phylogenetics-evidence-book/"
    )
    assert maintainer_urls["Documentation"].endswith(
        "/05-bijux-phylogenetics-maintain/"
    )


def test_package_readmes_define_installed_alias_and_governance_contracts() -> None:
    runtime = _read(REPO_ROOT / "packages" / "bijux-phylogenetics" / "README.md")
    alias = _read(REPO_ROOT / "packages" / "phylogenetic" / "README.md")
    maintainer = _read(REPO_ROOT / "packages" / "bijux-phylogenetics-dev" / "README.md")

    assert "Verify The Installed Artifact" in runtime
    assert "Select A Public Surface" in runtime
    assert "it is not an inference-accuracy result" in runtime
    assert "Verify The Alias End To End" in alias
    assert "Diagnose A Difference" in alias
    assert "It does not create a second implementation" in alias
    assert "What The Governed Commands Establish" in maintainer
    assert "Change-To-Check Map" in maintainer
    assert "make validate-evidence-book" in maintainer


def test_maintainer_handbook_defines_change_and_artifact_ownership() -> None:
    handbook = _read(MAINTAINER_ROOT / "index.md")
    ownership = _read(
        MAINTAINER_ROOT / "documentation-system" / "documentation-ownership.md"
    )
    contracts = _read(
        MAINTAINER_ROOT / "trust-surfaces" / "public-trust-surface-maintenance.md"
    )
    artifacts = _read(MAINTAINER_ROOT / "repository-operations" / "artifact-hygiene.md")
    release = _read(
        MAINTAINER_ROOT / "repository-operations" / "release-surface-alignment.md"
    )

    assert "Change Classification" in handbook
    assert "Maintained Invariants" in handbook
    assert "Canonicality Rules" in ownership
    assert "Change Impact Review" in contracts
    assert "Contract Change Classes" in contracts
    assert "Completeness Audit" in contracts
    assert "Classify Before Writing" in artifacts
    assert "make check-artifact-governance" in artifacts
    assert "Safe Cleanup" in artifacts
    assert "Alignment Matrix" in release
    assert "make check-release-readiness" in release
    assert "Release Blockers" in release
    assert "narrative alignment" not in release


def test_maintainer_indexes_define_operational_decisions_without_meta_prose() -> None:
    documentation = _read(MAINTAINER_ROOT / "documentation-system" / "index.md")
    trust = _read(MAINTAINER_ROOT / "trust-surfaces" / "index.md")
    operations = _read(MAINTAINER_ROOT / "repository-operations" / "index.md")

    assert "Resolve Ownership Before Editing" in documentation
    assert "Completion Criteria" in documentation
    assert "Detect Cross-Surface Drift" in trust
    assert "Release Review Questions" in trust
    assert "Classify Every Output" in operations
    assert "Operational Completion" in operations
    for page in (documentation, trust, operations):
        assert "Why This Section Matters" not in page
        assert "public story" not in page


def test_public_markdown_avoids_obsolete_meta_language() -> None:
    public_roots = (PRODUCT_ROOT, DOMAINS_ROOT, EVIDENCE_GUIDE_ROOT)
    forbidden = (
        "public story",
        "what this section has to protect",
        "this page exists so",
        "why this page matters",
    )

    for root in public_roots:
        for path in root.rglob("*.md"):
            text = _read(path).lower()
            for phrase in forbidden:
                assert phrase not in text, f"{path}: obsolete phrase {phrase!r}"
