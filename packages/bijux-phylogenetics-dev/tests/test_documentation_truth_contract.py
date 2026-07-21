from __future__ import annotations

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


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _load_toml(path: Path) -> dict[str, object]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


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


def test_site_landing_connects_first_run_to_scientific_review() -> None:
    landing = _read(DOCS_ROOT / "index.md")

    assert "Run The Installed Runtime" in landing
    assert "Read the\nmanifest before the reports" in landing
    assert "Follow A Scientific Question" in landing
    assert "it is not a numerical inference\nbenchmark" in landing


def test_foundation_and_architecture_define_owned_runtime_contracts() -> None:
    foundation = _read(PRODUCT_ROOT / "foundation" / "index.md")
    ownership = _read(PRODUCT_ROOT / "foundation" / "product-surface-and-ownership.md")
    architecture = _read(PRODUCT_ROOT / "architecture" / "index.md")
    system_map = _read(PRODUCT_ROOT / "architecture" / "runtime-system-map.md")

    assert "finite-state likelihood" in foundation
    assert "supported Bayesian DNA inference" in foundation
    assert "Native And External Computation" in _read(
        PRODUCT_ROOT / "foundation" / "repository-scope-and-limits.md"
    )
    assert "Ownership Matrix" in ownership
    assert "Dependency Direction" in architecture
    assert "Failure Attribution" in system_map


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


def test_parity_handbook_separates_live_reference_and_adapter_boundaries() -> None:
    parity_index = _read(PARITIES_ROOT / "index.md")
    r_ecosystem = _read(PARITIES_ROOT / "r-ecosystem" / "index.md")
    engines = _read(PARITIES_ROOT / "external-engines" / "index.md")
    correspondence = _read(PARITIES_ROOT / "correspondence" / "index.md")

    for tool in ("APE", "GEIGER", "phytools", "MAFFT", "IQ-TREE2", "BEAST"):
        assert tool in parity_index
    assert "RUN_LIVE_PARITY=1" in r_ecosystem
    assert "computational owners" in engines
    assert "matched with tolerance" in correspondence
    assert "--ape-shard load_tree" in correspondence


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
    assert "--manifest artifacts/parity/checked-manifest.json" in reproduction
    assert (
        "Preserve mismatches, skips, warnings, and incomplete outputs" in reproduction
    )


def test_quality_guides_separate_tests_benchmarks_and_release_evidence() -> None:
    test_strategy = _read(PRODUCT_ROOT / "quality" / "test-strategy.md")
    benchmark = _read(PRODUCT_ROOT / "quality" / "native-benchmark-evidence.md")
    release = _read(PRODUCT_ROOT / "quality" / "release-and-publication-gates.md")

    assert "engine_contract" in test_strategy
    assert "engine_real" in test_strategy
    assert "scientific_validation" in test_strategy
    assert "make package-verify" in test_strategy
    assert "native-bug" in benchmark
    assert "Missing And Failed Cases" in benchmark
    assert "make check-release-readiness" in release
    assert "Open Scientific Verdicts" in release


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


def test_maintainer_handbook_defines_change_and_artifact_ownership() -> None:
    handbook = _read(MAINTAINER_ROOT / "index.md")
    ownership = _read(
        MAINTAINER_ROOT / "documentation-system" / "documentation-ownership.md"
    )
    contracts = _read(
        MAINTAINER_ROOT / "trust-surfaces" / "public-trust-surface-maintenance.md"
    )
    artifacts = _read(MAINTAINER_ROOT / "repository-operations" / "artifact-hygiene.md")

    assert "Change Classification" in handbook
    assert "Maintained Invariants" in handbook
    assert "Canonicality Rules" in ownership
    assert "Change Impact Review" in contracts
    assert "artifacts/" in artifacts


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
