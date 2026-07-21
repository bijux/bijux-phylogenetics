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
    assert "Select An Execution Contract" in readme
    assert "the contract before choosing the command" in readme
    assert "Read The Demo As A Connected Record" in readme
    assert "phylo-inputs-report.json" in readme
    assert "Follow One Analysis From Data To Claim" in readme
    assert "capability proximity is not inherited evidence" in readme
    assert "Choose A Scientific Route" in readme
    assert "Start from the scientific object and required result" in readme


def test_site_landing_connects_first_run_to_scientific_review() -> None:
    landing = _read(DOCS_ROOT / "index.md")

    assert "Run The Installed Runtime" in landing
    assert "Read the\nmanifest before the reports" in landing
    assert "Follow A Scientific Question" in landing
    assert "it is not a numerical inference\nbenchmark" in landing
    assert "Locate The Authoritative Record" in landing
    assert "Build A Review Packet" in landing
    assert "missing outputs disappear from the denominator" in landing
    assert (
        "Runtime parity and Evidence Book status use different vocabularies" in landing
    )
    assert "Follow The Demo Records" in landing
    assert "Choose The Evidence Needed For The Decision" in landing
    assert "Trace A Result When Something Looks Wrong" in landing
    assert "changing a caption\nor omitting a row cannot resolve" in landing
    assert "Choose The Durable Record First" in landing
    assert "the selection denominator" in landing


def test_foundation_and_architecture_define_owned_runtime_contracts() -> None:
    foundation = _read(PRODUCT_ROOT / "foundation" / "index.md")
    scope = _read(PRODUCT_ROOT / "foundation" / "repository-scope-and-limits.md")
    claims = _read(PRODUCT_ROOT / "foundation" / "claim-reading-guide.md")
    evidence_boundary = _read(PRODUCT_ROOT / "foundation" / "evidence-boundary.md")
    ownership = _read(PRODUCT_ROOT / "foundation" / "product-surface-and-ownership.md")
    architecture = _read(PRODUCT_ROOT / "architecture" / "index.md")
    package_map = _read(PRODUCT_ROOT / "architecture" / "package-map.md")
    system_map = _read(PRODUCT_ROOT / "architecture" / "runtime-system-map.md")
    evidence_consumer = _read(
        PRODUCT_ROOT / "architecture" / "evidence-consumer-boundary.md"
    )

    assert "finite-state likelihood" in foundation
    assert "supported Bayesian DNA inference" in foundation
    assert "Contract Propagation" in foundation
    assert "Refusal Is Part Of The Contract" in foundation
    assert "Establish Object Identity Before Modeling" in foundation
    assert "Stop when" in foundation
    assert "Native And External Computation" in scope
    assert "Responsibility At Failure Boundaries" in scope
    assert "Capability And Assurance Are Independent Axes" in scope
    assert "Boundaries By Decision" in scope
    assert "without governed primary\nposterior outputs" in scope
    assert "Qualifiers That Change Meaning" in claims
    assert "Combining Proofs" in claims
    assert "Promotion Requirements" in evidence_boundary
    assert "Valid Boundary States" in evidence_boundary
    assert "Evidence Record Anatomy" in evidence_boundary
    assert "Do Not Promote By Proximity" in evidence_boundary
    assert "Changes That Reopen A Claim" in evidence_boundary
    assert "Ownership Matrix" in ownership
    assert "How To Recognize A Supported Surface" in ownership
    assert "Preserve The Computation Owner" in ownership
    assert "Dependency Direction" in architecture
    assert "Layer Contracts" in architecture
    assert "Scientific State Transformations" in architecture
    assert "Runtime Ownership Map" in architecture
    assert "Installed Identity Chain" in package_map
    assert "record both distribution versions" in package_map
    assert "Inspect Installed Identity" in package_map
    assert "Package Boundary Invariants" in package_map
    assert "What Ships In Each Distribution" in package_map
    assert "Reconstruct Distribution Identity" in package_map
    assert "Diagnose Boundary Failures" in package_map
    assert "Failure Attribution" in system_map
    assert "State Retained Across Boundaries" in system_map
    assert "Partial And Refused Runs" in system_map
    assert "Artifact Acceptance Gates" in system_map
    assert "The Identity Spine" in system_map
    assert "Consumer Precedence" in system_map
    assert "only successful outputs remain in the denominator" in system_map
    assert "Dependency Rule" in evidence_consumer
    assert "Consumption Contract" in evidence_consumer
    assert "Material Changes And Freshness" in evidence_consumer
    assert "Promotion Protocol" in evidence_consumer
    assert "Missing comparison data\nproduces `not_comparable`" in evidence_consumer


def test_product_handbook_routes_analysis_and_preserves_run_identity() -> None:
    product = _read(PRODUCT_ROOT / "index.md")
    operations = _read(PRODUCT_ROOT / "operations" / "index.md")

    assert "Select An Analysis Route" in product
    assert "A Complete Result Is More Than An Estimate" in product
    assert "Locate The Owning Contract" in product
    assert "parity registry, observation rows, and failure artifacts" in product
    assert "claim identifier, source, inputs, checks, freshness, verdict" in product
    assert "Make Four Decisions Before Execution" in product
    assert "Locate The Scientific Decision" in product
    assert "public interface starts a computation" in product
    assert "Define Completion Before Running" in product
    assert "A partial or refused result can be a complete record" in product
    assert "Open or blocked\n  resolution work" in product
    assert "Establish The Run Identity First" in operations
    assert "Minimum Handoff Record" in operations
    assert "Changing presentation or omitting the failed check" in operations
    assert "env inspect --json" in operations
    assert "Follow The Operational State Machine" in operations

    boundaries = _read(PRODUCT_ROOT / "operations" / "operational-boundaries.md")
    assert "Classify Before Retrying" in boundaries
    assert "Minimum Diagnostic Packet" in boundaries
    assert "Stop Conditions" in boundaries
    assert "widen tolerance after observing the result" in boundaries


def test_public_interfaces_define_supported_python_and_handoff_contracts() -> None:
    interfaces = PRODUCT_ROOT / "interfaces"
    index = _read(interfaces / "index.md")
    selection = _read(interfaces / "surface-selection.md")
    python = _read(interfaces / "python-surface.md")
    cli = _read(interfaces / "cli-surface.md")

    assert "Interface Invariants" in index
    assert "Preserve State Across Handoffs" in index
    assert "Silent coercion is a contract violation" in index
    assert "Worked Selections" in selection
    assert "Supported Import Layers" in python
    assert "Workflow API" in python
    assert "Compose A Native Workflow" in python
    assert "Compose An External-Engine Workflow" in python
    assert "Direct Scientific Modules" in python
    assert "Result Handling Contract" in python
    assert "Alias Imports" in python
    assert "Discover Without Guessing" in python
    assert "Snapshot A Reviewable Boundary" in python
    assert "Object identity is not artifact identity" in python
    assert "now also" not in python
    assert "bijux_phylogenetics.PhyloTree" not in python
    assert "Start With Live Discovery" in cli
    assert "Choose The Owning Family" in cli
    assert "Run An External Engine Without Erasing Ownership" in cli
    assert "Read Output Channels Deliberately" in cli
    assert "Keep Assurance Records Separate" in cli
    assert "Chain Commands Through Accepted Artifacts" in cli
    assert "later success erase the earlier failure denominator" in cli
    assert "commands --format json" in cli
    assert "comparative pgls" in cli
    assert "ancestral continuous" in cli
    assert "phylo run" in cli
    assert "adapter align" in cli
    assert "now lives" not in cli
    assert "now also" not in cli
    assert len(cli.splitlines()) < 400


def test_common_workflows_are_goal_oriented_and_match_registered_commands() -> None:
    operations = PRODUCT_ROOT / "operations"
    installation = _read(operations / "installation-and-setup.md")
    workflows = _read(operations / "common-workflows.md")

    assert "Run The Installed-Wheel Proof" in installation
    assert "phylo preflight" in installation
    assert "Admit The Environment Deliberately" in installation
    assert "installation health does not establish analysis readiness" in installation
    assert "Choose A Workflow" in workflows
    assert "comparative pgls" in workflows
    assert "ancestral continuous" in workflows
    assert "phylo preflight" in workflows
    assert "phylo run" in workflows
    assert "Read Status Before Values" in workflows
    assert "Put An Acceptance Gate Between Stages" in workflows
    assert "downstream success\ncannot repair an upstream contradiction" in workflows
    assert "ancestral reconstruct" not in workflows
    assert "env --json" not in workflows
    assert len(workflows.splitlines()) < 400


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
    assert "Select The Contract From The Consumer" in api_surface
    assert "ConfiguredPhyloWorkflowResult` deliberately has no" in api_surface
    assert "Failure And Recovery Contract" in api_surface
    assert "Cross Contracts By Meaning, Not Shape" in api_surface
    assert "matching keys" in api_surface
    assert "benchmark_native_maximum_likelihood_suite" in native_surface
    assert "expected-assumption-difference" in native_surface
    assert "Stable Native Boundary At A Glance" in native_surface
    assert "A Native Result Preserves The Decision Path" in native_surface
    assert "Choose The Narrowest Honest Claim" in native_surface
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
    assert "Semantic Compatibility" in contracts
    assert "External-Engine Provenance" in contracts
    assert "Producer Acceptance Criteria" in contracts
    assert "Establish Provenance Before Values" in consumption
    assert "Handle Incomplete Bundles" in consumption
    assert "Record The Consumer Acceptance Decision" in consumption
    assert "Do not edit the producer manifest" in consumption
    assert "Reconstruct The Denominator" in consumption
    assert "Join Ledgers By Declared Identity" in consumption
    assert "Produce A Review Handoff" in consumption


def test_native_inference_guides_preserve_search_and_chain_diagnostics() -> None:
    maximum_likelihood = _read(
        PRODUCT_ROOT / "operations" / "native-maximum-likelihood-workflows.md"
    )
    bayesian = _read(PRODUCT_ROOT / "operations" / "native-bayesian-workflows.md")
    boundaries = _read(PRODUCT_ROOT / "operations" / "operational-boundaries.md")

    assert "start_tree_count" in maximum_likelihood
    assert "accepted search trace" in maximum_likelihood
    assert "bootstrap_replicate_count" in maximum_likelihood
    assert "Choose The Inference Contract" in maximum_likelihood
    assert "load_nucleotide_maximum_likelihood_result_json" in maximum_likelihood
    assert "acceptance record" in maximum_likelihood
    assert "Understand The Search Pipeline" in maximum_likelihood
    assert "Read The Result Hierarchy" in maximum_likelihood
    assert "Compare Search Outcomes, Not Just Scores" in maximum_likelihood
    assert "run_bayesian_inference" in bayesian
    assert "effective sample size" in bayesian
    assert "distinct_topology_count" in bayesian
    assert "Read The Returned Record" in bayesian
    assert "Decide Whether The Run Is Interpretable" in bayesian
    assert "Observation Policy Changes The Likelihood" in bayesian
    assert "Read One Posterior Row" in bayesian
    assert "The Chain Report Owns Execution State" in bayesian
    assert "Failure Handling" in boundaries
    assert "Publication Threshold" in boundaries


def test_native_benchmark_review_exposes_executable_suite_and_dimensions() -> None:
    benchmark = _read(PRODUCT_ROOT / "operations" / "native-benchmark-review.md")

    assert "benchmark native-maximum-likelihood-suite" in benchmark
    assert "mkdir -p artifacts/native-benchmark" in benchmark
    assert "--manifest artifacts/native-benchmark/manifest.json" in benchmark
    assert "Read Suite Dimensions Separately" in benchmark
    assert "Minimum Publication Packet" in benchmark
    assert "native-bug" in benchmark


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
    ecology_domain = _read(
        DOMAINS_ROOT / "evolutionary-analysis" / "ecology-and-phylogeography.md"
    )
    parsimony_domain = _read(DOMAINS_ROOT / "inference" / "parsimony.md")
    distance_domain = _read(DOMAINS_ROOT / "inference" / "distance-based-inference.md")
    bayesian_domain = _read(DOMAINS_ROOT / "inference" / "bayesian-inference.md")
    simulation_domain = _read(
        DOMAINS_ROOT / "scientific-assurance" / "simulation-and-datasets.md"
    )
    assurance_domain = _read(
        DOMAINS_ROOT / "scientific-assurance" / "validation-benchmarks-and-reporting.md"
    )
    diagnostics_domain = _read(
        DOMAINS_ROOT / "scientific-assurance" / "diagnostics-and-sensitivity.md"
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
    assert "Runtime Locator Map" in tree_domain
    assert "reconciliation handoff" in tree_domain
    assert "Downstream Handoff Contracts" in tree_domain
    assert "Repair creates a new governed object" in tree_domain
    assert "Finite-State Likelihood" in inference_domain
    assert "Bayesian Inference" in inference_domain
    assert "Layered Runtime Surfaces" in inference_domain
    assert "Comparison Admissibility" in inference_domain
    assert "Enter At The Layer Your Question Requires" in inference_domain
    assert "Trace A Reported Likelihood" in inference_domain
    assert "PGLS And Trait Covariance" in comparative_domain
    assert "Ancestral Reconstruction" in comparative_domain
    assert "Method And Record Locators" in comparative_domain
    assert "Model Identity Before Values" in comparative_domain
    assert "Choose The Method From The Estimand" in comparative_domain
    assert "Publish A Comparative Result Packet" in comparative_domain
    assert "Host-Switching Records" in ecology_domain
    assert "summarize_continuous_phylogeography" in ecology_domain
    assert "Presentation Is Downstream" in ecology_domain
    assert "Searching Tree Space" in parsimony_domain
    assert "Support And Consensus" in parsimony_domain
    assert "Public Surface Map" in parsimony_domain
    assert "Minimum Parsimony Result Packet" in parsimony_domain
    assert "Establish Matrix Identity" in distance_domain
    assert "diagnose_distance_additivity" in distance_domain
    assert "write_distance_reproducibility_bundle" in distance_domain
    assert "Trace The Distance Result Lineage" in distance_domain
    assert "Refusal Boundaries" in distance_domain
    assert "Posterior Identity Begins Before Sampling" in bayesian_domain
    assert "Resume is continuation, not an independent replicate" in bayesian_domain
    assert "Native And External Boundaries" in bayesian_domain
    assert "Choose The Native Surface By State Space" in bayesian_domain
    assert "Convergence Is Quantity-Specific" in bayesian_domain
    assert "`evaluate_bayesian_posterior` evaluates" in bayesian_domain
    assert "not a convergence verdict" in bayesian_domain
    assert "Native Simulation Families" in simulation_domain
    assert "Public Simulation And Truth Surfaces" in simulation_domain
    assert "Keep Generative Truth Separate From Estimated Results" in simulation_domain
    assert "Preventing Simulation Circularity" in simulation_domain
    assert "Benchmark Design" in assurance_domain
    assert "Route A Claim To The Required Proof" in assurance_domain
    assert "Audit A Report From The Headline Backward" in assurance_domain
    assert "From Benchmark To Evidence" in assurance_domain
    assert "Tree Inspection, Validation, And Forensics" in diagnostics_domain
    assert "Diagnostic Consequence Matrix" in diagnostics_domain
    assert "Propagate Findings To Dependent Results" in diagnostics_domain
    assert "diagnose_tip_date_randomization" in diagnostics_domain
    assert "Conclusion Stability" in diagnostics_domain
    assert "Domain Handoffs" in domain_index
    assert "Cross-Domain Handoff Record" in domain_index
    assert "Attribute Failure To The Deciding Domain" in domain_index
    assert "Intake Contract" in data_index
    assert "Refuse Ambiguous Inputs" in data_index
    assert "The Identity Spine" in data_index
    assert "Handoff Requirements By Consumer" in data_index
    assert "Operations That Create A New Object" in data_index
    assert "Select The Inferential Contract" in inference_index
    assert "Review The Search, Not Only The Winner" in inference_index
    assert "Method Ownership Map" in inference_index
    assert "Distance-based inference" in domain_index
    assert "Bayesian inference" in domain_index
    assert "Interpret The Run State" in inference_index
    assert "Separate Objective, Exploration, And Support" in inference_index
    assert "Native And External Inference Paths" in inference_index
    assert "Search Evidence By Family" in inference_index
    assert "Completion is not a scientific verdict" in inference_index
    assert "Match The Question To The Model Family" in evolution_index
    assert "From Reconstructed State To Historical Claim" in evolution_index
    assert "Sensitivity Axes" in evolution_index
    assert "Establish The Observation Population" in evolution_index
    assert "Keep Three Result Layers Separate" in evolution_index
    assert "Refusal And Qualification Boundaries" in evolution_index
    assert "Choose Evidence That Matches The Claim" in assurance_index
    assert "Assurance Is A Decision Chain" in assurance_index
    assert "Preserve Negative Evidence" in assurance_index
    assert "Assurance Requires Independence And Coverage" in assurance_index
    assert "Denominator Contract" in assurance_index
    assert "Read Sensitivity As A Result" in assurance_index
    assert "Ecology and phylogeography" in domain_index
    assert "Diagnostics and sensitivity" in domain_index
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
    assert "Choose The Relationship From The Question" in parity_index
    assert "selected cases" in parity_index
    assert "Matched · mismatched" not in parity_index


def test_r_parity_guides_publish_registry_scope_and_non_claims() -> None:
    ecosystem = _read(PARITIES_ROOT / "r-ecosystem" / "index.md")
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
    assert "Read A Harness Result" in ecosystem
    assert "selected, passed, failed, and skipped" in ecosystem
    assert "Interpret The Case Status" in ape
    assert "case status and mismatch reason" in ape
    assert "Denominator And Model Identity" in geiger
    assert "They are not Evidence Book verdict names" in geiger
    assert "What A Passing Stochastic Case Means" in phytools
    assert "within-case replication denominators" in phytools
    assert "Registry Lifecycle" in ecosystem
    assert "Execution Identity" in ecosystem
    assert "Programmatic Entry Points" in ecosystem
    assert "Choose The Reference From The Operation" in ecosystem
    assert "Correspondence Key" in ecosystem
    assert "Programmatic Surface" in ape
    assert "Accept A Structural Comparison" in ape
    assert "Representative Function Correspondence" in ape
    assert "Read One Observation Row" in ape
    assert "Programmatic Surface" in geiger
    assert "Model Parameter Record" in geiger
    assert "Model Correspondence Key" in geiger
    assert "Diagnose The Optimizer Before The Parameter" in geiger
    assert "Programmatic Surface" in phytools
    assert "Preserve Both Denominators" in phytools
    assert "Function Correspondence Key" in phytools
    assert "The Stochastic Observation Unit" in phytools


def test_trust_guides_preserve_status_and_claim_ownership() -> None:
    trust_root = DOCS_ROOT / "05-bijux-phylogenetics-maintain" / "trust-surfaces"
    trust_index = _read(trust_root / "index.md")
    cross_surface = _read(trust_root / "public-trust-surface-maintenance.md")
    claim_maintenance = _read(trust_root / "trust-lane-and-study-guide-maintenance.md")

    assert "Questions Owned By Each Surface" in trust_index
    assert "Trust Decision Packet" in trust_index
    assert "Trace A Claim Across Surfaces" in trust_index
    assert "Failure Ownership Flow" in trust_index
    assert "Status Ownership" in cross_surface
    assert "Release Language Decision" in cross_surface
    assert "Keep Status Namespaces Separate" in claim_maintenance
    assert "Claim Dependency Record" in claim_maintenance
    assert "Promote A Result Into A Public Claim" in claim_maintenance
    assert "Reader-Facing Claim Anatomy" in claim_maintenance
    assert "Reconcile Conflicting Evidence Projections" in claim_maintenance
    assert "Handle Non-Comparable Claims" in claim_maintenance
    assert "contract-only wrapper" in claim_maintenance
    assert "intentionally open scientific verdict" not in cross_surface
    assert "public trust story" not in claim_maintenance


def test_engine_correspondence_guides_preserve_computational_ownership() -> None:
    engine_root = PARITIES_ROOT / "external-engines"
    index = _read(engine_root / "index.md")
    alignment = _read(engine_root / "alignment-and-trimming.md")
    inference = _read(engine_root / "tree-inference.md")
    bayesian = _read(engine_root / "bayesian-inference.md")

    assert "Three Kinds Of Record" in index
    assert "prepared analysis" in index
    assert "parsed artifact" in index
    assert "live workflow run" in index
    assert "Know What Was Actually Run" in index
    assert "Claim Ladder" in index
    assert "Preserve The Artifact Lineage" in index
    assert "an explicit parent identity" in index
    assert "MAFFT owns" in alignment
    assert "numerical parity claim" in alignment
    assert "Stage Acceptance And Handoff" in alignment
    assert "`--resume` is reuse under an identity check" in alignment
    assert "Codon-Aware Alignment Preserves The Translation Contract" in alignment
    assert "`gappyout`, `strict`, `strictplus`" in alignment
    assert "Preserve Sequence And Site Coordinates" in alignment
    assert "MAFFT column → trimAl result" in alignment
    assert "incomplete bundles" in inference
    assert "remain different statistical" in inference
    assert "quantities and are labeled" in inference
    assert "Engine Result Identity" in inference
    assert "Comparison Denominator" in inference
    assert "Choose The Command By Question" in inference
    assert "Robinson–Foulds distance" in inference
    assert "Reproducibility Has An Explicit Denominator" in inference
    assert "Retain Selection And Search Context" in inference
    assert "withdraw the paired conclusion" in inference
    assert "prepared-only analysis" in bayesian
    assert "does not dispatch through MrBayes or BEAST" in bayesian
    assert "Run Identity And Resume Safety" in bayesian
    assert "Resume therefore preserves execution lineage" in bayesian
    assert "Diagnostic Interpretation" in bayesian
    assert "Preparation, Execution, And Review Are Separate" in bayesian
    assert "Evidence Package Acceptance" in bayesian
    assert "Preserve Posterior Sample Lineage" in bayesian
    assert "resumed segments are not miscounted" in bayesian


def test_correspondence_matrix_distinguishes_verified_and_open_surfaces() -> None:
    study_references = _read(PARITIES_ROOT / "r-ecosystem" / "study-references.md")
    matrix = _read(PARITIES_ROOT / "correspondence" / "verification-matrix.md")
    reproduction = _read(
        PARITIES_ROOT / "correspondence" / "artifacts-and-reproduction.md"
    )

    assert "no standalone live `nlme` registry" in study_references
    assert "tip count of 75" in study_references
    assert "plot-only" in study_references
    assert "Record Origin Matrix" in study_references
    assert "Promotion Into A Study Verdict" in study_references
    assert "no posterior parity claim" in matrix
    assert "External adapter" in matrix
    assert "Mode, Outcome, And Verdict" in matrix
    assert "Verification Has Independent Axes" in matrix
    assert "Revision Rule" in matrix
    assert "not_comparable` verdicts and open posterior follow-through" in matrix
    assert "--manifest artifacts/parity/checked-manifest.json" in reproduction
    assert "Reconstruct The Denominator" in reproduction
    assert "Verify The Artifact Set As One Record" in reproduction
    assert "Reproduction Completion" in reproduction
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
    assert "Evaluate One Result" in quality_index
    assert "Assemble A Trust Packet" in quality_index
    assert "Select A Lane From The Claim" in lanes
    assert "mismatch_unexplained" in lanes
    assert "Trace One Native Inference Claim" in lanes
    assert "Propagate Failure To The Right Boundary" in lanes
    assert "Governed Verdict Meanings" in honesty
    assert "Observation Status And Bundle Verdict" in honesty
    assert "Match Language To The Record" in honesty
    assert "Handle Contradictions Explicitly" in honesty
    assert "engine_contract" in test_strategy
    assert "engine_real" in test_strategy
    assert "scientific_validation" in test_strategy
    assert "make package-verify" in test_strategy
    assert "Match The Lane To The Question" in test_strategy
    assert "Build Numerical Claims From Independent Checks" in test_strategy
    assert "Test Failure And Refusal Paths" in test_strategy
    assert "Read A Test Report Without Losing Scope" in test_strategy
    assert "native-bug" in benchmark
    assert "Missing And Failed Cases" in benchmark
    assert "Interpret The Six Dimensions" in benchmark
    assert "Respect Measurement Boundaries" in benchmark
    assert "Understand Current Wrapper Scope" in benchmark
    assert "Minimum Publication Packet" in benchmark
    assert "make check-release-readiness" in release
    assert "Non-Comparable Claims And Open Follow-Through" in release
    assert "may not rewrite the bundle verdict as “open,”" in release
    assert "Define Candidate Identity Before The Gate" in release
    assert "Read The Publish Scorecard" in release
    assert "publish-readiness.json" in release
    assert "Open Scientific Verdicts" not in release
    assert "Required Evidence Set" in review
    assert "Review Procedure" in review
    assert "Release Blockers" in review
    assert "Generate Reviewable Reports First" in review
    assert "Invalidate Stale Evidence" in review
    assert "Complete The Decision Packet" in review


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
    pcm1_index = _read(pcm1 / "index.md")
    preparation = _read(pcm1 / "data-and-tree-preparation.md")
    signal = _read(pcm1 / "signal-and-ancestral-results.md")
    models = _read(pcm2 / "model-comparisons-and-pgls.md")
    boundaries = _read(pcm2 / "signal-diagnostics-and-open-boundaries.md")

    assert "tree–trait join" in preparation
    assert "evidence-009" in preparation
    assert "All eight focused preparation bundles" in preparation
    assert (
        "aggregate `evidence-001` manifest is `matched_with_tolerance`" in preparation
    )
    assert "76 raw rows - 1 duplicate-species contraction" in preparation
    assert "Lagothrix_lagotricha" in preparation
    assert "Eulemur_fulvus" in preparation
    assert "Accept The Downstream Handoff" in preparation
    assert "matched_with_tolerance" in signal
    assert "not_comparable" in signal
    assert "31 scalar observations" in signal
    assert "The five seeded random datasets are parity fixtures" in signal
    assert "74 shared internal nodes" in signal
    assert "Bonobo–Gibbon MRCA Spot Check" in signal
    assert "Branch-Wise Increase Counts" in signal
    assert "What The Evidence Supports" in signal
    assert "Governed Observations" in pcm1_index
    assert "Analytical Population Dependency Graph" in pcm1_index
    assert "Citation Packet For PCM1" in pcm1_index
    assert "| estimated Pagel λ | 0.803 |" in pcm1_index
    assert "41.030408" in pcm1_index
    assert "520.078423" in pcm1_index
    assert "16 `matched` and 15" in pcm1_index
    assert "Pagel-λ PGLS" in models
    assert "likelihood-ratio" in models
    assert "Model Evidence Ledger" in models
    assert (
        "`evidence-008` | pairwise likelihood-ratio comparisons | `matched_with_tolerance`"
        in models
    )
    assert "56-row scalar ledger" in models
    assert "Focused Bundle Verdicts Do Not Erase Scalar Mismatches" in models
    assert "10.6718757865625" in models
    assert "15.170304341821" in models
    assert "direction or sign convention" in models
    assert "intercept-mode likelihood sweep" in boundaries
    assert "boundary register" in boundaries
    assert "0.095590314522858" in boundaries
    assert "Residual Diagnostics Use Named Scales" in boundaries
    assert "74 internal nodes" in boundaries
    pcm2_index = _read(pcm2 / "index.md")
    assert "Bundle And Scalar Verdicts" in pcm2_index
    assert "12 mismatch_unexplained" in pcm2_index
    assert "`evidence-009` | ancestral modes | `matched_with_tolerance`" in pcm2_index
    assert "Governed Findings" in pcm2_index
    assert "Unexplained Observation Inventory" in pcm2_index
    assert "Comparison Dependency Graph" in pcm2_index
    assert "Cite PCM2 Without Flattening It" in pcm2_index
    assert "10.6718757865625" in pcm2_index
    assert "15.170304341821" in pcm2_index


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
    assert "Observation-Level Mismatch Disclosure" in evidence_index
    assert "12 mismatch_unexplained" in evidence_index
    assert "Non-Comparable Evidence Is First-Class" in evidence_index
    assert "What The Evidence Supports At This Revision" in evidence_index
    assert "Cite The Record, Not This Summary" in evidence_index
    assert "Find The Record From The Statement" in evidence_index
    assert "Resolve A Value Through The Record Graph" in evidence_index
    assert "Read Coefficients From The Design Matrix" in pcm3
    assert "All eight bundle manifests record `not_comparable`" in pcm3
    assert "Verdict And Resolution State" in pcm3
    assert "No Numerical Headline Is Supported" in pcm3
    assert "Claim Closure Contract" in pcm3
    assert "parity verdicts remain open" not in pcm3
    pcm3_design = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm3-primate-continuous-and-categorical-regression"
        / "model-design-and-phylogenetic-regression.md"
    )
    assert "Governed Claim Ledger" in pcm3_design
    assert "Comparable Execution Requirements" in pcm3_design
    assert "every current bundle verdict\nis `not_comparable`" in pcm3_design
    assert "The Source Model Sequence" in pcm3_design
    assert "empty `primary_output_paths` list" in pcm3_design
    assert "Marginal And Sequential Tests Are Different Questions" in pcm3_design
    assert "What Can Be Cited Now" in pcm3_design
    assert "What Readers Should Watch Closely" not in pcm3_design
    assert "Model Components" in pcm4
    assert "The Preserved Source Analysis" in pcm4
    assert "m1a`, `m1b`, `m1c" in pcm4
    assert "Minimum Numerical Denominator" in pcm4
    assert "All seven bundle manifests record `not_comparable`" in pcm4
    assert "Verdict And Follow-Through" in pcm4
    assert "What Would Make The Posterior Comparable" in pcm4
    assert "Retained Record State" in pcm4
    assert "governed_primary_output_count: 0" in pcm4
    assert "bundle_contract_only" in pcm4
    assert "posterior parity verdicts\nremain open" not in pcm4
    assert "Preserve The Scale Of Every Quantity" in pcm5
    assert "The Preserved Source Model Sequence" in pcm5
    assert "m10PhyloSig" in pcm5
    assert "Transform Draws Before Summarizing" in pcm5
    assert "All nine bundle manifests record `not_comparable`" in pcm5
    assert "Comparison Contract By Model Family" in pcm5
    assert "Do not use open follow-through as the bundle verdict" in pcm5
    assert "Retained Record State" in pcm5
    assert "governed_primary_output_count: 0" in pcm5
    assert "bundle_contract_only" in pcm5
    assert "posterior parity verdicts remain open" not in pcm5
    pcm4_boundary = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models"
        / "mixed-model-evidence-boundary.md"
    )
    assert "Claim Closure Matrix" in pcm4_boundary
    assert "Closure Dependency Order" in pcm4_boundary
    assert "Distinguish Structural And Numerical Checks" in pcm4_boundary
    assert "Required Observation Ledger" in pcm4_boundary
    assert (
        "No PCM4 bundle currently retains a governed primary scientific output"
        in pcm4_boundary
    )
    pcm5_boundary = _read(
        EVIDENCE_GUIDE_ROOT
        / "pcm5-cichlid-phylogenetic-generalized-linear-mixed-models"
        / "generalized-mixed-model-evidence-boundary.md"
    )
    assert "Claim-Specific Closure Requirements" in pcm5_boundary
    assert "Closure Dependency Graph" in pcm5_boundary
    assert "Derived-Quantity Ledger" in pcm5_boundary
    assert "What Execution Would Change The Verdict State" in pcm5_boundary
    assert "Observation And Diagnostic Ledger" in pcm5_boundary
    assert "zero governed primary\n  outputs" in pcm5_boundary
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
    assert "Accept Or Refuse A Result" in runtime
    assert "Confirm Package Identity" in runtime
    assert "What Installation Does And Does Not Supply" in runtime
    assert "Publish A Reconstructable Result" in runtime
    assert "Trace The Installed Result Lineage" in runtime
    assert "Verify The Alias End To End" in alias
    assert "Diagnose A Difference" in alias
    assert "Resolve Both Distribution Identities" in alias
    assert "Compatibility Failure Classes" in alias
    assert "two scientific owners" in alias
    assert "It does not create a second implementation" in alias
    assert "Prove Delegation, Not Similarity" in alias
    assert "What The Governed Commands Establish" in maintainer
    assert "Change-To-Check Map" in maintainer
    assert "make validate-evidence-book" in maintainer
    assert "Governed Report Map" in maintainer
    assert "Read Report Dependencies" in maintainer
    assert "Handle A Failed Gate Without Bypassing It" in maintainer
    assert "Verification Handoff Record" in maintainer
    assert "Gate Failure Ownership" in maintainer
    assert "blocked follow-through belongs to the resolution workflow" in maintainer


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
    assert "Scale Review By Semantic Risk" in handbook
    assert "Resolve Contradictions Before Completion" in handbook
    assert "Maintained Invariants" in handbook
    assert "Keep Status In Its Owning Record" in handbook
    assert "Resolve Canonical Ownership" in handbook
    assert "Change-To-Check Routing" in handbook
    assert (
        "ready release can still publish explicit `not_comparable` boundaries"
        in handbook
    )
    assert "Canonicality Rules" in ownership
    assert "Change Impact Review" in contracts
    assert "Contract Change Classes" in contracts
    assert "Completeness Audit" in contracts
    assert "Classify Before Writing" in artifacts
    assert "Artifact Lifecycle Record" in artifacts
    assert "make check-artifact-governance" in artifacts
    assert "Safe Cleanup" in artifacts
    assert "Alignment Matrix" in release
    assert "Release Decision Record" in release
    assert "make check-release-readiness" in release
    assert "Release Blockers" in release
    assert "Limitation Versus Contradiction" in release
    assert "Build The Alignment Packet" in release
    assert "Compare Meaning Across Surfaces" in release
    assert "There is no single" in release
    assert "`release-readiness.json`" in release
    assert (
        "does not mean every scientific comparison has a favorable verdict" in release
    )
    assert "narrative alignment" not in release


def test_maintainer_indexes_define_operational_decisions_without_meta_prose() -> None:
    documentation = _read(MAINTAINER_ROOT / "documentation-system" / "index.md")
    trust = _read(MAINTAINER_ROOT / "trust-surfaces" / "index.md")
    public_interface = _read(
        MAINTAINER_ROOT
        / "trust-surfaces"
        / "public-interface-and-artifact-maintenance.md"
    )
    operations = _read(MAINTAINER_ROOT / "repository-operations" / "index.md")

    assert "Resolve Ownership Before Editing" in documentation
    assert "Separate Reader And Maintainer Language" in documentation
    assert "Propagate A Contract Change" in documentation
    assert "Documentation Review Packet" in documentation
    assert "Completion Criteria" in documentation
    assert "Detect Cross-Surface Drift" in trust
    assert "Release Review Questions" in trust
    assert "Contract Projection Matrix" in public_interface
    assert "Maintain The Artifact Handoff" in public_interface
    assert "Public Result Invariants" in public_interface
    assert "Review A Material Interface Change" in public_interface
    assert "Publication Decision" in public_interface
    assert "Classify Every Output" in operations
    assert "Follow The Repository State Machine" in operations
    assert "Decide Whether Generated Output Is Governed" in operations
    assert "Inspect Each Handoff" in operations
    assert "Repository Gate Routing" in operations
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
