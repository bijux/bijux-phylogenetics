from __future__ import annotations

from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_ou_traits,
    simulate_protein_alignment,
)
from bijux_phylogenetics.validation.validation_corpus.benchmark_dashboards import (
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_method_accuracy_dashboard,
    build_runtime_benchmark_dashboard,
    build_workflow_practical_limit_dashboard,
)
from bijux_phylogenetics.validation.validation_corpus.dataset_corpora import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
    default_fixtures_root as _default_fixtures_root,
    fixture as _fixture,
)
from bijux_phylogenetics.validation.validation_corpus.regression_corpus import (
    build_regression_dataset_corpus,
)
from bijux_phylogenetics.validation.validation_corpus.contracts import (
    BenchmarkCorpusReport,
    BenchmarkDashboardRow,
    CorpusDatasetCase,
    CorpusDatasetCaseResult,
    LargeAlignmentScalingBenchmarkDashboard,
    LargeTreeScalingBenchmarkDashboard,
    LargeTreeSetScalingBenchmarkDashboard,
    MemoryBenchmarkDashboard,
    MethodAccuracyDashboard,
    MethodAccuracyRow,
    MethodLimitationEntry,
    MethodLimitationRegistry,
    RegressionDatasetCaseResult,
    RegressionDatasetCorpusReport,
    RuntimeBenchmarkDashboard,
    ScientificValidationClaim,
    ScientificValidationReport,
    SimulationReproducibilityCase,
    SimulationReproducibilityReport,
    WorkflowPracticalLimitDashboard,
)

_CORPUS_EXPORT_TYPES = (
    BenchmarkCorpusReport,
    BenchmarkDashboardRow,
    CorpusDatasetCase,
    CorpusDatasetCaseResult,
    LargeAlignmentScalingBenchmarkDashboard,
    LargeTreeScalingBenchmarkDashboard,
    LargeTreeSetScalingBenchmarkDashboard,
    MemoryBenchmarkDashboard,
    MethodAccuracyDashboard,
    MethodAccuracyRow,
    MethodLimitationEntry,
    MethodLimitationRegistry,
    RegressionDatasetCaseResult,
    RegressionDatasetCorpusReport,
    RuntimeBenchmarkDashboard,
    ScientificValidationClaim,
    ScientificValidationReport,
    SimulationReproducibilityCase,
    SimulationReproducibilityReport,
    WorkflowPracticalLimitDashboard,
)

_CORPUS_EXPORT_BUILDERS = (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_messy_benchmark_corpus,
    build_regression_dataset_corpus,
    build_method_accuracy_dashboard,
    build_runtime_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_workflow_practical_limit_dashboard,
)

def build_method_limitation_registry() -> MethodLimitationRegistry:
    """Enumerate major method families with explicit assumptions and trust boundaries."""
    return MethodLimitationRegistry(
        goal_id=250,
        entries=[
            MethodLimitationEntry(
                method="tree-validation",
                status="validated",
                validated_by=["level1-reference-validation", "broken-benchmark-corpus"],
                assumptions=[
                    "input tree syntax is parseable",
                    "tip labels are intended to identify taxa",
                ],
                invalid_inputs=[
                    "duplicate tip labels",
                    "negative branch lengths",
                    "malformed rootedness assumptions",
                ],
                limitations=[
                    "validation does not by itself prove biological interpretation or engine compatibility beyond audited surfaces"
                ],
            ),
            MethodLimitationEntry(
                method="dataset-audit",
                status="validated",
                validated_by=[
                    "clean-benchmark-corpus",
                    "messy-benchmark-corpus",
                    "regression-dataset-corpus",
                ],
                assumptions=[
                    "metadata and trait tables use stable taxon keys",
                    "caller-provided files correspond to one biological dataset",
                ],
                invalid_inputs=[
                    "missing tree taxa in metadata or traits",
                    "unsafe alignments for inference",
                    "invalid tip dates or calibrations",
                ],
                limitations=[
                    "dataset readiness is reviewer-facing triage, not a substitute for downstream method-specific validation"
                ],
            ),
            MethodLimitationEntry(
                method="distance-methods",
                status="validated",
                validated_by=[
                    "distance reference fixtures",
                    "runtime benchmark dashboard",
                ],
                assumptions=[
                    "aligned homologous sequences",
                    "distance model matches alphabet and saturation regime",
                ],
                invalid_inputs=[
                    "too few comparable sites",
                    "severe saturation",
                    "unsafe ambiguity handling for the chosen policy",
                ],
                limitations=[
                    "distance trees remain approximations and can disagree with likelihood or Bayesian inference"
                ],
            ),
            MethodLimitationEntry(
                method="comparative-models",
                status="validated",
                validated_by=[
                    "comparative reference fixtures",
                    "level1-reference-validation",
                ],
                assumptions=[
                    "tree/taxon linkage is correct",
                    "traits satisfy declared model assumptions",
                ],
                invalid_inputs=[
                    "overfit models",
                    "non-identifiable OU settings",
                    "missing trait coverage after pruning",
                ],
                limitations=[
                    "validated examples do not remove the need for biological judgment about causality or model adequacy"
                ],
            ),
            MethodLimitationEntry(
                method="ancestral-reconstruction",
                status="experimental",
                validated_by=[
                    "simulation validation surfaces",
                    "ancestral reference examples",
                ],
                assumptions=[
                    "chosen transition or continuous model is defensible",
                    "tree uncertainty has been considered",
                ],
                invalid_inputs=[
                    "impossible discrete coding",
                    "low-information or unstable internal nodes",
                    "unsafe pruning sensitivity",
                ],
                limitations=[
                    "external tool comparison and maturity gates are still incomplete for every reconstruction mode"
                ],
            ),
            MethodLimitationEntry(
                method="bayesian-time-tree",
                status="experimental",
                validated_by=[
                    "tip-date validation",
                    "calibration validation",
                    "BEAST and MrBayes workflow surfaces",
                ],
                assumptions=[
                    "valid dates, calibrations, priors, and convergence diagnostics",
                    "posterior sampling has mixed adequately",
                ],
                invalid_inputs=[
                    "invalid tip dates",
                    "invalid calibrations",
                    "low ESS or conflicting independent runs",
                ],
                limitations=[
                    "workflow support exists, but cross-environment reproducibility and broader benchmark validation remain incomplete"
                ],
            ),
        ],
        limitations=[
            "registry statuses summarize the current checked-in evidence and should move only when new validation surfaces are added or removed",
        ],
    )


def build_scientific_validation_report(
    *, fixtures_root: Path | None = None
) -> ScientificValidationReport:
    """Separate validated, unvalidated, experimental, and unsafe claims for reviewers."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    accuracy = build_method_accuracy_dashboard(fixtures_root=root)
    registry = build_method_limitation_registry()
    claims = [
        ScientificValidationClaim(
            status="validated",
            claim="checked-in Level 1 reference fixtures and benchmark corpora currently pass as expected",
            evidence=[
                f"{row.surface}:{row.passed_count}/{row.coverage_count}"
                for row in accuracy.rows
            ],
        ),
        ScientificValidationClaim(
            status="validated",
            claim="clean, broken, messy, and regression dataset corpora all preserve their expected trust signals",
            evidence=[
                "clean-benchmark-corpus",
                "broken-benchmark-corpus",
                "messy-benchmark-corpus",
                "regression-dataset-corpus",
            ],
        ),
        ScientificValidationClaim(
            status="experimental",
            claim="ancestral reconstruction and Bayesian time-tree workflows have substantive support but still carry incomplete maturity evidence",
            evidence=[
                entry.method
                for entry in registry.entries
                if entry.status == "experimental"
            ],
        ),
        ScientificValidationClaim(
            status="unvalidated",
            claim="cross-environment reproducibility is not yet claimed, and external R ecosystem comparisons remain narrower than the broader validation roadmap",
            evidence=[
                "goal-252 remains outside this iteration",
                "evidence-book/studies/primate-longevity-signal/evidence-001 demonstrates one checked-in R comparison bundle",
            ],
        ),
        ScientificValidationClaim(
            status="unsafe",
            claim="publication-grade time-tree conclusions remain unsafe when tip dates, calibrations, or convergence diagnostics fail",
            evidence=[
                "dataset-audit blockers include invalid tip dates and invalid calibrations",
                "bayesian-time-tree is still marked experimental in the limitation registry",
            ],
        ),
    ]
    return ScientificValidationReport(
        goal_id=249,
        claims=claims,
        limitations=[
            "the report summarizes current checked-in evidence; it does not replace method-specific diagnostics on a new biological dataset",
        ],
    )


def _normalize_jsonable(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: _normalize_jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _normalize_jsonable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_jsonable(item) for item in value]
    return value


def write_validation_corpus_json(path: Path, report: object) -> Path:
    """Write a validation-corpus or dashboard report as deterministic JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _normalize_jsonable(report)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def validate_simulation_reproducibility(
    *, fixtures_root: Path | None = None
) -> SimulationReproducibilityReport:
    """Verify that repeated simulations with the same seed produce identical structured results."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    tree_path = _fixture(root, "trees", "example_tree.nwk")

    def digest(payload: object) -> tuple[str, str]:
        normalized = _normalize_jsonable(payload)
        encoded = json.dumps(normalized, sort_keys=True).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest(), encoded.decode("utf-8")

    cases: list[tuple[str, object]] = [
        (
            "birth-death-tree-set",
            lambda: (
                simulate_birth_death_trees(tree_count=3, tip_count=4, seed=7)[1].records
            ),
        ),
        (
            "brownian-traits",
            lambda: simulate_brownian_traits(
                tree_path, root_state=1.0, sigma=0.5, seed=7
            ),
        ),
        (
            "ou-traits",
            lambda: simulate_ou_traits(
                tree_path, root_state=1.0, sigma=0.5, alpha=0.7, theta=0.2, seed=7
            ),
        ),
        (
            "discrete-traits",
            lambda: simulate_discrete_traits(
                tree_path, states=["north", "south"], transition_rate=0.8, seed=7
            ),
        ),
        (
            "dna-alignment",
            lambda: simulate_dna_alignment(
                tree_path, sequence_length=16, substitution_rate=0.9, seed=7
            ),
        ),
        (
            "protein-alignment",
            lambda: simulate_protein_alignment(
                tree_path, sequence_length=12, substitution_rate=0.6, seed=7
            ),
        ),
    ]

    results: list[SimulationReproducibilityCase] = []
    for surface, callback in cases:
        first = callback()
        second = callback()
        first_digest, first_payload = digest(first)
        second_digest, second_payload = digest(second)
        notes: list[str] = []
        if first_digest != second_digest or first_payload != second_payload:
            notes.append("same-seed simulation output drifted between repeated runs")
        results.append(
            SimulationReproducibilityCase(
                surface=surface,
                passed=not notes,
                digest=first_digest,
                notes=notes,
            )
        )
    return SimulationReproducibilityReport(
        goal_id=251,
        passed=all(case.passed for case in results),
        cases=results,
        limitations=[
            "the reproducibility check covers deterministic seeded library surfaces and does not yet assert cross-environment bit-for-bit stability",
        ],
    )
