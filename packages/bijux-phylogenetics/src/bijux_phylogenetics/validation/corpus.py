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
from bijux_phylogenetics.validation.validation_corpus.scientific_review import (
    build_method_limitation_registry as _build_method_limitation_registry,
    build_scientific_validation_report as _build_scientific_validation_report,
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
    return _build_method_limitation_registry()


def build_scientific_validation_report(
    *, fixtures_root: Path | None = None
) -> ScientificValidationReport:
    """Separate validated, unvalidated, experimental, and unsafe claims for reviewers."""
    return _build_scientific_validation_report(fixtures_root=fixtures_root)


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
