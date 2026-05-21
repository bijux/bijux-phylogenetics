from __future__ import annotations

from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.benchmark import (
    benchmark_alignment_site_scaling,
    benchmark_large_alignment_scaling,
    benchmark_large_tree_scaling,
    benchmark_large_tree_set_scaling,
    benchmark_tree_comparison,
    benchmark_tree_set_consensus,
    benchmark_tree_validation,
    benchmark_workflow_practical_limits,
)
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_brownian_traits,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_ou_traits,
    simulate_protein_alignment,
)
from bijux_phylogenetics.validation.reference import (
    build_core_workflow_validation_report,
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


def build_method_accuracy_dashboard(
    *, fixtures_root: Path | None = None
) -> MethodAccuracyDashboard:
    """Summarize validation accuracy, error counts, and coverage across benchmark surfaces."""
    root = _default_fixtures_root() if fixtures_root is None else fixtures_root
    core = build_core_workflow_validation_report(fixtures_root=root)
    clean = build_clean_benchmark_corpus(fixtures_root=root)
    broken = build_broken_benchmark_corpus(fixtures_root=root)
    messy = build_messy_benchmark_corpus(fixtures_root=root)
    regression = build_regression_dataset_corpus(fixtures_root=root)

    def row(
        surface: str,
        passed_count: int,
        failed_count: int,
        coverage_count: int,
        bias_notes: list[str],
        error_notes: list[str],
    ) -> MethodAccuracyRow:
        total = max(coverage_count, 1)
        return MethodAccuracyRow(
            surface=surface,
            accuracy=round(passed_count / total, 15),
            passed_count=passed_count,
            failed_count=failed_count,
            coverage_count=coverage_count,
            bias_notes=bias_notes,
            error_notes=error_notes,
        )

    rows = [
        row(
            "level1-reference-validation",
            core.passed_fixture_count,
            core.failed_fixture_count,
            core.total_fixture_count,
            core.limitations,
            [case.fixture_name for case in core.failure_gallery if not case.passed],
        ),
        row(
            "clean-benchmark-corpus",
            clean.passed_case_count,
            clean.failed_case_count,
            clean.case_count,
            clean.limitations,
            [case.name for case in clean.cases if not case.passed],
        ),
        row(
            "broken-benchmark-corpus",
            broken.passed_case_count,
            broken.failed_case_count,
            broken.case_count,
            broken.limitations,
            [case.name for case in broken.cases if not case.passed],
        ),
        row(
            "messy-benchmark-corpus",
            messy.passed_case_count,
            messy.failed_case_count,
            messy.case_count,
            messy.limitations,
            [case.name for case in messy.cases if not case.passed],
        ),
        row(
            "regression-dataset-corpus",
            regression.passed_case_count,
            regression.failed_case_count,
            regression.case_count,
            regression.limitations,
            [case.name for case in regression.cases if not case.passed],
        ),
    ]
    return MethodAccuracyDashboard(
        goal_id=246,
        rows=rows,
        limitations=[
            "accuracy currently summarizes checked-in fixture and corpus pass rates; it does not yet replace external software comparison studies",
        ],
    )


def build_runtime_benchmark_dashboard(
    *, replicates: int = 1
) -> RuntimeBenchmarkDashboard:
    """Summarize runtime scaling across taxa, sites, tree counts, and posterior-like samples."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(
                replicates=replicates
            ).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(
                replicates=replicates
            ).observations,
        ),
    ]
    return RuntimeBenchmarkDashboard(
        goal_id=247,
        rows=rows,
        limitations=[
            "runtime summaries measure local benchmark fixtures and should be re-run on target hardware before operational promises are made",
        ],
    )


def build_memory_benchmark_dashboard(
    *, replicates: int = 1
) -> MemoryBenchmarkDashboard:
    """Summarize peak memory scaling across the main benchmark axes."""
    rows = [
        BenchmarkDashboardRow(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=benchmark_tree_validation(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=benchmark_tree_comparison(replicates=replicates).observations,
        ),
        BenchmarkDashboardRow(
            workflow="alignment-diagnostics",
            scaling_axis="sites",
            observations=benchmark_alignment_site_scaling(
                replicates=replicates
            ).observations,
        ),
        BenchmarkDashboardRow(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=benchmark_tree_set_consensus(
                replicates=replicates
            ).observations,
        ),
    ]
    return MemoryBenchmarkDashboard(
        goal_id=248,
        rows=rows,
        limitations=[
            "memory summaries capture Python-side peak allocations during benchmark runs and do not model every external-engine workflow",
        ],
    )


def build_large_tree_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    tip_counts: list[int] | None = None,
) -> LargeTreeScalingBenchmarkDashboard:
    """Summarize realistic large-tree scaling across validation and review workflows."""
    report = benchmark_large_tree_scaling(
        replicates=replicates,
        tip_counts=tip_counts,
    )
    return LargeTreeScalingBenchmarkDashboard(
        goal_id=221,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_large_alignment_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeAlignmentScalingBenchmarkDashboard:
    """Summarize realistic large-alignment scaling across review workflows."""
    report = benchmark_large_alignment_scaling(
        replicates=replicates,
        size_classes=size_classes,
    )
    return LargeAlignmentScalingBenchmarkDashboard(
        goal_id=222,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_large_tree_set_scaling_benchmark_dashboard(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeTreeSetScalingBenchmarkDashboard:
    """Summarize realistic large-tree-set scaling across posterior review workflows."""
    report = benchmark_large_tree_set_scaling(
        replicates=replicates,
        size_classes=size_classes,
    )
    return LargeTreeSetScalingBenchmarkDashboard(
        goal_id=223,
        workflows=report.workflows,
        limitations=report.limitations,
    )


def build_workflow_practical_limit_dashboard(
    *,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> WorkflowPracticalLimitDashboard:
    """Summarize the largest governed workflow classes currently exercised in benchmark and stress lanes."""
    report = benchmark_workflow_practical_limits(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    return WorkflowPracticalLimitDashboard(
        goal_id=224,
        entries=report.entries,
        limitations=report.limitations,
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
