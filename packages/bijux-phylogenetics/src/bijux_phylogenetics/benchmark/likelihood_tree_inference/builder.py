from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile
import time
import tracemalloc

from bijux_phylogenetics.benchmark._fixtures import (
    NATIVE_MAXIMUM_LIKELIHOOD_SPEED_SIZE_CLASSES,
    write_large_nucleotide_alignment,
)
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
    infer_nucleotide_likelihood_tree,
    infer_nucleotide_likelihood_tree_from_alignment,
    optimize_fixed_topology_nucleotide_branch_lengths,
    validate_nucleotide_likelihood_tree_inference_search_method,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_bipartition,
    canonical_clade_id,
    informative_unrooted_splits,
    robinson_foulds_metrics,
    split_sort_key,
)

from .contracts import (
    MaximumLikelihoodWrapperBenchmarkBundle,
    MaximumLikelihoodWrapperBenchmarkObservation,
    MaximumLikelihoodWrapperBenchmarkReport,
    MaximumLikelihoodWrapperBenchmarkSummaryRow,
    MaximumLikelihoodWrapperSupportComparison,
    NativeMaximumLikelihoodAccuracyBenchmarkBundle,
    NativeMaximumLikelihoodAccuracyBenchmarkObservation,
    NativeMaximumLikelihoodAccuracyBenchmarkReport,
    NativeMaximumLikelihoodMemoryBenchmarkBundle,
    NativeMaximumLikelihoodMemoryBenchmarkObservation,
    NativeMaximumLikelihoodMemoryBenchmarkReport,
    NativeMaximumLikelihoodMemoryTrendRow,
    NativeMaximumLikelihoodSuiteObservation,
    NativeMaximumLikelihoodSuiteReport,
    NativeMaximumLikelihoodSuiteSummaryRow,
    NativeMaximumLikelihoodSpeedBenchmarkBundle,
    NativeMaximumLikelihoodSpeedBenchmarkObservation,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
    NativeMaximumLikelihoodSpeedTrendRow,
)
from .accuracy_fixtures import (
    build_native_maximum_likelihood_accuracy_truth_tree,
    governed_native_maximum_likelihood_accuracy_cases,
    simulate_native_maximum_likelihood_accuracy_alignment,
)
from .wrapper_references import MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS
from .artifact_outputs import (
    write_maximum_likelihood_wrapper_comparison_observation_table,
    write_maximum_likelihood_wrapper_comparison_run_json,
    write_maximum_likelihood_wrapper_comparison_summary_table,
    write_native_maximum_likelihood_accuracy_observation_table,
    write_native_maximum_likelihood_accuracy_run_json,
    write_native_maximum_likelihood_accuracy_summary_table,
    write_native_maximum_likelihood_accuracy_wrapper_table,
    write_native_maximum_likelihood_memory_observation_table,
    write_native_maximum_likelihood_memory_run_json,
    write_native_maximum_likelihood_memory_summary_table,
    write_native_maximum_likelihood_memory_trend_table,
    write_native_maximum_likelihood_speed_observation_table,
    write_native_maximum_likelihood_speed_run_json,
    write_native_maximum_likelihood_speed_summary_table,
    write_native_maximum_likelihood_speed_trend_table,
)

_MAXIMUM_LIKELIHOOD_WRAPPER_BENCHMARK_STATUSES = (
    "tolerance-match",
    "expected-model-assumption-difference",
    "unsupported-case",
    "native-bug",
)
_LIKELIHOOD_TOLERANCE = 1e-3
_SUPPORT_TOLERANCE = 1e-9
_SUPPORTED_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES = (
    "jc69",
    "k80",
    "f81",
    "hky85",
)
_DEFAULT_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES = (
    "jc69",
    "k80",
)
_NATIVE_MAXIMUM_LIKELIHOOD_SUITE_STATUSES = (
    "native-advantage",
    "native-bug",
    "expected-assumption-difference",
    "unsupported-case",
)
_NATIVE_MAXIMUM_LIKELIHOOD_ACCURACY_SHORTFALL_TOLERANCE = 0.25
_ELEVATED_MEMORY_RATIO_THRESHOLD = 1.25
_ELEVATED_MEMORY_DELTA_THRESHOLD_BYTES = 100_000
_MEMORY_BLOWUP_RATIO_THRESHOLD = 2.0
_MEMORY_BLOWUP_DELTA_THRESHOLD_BYTES = 500_000


@dataclass(frozen=True, slots=True)
class _MeasuredNativeTreeInference:
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    log_likelihood: float
    tree_newick: str
    topology_signature: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _MeasuredNativeBootstrapSupport:
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    log_likelihood: float
    tree_newick: str
    topology_signature: tuple[str, ...]
    support_by_split_id: dict[str, float]


@dataclass(frozen=True, slots=True)
class _MeasuredNativeMaximumLikelihoodLane:
    size_class: str
    model_name: str
    search_method: str
    taxon_count: int
    site_count: int
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    selected_model_name: str
    best_run_source_label: str
    search_algorithm: str
    final_log_likelihood: float


def benchmark_native_maximum_likelihood_accuracy(
    *,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodAccuracyBenchmarkReport:
    """Benchmark native ML recovery against governed simulated truth cases."""
    if start_tree_count < 1:
        raise ValueError(f"start_tree_count must be at least 1, got {start_tree_count}")
    resolved_model_names = _resolve_native_maximum_likelihood_accuracy_model_names(
        model_names
    )
    resolved_search_methods = [
        validate_nucleotide_likelihood_tree_inference_search_method(search_method)
        for search_method in (search_methods or ("nni",))
    ]
    cases = [
        case
        for case in governed_native_maximum_likelihood_accuracy_cases()
        if case.model_name in resolved_model_names
    ]
    observations: list[NativeMaximumLikelihoodAccuracyBenchmarkObservation] = []
    for case in cases:
        records = simulate_native_maximum_likelihood_accuracy_alignment(case)
        true_topology_log_likelihood = _optimize_native_accuracy_truth_topology(
            case,
            records,
            upper_branch_length_bound=upper_branch_length_bound,
        )
        for search_method in resolved_search_methods:
            inference_report = infer_nucleotide_likelihood_tree(
                records,
                model_name=case.model_name,
                search_method=search_method,
                start_tree_count=start_tree_count,
                start_tree_seed=start_tree_seed,
                upper_branch_length_bound=upper_branch_length_bound,
            )
            observations.append(
                _build_native_maximum_likelihood_accuracy_observation(
                    case,
                    search_method=search_method,
                    inference_report=inference_report,
                    true_topology_log_likelihood=true_topology_log_likelihood,
                )
            )
    wrapper_report = benchmark_maximum_likelihood_wrapper_comparison()
    topology_recovery_case_count = sum(
        1 for row in observations if row.topology_recovered
    )
    parameter_recovery_case_count = sum(
        1 for row in observations if row.parameter_recovery_status == "parameter-recovered"
    )
    parameter_not_applicable_case_count = sum(
        1 for row in observations if row.parameter_recovery_status == "parameter-not-applicable"
    )
    likelihood_quality_case_count = sum(
        1
        for row in observations
        if row.likelihood_quality_status == "no-likelihood-shortfall"
    )
    wrapper_clear_case_count = sum(
        1 for row in wrapper_report.observations if row.supported and not row.blocking
    )
    all_truth_cases_clear = all(
        row.topology_recovered
        and row.likelihood_quality_status == "no-likelihood-shortfall"
        and row.parameter_recovery_status
        in {"parameter-recovered", "parameter-not-applicable"}
        for row in observations
    )
    return NativeMaximumLikelihoodAccuracyBenchmarkReport(
        observations=observations,
        wrapper_observations=wrapper_report.observations,
        model_names=resolved_model_names,
        search_methods=resolved_search_methods,
        topology_recovery_case_count=topology_recovery_case_count,
        parameter_recovery_case_count=parameter_recovery_case_count,
        parameter_not_applicable_case_count=parameter_not_applicable_case_count,
        likelihood_quality_case_count=likelihood_quality_case_count,
        wrapper_supported_case_count=wrapper_report.supported_case_count,
        wrapper_clear_case_count=wrapper_clear_case_count,
        all_truth_cases_clear=all_truth_cases_clear,
        all_supported_wrapper_cases_clear=wrapper_report.all_supported_cases_clear,
        limitations=[
            "governed accuracy fixtures use deterministic native simulation seeds so topology and parameter recovery can be compared to known truth instead of anecdotal ad hoc examples",
            "the default governed accuracy roster runs the K80 truth lane; the F81 composition-bias lane remains available through explicit model overrides when additional parameter-recovery coverage is needed",
            "parameter recovery is scored from the fixed-model selection surface that seeds tree search, because the current native tree-inference workflow keeps substitution parameters fixed during local topology search",
            "likelihood quality is measured as shortfall against the true simulated topology scored under the known generating substitution parameters with native branch-length optimization, so search shortfalls remain explicit without turning the benchmark into a second full fitting workflow",
            "wrapper evidence remains limited to the cached governed JC69 comparison lanes already owned in-repository, so external-engine coverage is narrower than the native simulated truth roster",
        ],
    )


def benchmark_native_maximum_likelihood_speed(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodSpeedBenchmarkReport:
    """Benchmark native ML runtime across governed model, size, and search lanes."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    if start_tree_count < 1:
        raise ValueError(f"start_tree_count must be at least 1, got {start_tree_count}")
    resolved_size_classes = _validate_native_maximum_likelihood_benchmark_size_classes(
        size_classes or list(NATIVE_MAXIMUM_LIKELIHOOD_SPEED_SIZE_CLASSES)
    )
    resolved_model_names = _validate_native_maximum_likelihood_benchmark_model_names(
        model_names or _DEFAULT_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES
    )
    resolved_search_methods = [
        validate_nucleotide_likelihood_tree_inference_search_method(search_method)
        for search_method in (
            search_methods
            or ("nni", "spr", "tbr")
        )
    ]

    measured_lanes = _measure_native_maximum_likelihood_benchmark_lanes(
        size_classes=resolved_size_classes,
        model_names=resolved_model_names,
        search_methods=resolved_search_methods,
        replicates=replicates,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
        temporary_directory_prefix="bijux-native-maximum-likelihood-speed-",
    )
    observations = [
        _build_native_maximum_likelihood_speed_observation(lane)
        for lane in measured_lanes
    ]

    size_class_labels = [label for label, _, _ in resolved_size_classes]
    trend_rows = _build_native_maximum_likelihood_speed_trend_rows(
        observations,
        size_class_labels=size_class_labels,
        model_names=resolved_model_names,
        search_methods=resolved_search_methods,
    )
    return NativeMaximumLikelihoodSpeedBenchmarkReport(
        replicates=replicates,
        observations=observations,
        trend_rows=trend_rows,
        size_class_labels=size_class_labels,
        model_names=list(resolved_model_names),
        search_methods=resolved_search_methods,
        limitations=[
            "native maximum-likelihood speed lanes use deterministic synthetic nucleotide alignments so runtime trends measure owned search cost without external engine variance",
            "replicates reuse the same governed start-tree seed so runtime comparisons stay lane-stable rather than mixing stochastic workload differences into the benchmark",
            "the default governed speed roster samples jc69 and k80; f81 and hky85 remain available through explicit model overrides when broader model timing is needed",
            "the governed native speed roster currently excludes fixed-parameter GTR because fixed-tree native inference still requires caller-supplied exchangeabilities for that lane",
            "peak memory is recorded with python tracemalloc as supporting context, but this benchmark treats runtime trend reporting as the primary completion signal",
        ],
    )


def benchmark_native_maximum_likelihood_memory(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodMemoryBenchmarkReport:
    """Benchmark native ML peak memory across governed model, size, and search lanes."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    if start_tree_count < 1:
        raise ValueError(f"start_tree_count must be at least 1, got {start_tree_count}")
    resolved_size_classes = _validate_native_maximum_likelihood_benchmark_size_classes(
        size_classes or list(NATIVE_MAXIMUM_LIKELIHOOD_SPEED_SIZE_CLASSES)
    )
    resolved_model_names = _validate_native_maximum_likelihood_benchmark_model_names(
        model_names or _DEFAULT_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES
    )
    resolved_search_methods = [
        validate_nucleotide_likelihood_tree_inference_search_method(search_method)
        for search_method in (
            search_methods
            or ("nni", "spr", "tbr")
        )
    ]

    measured_lanes = _measure_native_maximum_likelihood_benchmark_lanes(
        size_classes=resolved_size_classes,
        model_names=resolved_model_names,
        search_methods=resolved_search_methods,
        replicates=replicates,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
        temporary_directory_prefix="bijux-native-maximum-likelihood-memory-",
    )
    observations = [
        _build_native_maximum_likelihood_memory_observation(lane)
        for lane in measured_lanes
    ]
    size_class_labels = [label for label, _, _ in resolved_size_classes]
    trend_rows = _build_native_maximum_likelihood_memory_trend_rows(
        observations,
        size_class_labels=size_class_labels,
        model_names=resolved_model_names,
        search_methods=resolved_search_methods,
    )
    stable_memory_case_count = sum(
        1
        for row in trend_rows
        if row.memory_growth_status == "stable-memory-growth"
    )
    elevated_memory_case_count = sum(
        1
        for row in trend_rows
        if row.memory_growth_status == "elevated-memory-growth"
    )
    memory_blowup_case_count = sum(
        1 for row in trend_rows if row.memory_growth_status == "memory-blowup"
    )
    return NativeMaximumLikelihoodMemoryBenchmarkReport(
        replicates=replicates,
        observations=observations,
        trend_rows=trend_rows,
        size_class_labels=size_class_labels,
        model_names=list(resolved_model_names),
        search_methods=resolved_search_methods,
        stable_memory_case_count=stable_memory_case_count,
        elevated_memory_case_count=elevated_memory_case_count,
        memory_blowup_case_count=memory_blowup_case_count,
        limitations=[
            "native maximum-likelihood memory lanes use python tracemalloc peak allocations, so they reflect owned likelihood context and cache growth rather than full process resident memory",
            "memory growth status is classified from both peak-memory ratio and absolute byte growth between governed size classes so trivial allocation noise does not appear as a blowup",
            "the default governed memory roster samples jc69 and k80; f81 and hky85 remain available through explicit model overrides when broader model memory review is needed",
            "the governed native memory roster currently excludes fixed-parameter gtr because fixed-tree native inference still requires caller-supplied exchangeabilities for that lane",
        ],
    )


def benchmark_maximum_likelihood_wrapper_comparison() -> (
    MaximumLikelihoodWrapperBenchmarkReport
):
    """Benchmark native ML inference against cached IQ-TREE/FastTree/RAxML-style evidence."""
    native_tree = _measure_native_tree_inference()
    native_support = _measure_native_bootstrap_support()
    observations = [
        _build_iqtree_maximum_likelihood_observation(native_tree),
        _build_iqtree_bootstrap_support_observation(native_support),
        _build_fasttree_maximum_likelihood_observation(native_tree),
        _build_raxml_unsupported_observation(),
    ]
    summary_rows = _build_summary_rows(observations)
    tolerance_match_case_count = sum(
        1 for row in observations if row.status == "tolerance-match"
    )
    expected_model_assumption_difference_case_count = sum(
        1
        for row in observations
        if row.status == "expected-model-assumption-difference"
    )
    unsupported_case_count = sum(
        1 for row in observations if row.status == "unsupported-case"
    )
    native_bug_case_count = sum(1 for row in observations if row.status == "native-bug")
    blocking_case_count = sum(1 for row in observations if row.blocking)
    supported_case_count = sum(1 for row in observations if row.supported)
    return MaximumLikelihoodWrapperBenchmarkReport(
        observations=observations,
        summary_rows=summary_rows,
        case_count=len(observations),
        supported_case_count=supported_case_count,
        tolerance_match_case_count=tolerance_match_case_count,
        expected_model_assumption_difference_case_count=(
            expected_model_assumption_difference_case_count
        ),
        unsupported_case_count=unsupported_case_count,
        native_bug_case_count=native_bug_case_count,
        blocking_case_count=blocking_case_count,
        all_supported_cases_clear=blocking_case_count == 0,
        limitations=[
            "native runtime and peak memory are measured in-process with python tracemalloc, while cached wrapper peak memory comes from macos /usr/bin/time maximum resident set size; absolute memory values therefore remain lane-specific rather than same-instrument parity",
            "the iqtree support lane compares one canonical unrooted split against rooted native bootstrap clades so equivalent internal edges match even when rooting conventions differ",
            "the fasttree lane remains an expected model-assumption difference because FastTree uses an approximate JC+CAT likelihood surface and the cached tiny-alignment support output does not expose a separate comparable split-support table",
            "the raxml lane remains explicitly unsupported until one governed small-reference corpus with runtime, memory, likelihood, topology, and support evidence is checked in",
        ],
    )


def benchmark_native_maximum_likelihood_suite(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
    scaling_model_names: list[str] | tuple[str, ...] | None = None,
    scaling_search_methods: list[str] | tuple[str, ...] | None = None,
    accuracy_model_names: list[str] | tuple[str, ...] | None = None,
    accuracy_search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodSuiteReport:
    """Benchmark native ML scaling, wrapper parity, and truth recovery together."""
    speed_report = benchmark_native_maximum_likelihood_speed(
        replicates=replicates,
        size_classes=size_classes,
        model_names=scaling_model_names,
        search_methods=scaling_search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    memory_report = benchmark_native_maximum_likelihood_memory(
        replicates=replicates,
        size_classes=size_classes,
        model_names=scaling_model_names,
        search_methods=scaling_search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    accuracy_report = benchmark_native_maximum_likelihood_accuracy(
        model_names=accuracy_model_names,
        search_methods=accuracy_search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    wrapper_report = benchmark_maximum_likelihood_wrapper_comparison()
    observations = _build_native_maximum_likelihood_suite_observations(
        wrapper_report,
        accuracy_report,
    )
    summary_rows = _build_native_maximum_likelihood_suite_summary_rows(observations)
    native_advantage_case_count = sum(
        1 for row in observations if row.status == "native-advantage"
    )
    native_bug_case_count = sum(
        1 for row in observations if row.status == "native-bug"
    )
    expected_assumption_difference_case_count = sum(
        1
        for row in observations
        if row.status == "expected-assumption-difference"
    )
    unsupported_case_count = sum(
        1 for row in observations if row.status == "unsupported-case"
    )
    return NativeMaximumLikelihoodSuiteReport(
        speed_report=speed_report,
        memory_report=memory_report,
        accuracy_report=accuracy_report,
        wrapper_report=wrapper_report,
        observations=observations,
        summary_rows=summary_rows,
        native_advantage_case_count=native_advantage_case_count,
        native_bug_case_count=native_bug_case_count,
        expected_assumption_difference_case_count=(
            expected_assumption_difference_case_count
        ),
        unsupported_case_count=unsupported_case_count,
        all_supported_cases_clear=native_bug_case_count == 0,
        limitations=[
            "the suite uses native scaling reports as supporting evidence while the suite comparison rows stay anchored to wrapper and truth cases whose measurements are governed together",
            "speed rows treat slower native execution against compiled or approximate wrapper tools as an expected assumption difference rather than an automatic correctness bug",
            *speed_report.limitations,
            *memory_report.limitations,
            *accuracy_report.limitations,
            *wrapper_report.limitations,
        ],
    )


def write_native_maximum_likelihood_speed_bundle(
    output_root: Path,
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodSpeedBenchmarkBundle:
    """Write one governed artifact bundle for the native ML speed benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_native_maximum_likelihood_speed(
        replicates=replicates,
        size_classes=size_classes,
        model_names=model_names,
        search_methods=search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    summary_path = write_native_maximum_likelihood_speed_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = write_native_maximum_likelihood_speed_observation_table(
        output_root / "observations.tsv",
        report,
    )
    trend_table_path = write_native_maximum_likelihood_speed_trend_table(
        output_root / "trends.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_speed_run_json(
        output_root / "run.json",
        report,
    )
    return NativeMaximumLikelihoodSpeedBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
        trend_table_path=trend_table_path,
        run_json_path=run_json_path,
    )


def write_native_maximum_likelihood_accuracy_bundle(
    output_root: Path,
    *,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodAccuracyBenchmarkBundle:
    """Write one governed artifact bundle for the native ML accuracy benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_native_maximum_likelihood_accuracy(
        model_names=model_names,
        search_methods=search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    summary_path = write_native_maximum_likelihood_accuracy_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = write_native_maximum_likelihood_accuracy_observation_table(
        output_root / "observations.tsv",
        report,
    )
    wrapper_observation_table_path = write_native_maximum_likelihood_accuracy_wrapper_table(
        output_root / "wrapper_observations.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_accuracy_run_json(
        output_root / "run.json",
        report,
    )
    return NativeMaximumLikelihoodAccuracyBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
        wrapper_observation_table_path=wrapper_observation_table_path,
        run_json_path=run_json_path,
    )


def write_native_maximum_likelihood_memory_bundle(
    output_root: Path,
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
    model_names: list[str] | tuple[str, ...] | None = None,
    search_methods: list[str] | tuple[str, ...] | None = None,
    start_tree_count: int = 1,
    start_tree_seed: int = 17,
    upper_branch_length_bound: float = 1.0,
) -> NativeMaximumLikelihoodMemoryBenchmarkBundle:
    """Write one governed artifact bundle for the native ML memory benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_native_maximum_likelihood_memory(
        replicates=replicates,
        size_classes=size_classes,
        model_names=model_names,
        search_methods=search_methods,
        start_tree_count=start_tree_count,
        start_tree_seed=start_tree_seed,
        upper_branch_length_bound=upper_branch_length_bound,
    )
    summary_path = write_native_maximum_likelihood_memory_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = write_native_maximum_likelihood_memory_observation_table(
        output_root / "observations.tsv",
        report,
    )
    trend_table_path = write_native_maximum_likelihood_memory_trend_table(
        output_root / "trends.tsv",
        report,
    )
    run_json_path = write_native_maximum_likelihood_memory_run_json(
        output_root / "run.json",
        report,
    )
    return NativeMaximumLikelihoodMemoryBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
        trend_table_path=trend_table_path,
        run_json_path=run_json_path,
    )


def write_maximum_likelihood_wrapper_comparison_bundle(
    output_root: Path,
) -> MaximumLikelihoodWrapperBenchmarkBundle:
    """Write one governed artifact bundle for the maximum-likelihood wrapper benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    report = benchmark_maximum_likelihood_wrapper_comparison()
    summary_path = write_maximum_likelihood_wrapper_comparison_summary_table(
        output_root / "summary.tsv",
        report,
    )
    observation_table_path = (
        write_maximum_likelihood_wrapper_comparison_observation_table(
            output_root / "observations.tsv",
            report,
        )
    )
    run_json_path = write_maximum_likelihood_wrapper_comparison_run_json(
        output_root / "run.json",
        report,
    )
    return MaximumLikelihoodWrapperBenchmarkBundle(
        output_root=output_root,
        summary_path=summary_path,
        observation_table_path=observation_table_path,
        run_json_path=run_json_path,
    )


def _measure_native_tree_inference() -> _MeasuredNativeTreeInference:
    alignment_path = _alignment_path()

    def run():
        return infer_nucleotide_likelihood_tree_from_alignment(
            alignment_path,
            model_name="jc69",
            search_method="nni",
            start_tree_count=3,
            start_tree_seed=17,
            upper_branch_length_bound=10.0,
        )

    report, runtime_seconds, peak_memory_bytes = _measure_callback(run)
    return _MeasuredNativeTreeInference(
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        log_likelihood=report.best_final_log_likelihood,
        tree_newick=report.best_final_tree_newick,
        topology_signature=_unrooted_topology_signature(report.best_final_tree_newick),
    )


def _measure_native_bootstrap_support() -> _MeasuredNativeBootstrapSupport:
    alignment_path = _alignment_path()

    def run():
        return bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
            alignment_path,
            model_name="jc69",
            search_method="nni",
            start_tree_count=2,
            start_tree_seed=5,
            replicate_count=1,
            bootstrap_seed=9,
            upper_branch_length_bound=10.0,
        )

    report, runtime_seconds, peak_memory_bytes = _measure_callback(run)
    return _MeasuredNativeBootstrapSupport(
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        log_likelihood=report.reference_log_likelihood,
        tree_newick=report.reference_tree_newick,
        topology_signature=_unrooted_topology_signature(report.reference_tree_newick),
        support_by_split_id=_native_support_by_split_id(report),
    )


def _measure_callback(callback):
    tracemalloc.start()
    started = time.perf_counter()
    try:
        report = callback()
    finally:
        elapsed_seconds = time.perf_counter() - started
        _, peak_memory_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    return report, round(elapsed_seconds, 15), peak_memory_bytes


def _build_native_maximum_likelihood_accuracy_observation(
    case,
    *,
    search_method: str,
    inference_report,
    true_topology_log_likelihood: float,
) -> NativeMaximumLikelihoodAccuracyBenchmarkObservation:
    topology_recovered, rf_distance, normalized_rf = _compare_topology(
        case.tree_newick,
        inference_report.best_final_tree_newick,
    )
    estimated_parameter_values = _extract_native_accuracy_parameter_estimates(
        inference_report.model_selection_report.rows[0].parameter_values,
        parameter_names=tuple(case.parameter_truth),
    )
    parameter_absolute_deltas = {
        parameter_name: round(
            abs(estimated_parameter_values[parameter_name] - truth_value),
            15,
        )
        for parameter_name, truth_value in case.parameter_truth.items()
    }
    return NativeMaximumLikelihoodAccuracyBenchmarkObservation(
        case_id=case.case_id,
        model_name=case.model_name,
        search_method=search_method,
        taxon_count=inference_report.taxon_count,
        site_count=inference_report.site_count,
        simulation_seed=case.simulation_seed,
        selected_model_name=inference_report.selected_model_name,
        true_topology_fingerprint=_topology_signature_label(case.tree_newick),
        inferred_topology_fingerprint=_topology_signature_label(
            inference_report.best_final_tree_newick
        ),
        topology_recovered=topology_recovered,
        topology_recovery_status=(
            "topology-recovered" if topology_recovered else "topology-miss"
        ),
        robinson_foulds_distance=rf_distance,
        normalized_robinson_foulds=normalized_rf,
        parameter_truth=dict(case.parameter_truth),
        estimated_parameter_values=estimated_parameter_values,
        parameter_absolute_deltas=parameter_absolute_deltas,
        parameter_recovery_status=_classify_native_accuracy_parameter_recovery(
            parameter_absolute_deltas,
            parameter_tolerances=case.parameter_tolerances,
        ),
        true_topology_log_likelihood=true_topology_log_likelihood,
        inferred_log_likelihood=inference_report.best_final_log_likelihood,
        likelihood_shortfall=round(
            max(
                0.0,
                true_topology_log_likelihood - inference_report.best_final_log_likelihood,
            ),
            15,
        ),
        likelihood_quality_status=(
            "no-likelihood-shortfall"
            if max(
                0.0,
                true_topology_log_likelihood - inference_report.best_final_log_likelihood,
            )
            <= _NATIVE_MAXIMUM_LIKELIHOOD_ACCURACY_SHORTFALL_TOLERANCE
            else "likelihood-shortfall"
        ),
    )


def _resolve_native_maximum_likelihood_accuracy_model_names(
    model_names: list[str] | tuple[str, ...] | None,
) -> list[str]:
    roster = model_names or ("k80",)
    resolved = _validate_native_maximum_likelihood_benchmark_model_names(roster)
    available = {
        case.model_name for case in governed_native_maximum_likelihood_accuracy_cases()
    }
    unsupported = [model_name for model_name in resolved if model_name not in available]
    if unsupported:
        supported = ", ".join(sorted(available))
        raise ValueError(
            "native ML accuracy benchmark supports simulated truth cases for "
            f"{supported}; got {', '.join(unsupported)}"
        )
    return resolved


def _optimize_native_accuracy_truth_topology(
    case,
    records,
    *,
    upper_branch_length_bound: float,
) -> float:
    truth_tree = build_native_maximum_likelihood_accuracy_truth_tree(case)
    optimization_report = optimize_fixed_topology_nucleotide_branch_lengths(
        truth_tree,
        records,
        model_name=case.model_name,
        kappa=case.parameter_truth.get("kappa"),
        base_frequencies=_native_accuracy_base_frequencies(case),
        upper_branch_length_bound=upper_branch_length_bound,
    )
    return optimization_report.optimized_log_likelihood


def _native_accuracy_base_frequencies(case) -> dict[str, float] | None:
    if case.model_name in {"jc69", "k80"}:
        return None
    return {
        "A": case.parameter_truth["base_frequency_a"],
        "C": case.parameter_truth["base_frequency_c"],
        "G": case.parameter_truth["base_frequency_g"],
        "T": case.parameter_truth["base_frequency_t"],
    }


def _extract_native_accuracy_parameter_estimates(
    parameter_values: dict[str, float],
    *,
    parameter_names: tuple[str, ...],
) -> dict[str, float]:
    return {
        parameter_name: float(parameter_values[parameter_name])
        for parameter_name in parameter_names
    }


def _classify_native_accuracy_parameter_recovery(
    parameter_absolute_deltas: dict[str, float],
    *,
    parameter_tolerances: dict[str, float],
) -> str:
    if not parameter_absolute_deltas:
        return "parameter-not-applicable"
    return (
        "parameter-recovered"
        if all(
            parameter_absolute_deltas[parameter_name]
            <= parameter_tolerances[parameter_name]
            for parameter_name in parameter_absolute_deltas
        )
        else "parameter-miss"
    )


def _topology_signature_label(tree_newick: str) -> str:
    return " | ".join(_unrooted_topology_signature(tree_newick))


def _measure_native_maximum_likelihood_benchmark_lanes(
    *,
    size_classes: list[tuple[str, int, int]],
    model_names: list[str],
    search_methods: list[str],
    replicates: int,
    start_tree_count: int,
    start_tree_seed: int,
    upper_branch_length_bound: float,
    temporary_directory_prefix: str,
) -> list[_MeasuredNativeMaximumLikelihoodLane]:
    measured_lanes: list[_MeasuredNativeMaximumLikelihoodLane] = []
    with tempfile.TemporaryDirectory(prefix=temporary_directory_prefix) as tmpdir:
        tmp_path = Path(tmpdir)
        for size_class, taxon_count, site_count in size_classes:
            alignment_path = write_large_nucleotide_alignment(
                tmp_path / f"{size_class}.fasta",
                sequence_count=taxon_count,
                sequence_length=site_count,
            )
            for model_name in model_names:
                for search_method in search_methods:
                    measured_lanes.append(
                        _measure_native_maximum_likelihood_benchmark_lane(
                            alignment_path,
                            size_class=size_class,
                            model_name=model_name,
                            search_method=search_method,
                            taxon_count=taxon_count,
                            site_count=site_count,
                            replicates=replicates,
                            start_tree_count=start_tree_count,
                            start_tree_seed=start_tree_seed,
                            upper_branch_length_bound=upper_branch_length_bound,
                        )
                    )
    return measured_lanes


def _measure_native_maximum_likelihood_benchmark_lane(
    alignment_path: Path,
    *,
    size_class: str,
    model_name: str,
    search_method: str,
    taxon_count: int,
    site_count: int,
    replicates: int,
    start_tree_count: int,
    start_tree_seed: int,
    upper_branch_length_bound: float,
) -> NativeMaximumLikelihoodSpeedBenchmarkObservation:
    runtimes: list[float] = []
    peak_memory_bytes = 0
    report = None
    for _ in range(replicates):
        report, runtime_seconds, observed_peak_memory_bytes = _measure_callback(
            lambda: infer_nucleotide_likelihood_tree_from_alignment(
                alignment_path,
                model_name=model_name,
                search_method=search_method,
                start_tree_count=start_tree_count,
                start_tree_seed=start_tree_seed,
                upper_branch_length_bound=upper_branch_length_bound,
            )
        )
        runtimes.append(runtime_seconds)
        peak_memory_bytes = max(peak_memory_bytes, observed_peak_memory_bytes)
    assert report is not None
    return _MeasuredNativeMaximumLikelihoodLane(
        size_class=size_class,
        model_name=model_name,
        search_method=search_method,
        taxon_count=taxon_count,
        site_count=site_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory_bytes,
        memory_observation_kind="python-tracemalloc",
        selected_model_name=report.selected_model_name,
        best_run_source_label=report.best_run_source_label,
        search_algorithm=report.best_search_report.algorithm,
        final_log_likelihood=report.best_final_log_likelihood,
    )


def _build_native_maximum_likelihood_speed_observation(
    lane: _MeasuredNativeMaximumLikelihoodLane,
) -> NativeMaximumLikelihoodSpeedBenchmarkObservation:
    return NativeMaximumLikelihoodSpeedBenchmarkObservation(
        label=(
            f"{lane.size_class}-{lane.model_name}-{lane.search_method}-"
            f"{lane.taxon_count}-taxa-{lane.site_count}-sites"
        ),
        size_class=lane.size_class,
        model_name=lane.model_name,
        search_method=lane.search_method,
        taxon_count=lane.taxon_count,
        site_count=lane.site_count,
        runtime_seconds=lane.runtime_seconds,
        peak_memory_bytes=lane.peak_memory_bytes,
        selected_model_name=lane.selected_model_name,
        best_run_source_label=lane.best_run_source_label,
        search_algorithm=lane.search_algorithm,
        final_log_likelihood=lane.final_log_likelihood,
    )


def _build_native_maximum_likelihood_memory_observation(
    lane: _MeasuredNativeMaximumLikelihoodLane,
) -> NativeMaximumLikelihoodMemoryBenchmarkObservation:
    return NativeMaximumLikelihoodMemoryBenchmarkObservation(
        label=(
            f"{lane.size_class}-{lane.model_name}-{lane.search_method}-"
            f"{lane.taxon_count}-taxa-{lane.site_count}-sites"
        ),
        size_class=lane.size_class,
        model_name=lane.model_name,
        search_method=lane.search_method,
        taxon_count=lane.taxon_count,
        site_count=lane.site_count,
        runtime_seconds=lane.runtime_seconds,
        peak_memory_bytes=lane.peak_memory_bytes,
        memory_observation_kind=lane.memory_observation_kind,
        selected_model_name=lane.selected_model_name,
        best_run_source_label=lane.best_run_source_label,
        search_algorithm=lane.search_algorithm,
        final_log_likelihood=lane.final_log_likelihood,
    )


def _build_native_maximum_likelihood_speed_trend_rows(
    observations: list[NativeMaximumLikelihoodSpeedBenchmarkObservation],
    *,
    size_class_labels: list[str],
    model_names: list[str] | tuple[str, ...],
    search_methods: list[str],
) -> list[NativeMaximumLikelihoodSpeedTrendRow]:
    observation_by_key = {
        (row.model_name, row.search_method, row.size_class): row for row in observations
    }
    trend_rows: list[NativeMaximumLikelihoodSpeedTrendRow] = []
    for baseline_size_class, comparison_size_class in zip(
        size_class_labels,
        size_class_labels[1:],
        strict=False,
    ):
        for model_name in model_names:
            for search_method in search_methods:
                baseline = observation_by_key[
                    (model_name, search_method, baseline_size_class)
                ]
                comparison = observation_by_key[
                    (model_name, search_method, comparison_size_class)
                ]
                trend_rows.append(
                    NativeMaximumLikelihoodSpeedTrendRow(
                        model_name=model_name,
                        search_method=search_method,
                        baseline_size_class=baseline_size_class,
                        comparison_size_class=comparison_size_class,
                        baseline_taxon_count=baseline.taxon_count,
                        comparison_taxon_count=comparison.taxon_count,
                        baseline_site_count=baseline.site_count,
                        comparison_site_count=comparison.site_count,
                        baseline_runtime_seconds=baseline.runtime_seconds,
                        comparison_runtime_seconds=comparison.runtime_seconds,
                        runtime_delta_seconds=round(
                            comparison.runtime_seconds - baseline.runtime_seconds,
                            15,
                        ),
                        runtime_ratio=round(
                            comparison.runtime_seconds / baseline.runtime_seconds,
                            15,
                        ),
                    )
                )
    return trend_rows


def _build_native_maximum_likelihood_memory_trend_rows(
    observations: list[NativeMaximumLikelihoodMemoryBenchmarkObservation],
    *,
    size_class_labels: list[str],
    model_names: list[str] | tuple[str, ...],
    search_methods: list[str],
) -> list[NativeMaximumLikelihoodMemoryTrendRow]:
    observation_by_key = {
        (row.model_name, row.search_method, row.size_class): row for row in observations
    }
    trend_rows: list[NativeMaximumLikelihoodMemoryTrendRow] = []
    for baseline_size_class, comparison_size_class in zip(
        size_class_labels,
        size_class_labels[1:],
        strict=False,
    ):
        for model_name in model_names:
            for search_method in search_methods:
                baseline = observation_by_key[
                    (model_name, search_method, baseline_size_class)
                ]
                comparison = observation_by_key[
                    (model_name, search_method, comparison_size_class)
                ]
                peak_memory_ratio = round(
                    comparison.peak_memory_bytes
                    / max(baseline.peak_memory_bytes, 1),
                    15,
                )
                peak_memory_delta_bytes = (
                    comparison.peak_memory_bytes - baseline.peak_memory_bytes
                )
                trend_rows.append(
                    NativeMaximumLikelihoodMemoryTrendRow(
                        model_name=model_name,
                        search_method=search_method,
                        baseline_size_class=baseline_size_class,
                        comparison_size_class=comparison_size_class,
                        baseline_taxon_count=baseline.taxon_count,
                        comparison_taxon_count=comparison.taxon_count,
                        baseline_site_count=baseline.site_count,
                        comparison_site_count=comparison.site_count,
                        baseline_peak_memory_bytes=baseline.peak_memory_bytes,
                        comparison_peak_memory_bytes=comparison.peak_memory_bytes,
                        peak_memory_delta_bytes=peak_memory_delta_bytes,
                        peak_memory_ratio=peak_memory_ratio,
                        memory_growth_status=_classify_native_maximum_likelihood_memory_growth(
                            peak_memory_ratio=peak_memory_ratio,
                            peak_memory_delta_bytes=peak_memory_delta_bytes,
                        ),
                    )
                )
    return trend_rows


def _classify_native_maximum_likelihood_memory_growth(
    *,
    peak_memory_ratio: float,
    peak_memory_delta_bytes: int,
) -> str:
    if (
        peak_memory_ratio >= _MEMORY_BLOWUP_RATIO_THRESHOLD
        and peak_memory_delta_bytes >= _MEMORY_BLOWUP_DELTA_THRESHOLD_BYTES
    ):
        return "memory-blowup"
    if (
        peak_memory_ratio >= _ELEVATED_MEMORY_RATIO_THRESHOLD
        and peak_memory_delta_bytes >= _ELEVATED_MEMORY_DELTA_THRESHOLD_BYTES
    ):
        return "elevated-memory-growth"
    return "stable-memory-growth"


def _validate_native_maximum_likelihood_benchmark_size_classes(
    size_classes: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    if not size_classes:
        raise ValueError("size_classes must contain at least one size class")
    validated_size_classes: list[tuple[str, int, int]] = []
    seen_labels: set[str] = set()
    for size_class, taxon_count, site_count in size_classes:
        if not size_class.strip():
            raise ValueError("size class labels must be non-empty")
        if size_class in seen_labels:
            raise ValueError(f"duplicate size class label: {size_class}")
        if taxon_count < 4:
            raise ValueError("native ML speed benchmark requires at least 4 taxa")
        if site_count < 2:
            raise ValueError("native ML speed benchmark requires at least 2 sites")
        seen_labels.add(size_class)
        validated_size_classes.append((size_class, taxon_count, site_count))
    return validated_size_classes


def _validate_native_maximum_likelihood_benchmark_model_names(
    model_names: list[str] | tuple[str, ...],
) -> list[str]:
    if not model_names:
        raise ValueError("model_names must contain at least one fixed model")
    validated_model_names: list[str] = []
    for model_name in model_names:
        normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
        if normalized_model_name not in _SUPPORTED_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES:
            supported = ", ".join(_SUPPORTED_NATIVE_MAXIMUM_LIKELIHOOD_SPEED_MODEL_NAMES)
            raise ValueError(
                "native ML speed benchmark supports fixed models "
                f"{supported}; got {model_name!r}"
            )
        validated_model_names.append(normalized_model_name)
    return validated_model_names


def _build_iqtree_maximum_likelihood_observation(
    native_tree: _MeasuredNativeTreeInference,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "iqtree-jc69-maximum-likelihood"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    topology_equal, rf_distance, normalized_rf = _compare_topology(
        native_tree.tree_newick,
        wrapper_tree_newick,
    )
    wrapper_log_likelihood = float(reference["log_likelihood"])
    log_likelihood_delta = native_tree.log_likelihood - wrapper_log_likelihood
    status = (
        "tolerance-match"
        if topology_equal and abs(log_likelihood_delta) <= _LIKELIHOOD_TOLERANCE
        else "native-bug"
    )
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="iqtree-jc69-maximum-likelihood",
        wrapper_engine="IQ-TREE-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="cached IQ-TREE JC maximum-likelihood reference",
        comparison_policy="topology-and-likelihood-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        native_runtime_seconds=native_tree.runtime_seconds,
        native_peak_memory_bytes=native_tree.peak_memory_bytes,
        native_memory_observation_kind=native_tree.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_tree.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_tree.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_tree.log_likelihood,
        wrapper_log_likelihood=wrapper_log_likelihood,
        log_likelihood_delta=log_likelihood_delta,
        native_tree_newick=native_tree.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_tree.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=rf_distance,
        normalized_robinson_foulds=normalized_rf,
        topology_equal=topology_equal,
        support_comparisons=[],
        notes=list(reference["notes"]),
    )


def _build_iqtree_bootstrap_support_observation(
    native_support: _MeasuredNativeBootstrapSupport,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "iqtree-jc69-ultrafast-bootstrap-support"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    support_comparisons = _support_comparisons(
        native_support.support_by_split_id,
        reference["support_rows"],
    )
    max_support_delta = max(
        (abs(row.support_percent_delta) for row in support_comparisons),
        default=0.0,
    )
    status = (
        "tolerance-match"
        if max_support_delta <= _SUPPORT_TOLERANCE and support_comparisons
        else "native-bug"
    )
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="iqtree-jc69-ultrafast-bootstrap-support",
        wrapper_engine="IQ-TREE-style",
        native_surface="bootstrap_nucleotide_likelihood_tree_inference_from_alignment",
        wrapper_surface="cached IQ-TREE JC ultrafast-bootstrap reference",
        comparison_policy="canonical-split-support-tolerance",
        status=status,
        supported=True,
        blocking=status == "native-bug",
        native_runtime_seconds=native_support.runtime_seconds,
        native_peak_memory_bytes=native_support.peak_memory_bytes,
        native_memory_observation_kind=native_support.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_support.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_support.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_support.log_likelihood,
        wrapper_log_likelihood=float(reference["log_likelihood"]),
        log_likelihood_delta=(
            native_support.log_likelihood - float(reference["log_likelihood"])
        ),
        native_tree_newick=native_support.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_support.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=0,
        normalized_robinson_foulds=0.0,
        topology_equal=True,
        support_comparisons=support_comparisons,
        notes=list(reference["notes"]),
    )


def _build_fasttree_maximum_likelihood_observation(
    native_tree: _MeasuredNativeTreeInference,
) -> MaximumLikelihoodWrapperBenchmarkObservation:
    reference = MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
        "fasttree-jc69-approximate-maximum-likelihood"
    ]
    wrapper_tree_newick = str(reference["tree_newick"])
    topology_equal, rf_distance, normalized_rf = _compare_topology(
        native_tree.tree_newick,
        wrapper_tree_newick,
    )
    wrapper_log_likelihood = float(reference["log_likelihood"])
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="fasttree-jc69-approximate-maximum-likelihood",
        wrapper_engine="FastTree-style",
        native_surface="infer_nucleotide_likelihood_tree_from_alignment",
        wrapper_surface="cached FastTree JC+CAT approximate ML reference",
        comparison_policy="expected-assumption-difference",
        status="expected-model-assumption-difference",
        supported=True,
        blocking=False,
        native_runtime_seconds=native_tree.runtime_seconds,
        native_peak_memory_bytes=native_tree.peak_memory_bytes,
        native_memory_observation_kind=native_tree.memory_observation_kind,
        wrapper_runtime_seconds=float(reference["runtime_seconds"]),
        wrapper_peak_memory_bytes=int(reference["peak_memory_bytes"]),
        wrapper_memory_observation_kind=str(reference["memory_observation_kind"]),
        runtime_delta_seconds=(
            native_tree.runtime_seconds - float(reference["runtime_seconds"])
        ),
        peak_memory_delta_bytes=(
            native_tree.peak_memory_bytes - int(reference["peak_memory_bytes"])
        ),
        native_log_likelihood=native_tree.log_likelihood,
        wrapper_log_likelihood=wrapper_log_likelihood,
        log_likelihood_delta=native_tree.log_likelihood - wrapper_log_likelihood,
        native_tree_newick=native_tree.tree_newick,
        wrapper_tree_newick=wrapper_tree_newick,
        native_topology_signature=native_tree.topology_signature,
        wrapper_topology_signature=_unrooted_topology_signature(wrapper_tree_newick),
        robinson_foulds_distance=rf_distance,
        normalized_robinson_foulds=normalized_rf,
        topology_equal=topology_equal,
        support_comparisons=[],
        notes=list(reference["notes"]),
    )


def _build_raxml_unsupported_observation() -> (
    MaximumLikelihoodWrapperBenchmarkObservation
):
    return MaximumLikelihoodWrapperBenchmarkObservation(
        case_id="raxml-small-maximum-likelihood-reference",
        wrapper_engine="RAxML-style",
        native_surface="none",
        wrapper_surface="governed RAxML small ML benchmark reference",
        comparison_policy="unsupported",
        status="unsupported-case",
        supported=False,
        blocking=False,
        native_runtime_seconds=None,
        native_peak_memory_bytes=None,
        native_memory_observation_kind=None,
        wrapper_runtime_seconds=None,
        wrapper_peak_memory_bytes=None,
        wrapper_memory_observation_kind=None,
        runtime_delta_seconds=None,
        peak_memory_delta_bytes=None,
        native_log_likelihood=None,
        wrapper_log_likelihood=None,
        log_likelihood_delta=None,
        native_tree_newick=None,
        wrapper_tree_newick=None,
        native_topology_signature=(),
        wrapper_topology_signature=(),
        robinson_foulds_distance=None,
        normalized_robinson_foulds=None,
        topology_equal=None,
        support_comparisons=[],
        notes=[
            "the repository does not yet ship one governed RAxML-style small-reference corpus with runtime, memory, likelihood, topology, and support evidence",
            "this lane remains explicitly unsupported until that cached reference corpus exists",
        ],
    )


def _support_comparisons(
    native_support_by_split_id: dict[str, float],
    wrapper_support_rows: object,
) -> list[MaximumLikelihoodWrapperSupportComparison]:
    wrapper_support_by_split_id = {
        canonical_clade_id(frozenset(row["split_taxa"])): float(row["support_percent"])
        for row in list(wrapper_support_rows)
    }
    split_ids = sorted(
        set(native_support_by_split_id) | set(wrapper_support_by_split_id),
        key=lambda value: split_sort_key(frozenset(value.split("|"))),
    )
    comparisons: list[MaximumLikelihoodWrapperSupportComparison] = []
    for split_id in split_ids:
        native_support_percent = native_support_by_split_id.get(split_id)
        wrapper_support_percent = wrapper_support_by_split_id.get(split_id)
        support_percent_delta = (
            None
            if native_support_percent is None or wrapper_support_percent is None
            else native_support_percent - wrapper_support_percent
        )
        comparisons.append(
            MaximumLikelihoodWrapperSupportComparison(
                split_taxa=split_id.split("|"),
                native_support_percent=native_support_percent,
                wrapper_support_percent=wrapper_support_percent,
                support_percent_delta=support_percent_delta,
            )
        )
    return comparisons


def _native_support_by_split_id(report) -> dict[str, float]:
    universe = set(loads_newick(report.reference_tree_newick).tip_names)
    support_by_split_id: dict[str, float] = {}
    for row in report.clade_support_rows:
        split = canonical_bipartition(set(row.descendant_taxa), universe)
        split_id = canonical_clade_id(split)
        existing = support_by_split_id.get(split_id)
        if existing is None:
            support_by_split_id[split_id] = row.support_percent
            continue
        if not math.isclose(existing, row.support_percent, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(
                f"native bootstrap reported inconsistent support for split {split_id}"
            )
    return support_by_split_id


def _compare_topology(
    native_tree_newick: str,
    wrapper_tree_newick: str,
) -> tuple[bool, int, float]:
    metrics = robinson_foulds_metrics(
        loads_newick(native_tree_newick),
        loads_newick(wrapper_tree_newick),
        shared_taxa=set(loads_newick(native_tree_newick).tip_names),
        rf_mode="unrooted",
    )
    return metrics.distance == 0, metrics.distance, metrics.normalized_distance


def _build_native_maximum_likelihood_suite_observations(
    wrapper_report: MaximumLikelihoodWrapperBenchmarkReport,
    accuracy_report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> list[NativeMaximumLikelihoodSuiteObservation]:
    observations: list[NativeMaximumLikelihoodSuiteObservation] = []
    for wrapper_row in wrapper_report.observations:
        observations.extend(
            _build_native_maximum_likelihood_suite_wrapper_rows(wrapper_row)
        )
    observations.extend(
        _build_native_maximum_likelihood_suite_accuracy_rows(accuracy_report)
    )
    return observations


def _build_native_maximum_likelihood_suite_wrapper_rows(
    row: MaximumLikelihoodWrapperBenchmarkObservation,
) -> list[NativeMaximumLikelihoodSuiteObservation]:
    observations: list[NativeMaximumLikelihoodSuiteObservation] = []
    if (
        row.supported
        and row.native_runtime_seconds is not None
        and row.wrapper_runtime_seconds is not None
    ):
        speed_status = (
            "native-advantage"
            if row.native_runtime_seconds <= row.wrapper_runtime_seconds
            else "expected-assumption-difference"
        )
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="speed",
                status=speed_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit="seconds",
                native_value=row.native_runtime_seconds,
                comparison_value=row.wrapper_runtime_seconds,
                delta_value=row.runtime_delta_seconds,
                notes=list(row.notes),
            )
        )
    if (
        row.supported
        and row.native_peak_memory_bytes is not None
        and row.wrapper_peak_memory_bytes is not None
    ):
        memory_status = (
            "native-advantage"
            if row.native_peak_memory_bytes <= row.wrapper_peak_memory_bytes
            else "expected-assumption-difference"
        )
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="memory",
                status=memory_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit="bytes",
                native_value=float(row.native_peak_memory_bytes),
                comparison_value=float(row.wrapper_peak_memory_bytes),
                delta_value=(
                    None
                    if row.peak_memory_delta_bytes is None
                    else float(row.peak_memory_delta_bytes)
                ),
                notes=list(row.notes),
            )
        )
    suite_status = _wrapper_status_to_suite_status(row.status)
    if row.native_log_likelihood is not None and row.wrapper_log_likelihood is not None:
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="likelihood",
                status=suite_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit="log-likelihood",
                native_value=row.native_log_likelihood,
                comparison_value=row.wrapper_log_likelihood,
                delta_value=row.log_likelihood_delta,
                notes=list(row.notes),
            )
        )
    elif not row.supported:
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="likelihood",
                status=suite_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit=None,
                native_value=None,
                comparison_value=None,
                delta_value=None,
                notes=list(row.notes),
            )
        )
    if row.topology_equal is not None:
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="topology",
                status=suite_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit="robinson-foulds-distance",
                native_value=(
                    None
                    if row.robinson_foulds_distance is None
                    else float(row.robinson_foulds_distance)
                ),
                comparison_value=0.0,
                delta_value=row.normalized_robinson_foulds,
                notes=list(row.notes),
            )
        )
    if row.support_comparisons:
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="support",
                status=suite_status,
                native_surface=row.native_surface,
                comparison_surface=row.wrapper_surface,
                measurement_unit="support-percent-delta",
                native_value=None,
                comparison_value=None,
                delta_value=_suite_max_support_percent_delta(row),
                notes=list(row.notes),
            )
        )
    return observations


def _build_native_maximum_likelihood_suite_accuracy_rows(
    report: NativeMaximumLikelihoodAccuracyBenchmarkReport,
) -> list[NativeMaximumLikelihoodSuiteObservation]:
    observations: list[NativeMaximumLikelihoodSuiteObservation] = []
    for row in report.observations:
        status = (
            "native-advantage"
            if row.topology_recovered
            and row.likelihood_quality_status == "no-likelihood-shortfall"
            and row.parameter_recovery_status
            in {"parameter-recovered", "parameter-not-applicable"}
            else "native-bug"
        )
        observations.append(
            NativeMaximumLikelihoodSuiteObservation(
                case_id=row.case_id,
                dimension="accuracy",
                status=status,
                native_surface="benchmark_native_maximum_likelihood_accuracy",
                comparison_surface="governed simulated truth",
                measurement_unit="normalized-robinson-foulds",
                native_value=row.normalized_robinson_foulds,
                comparison_value=0.0,
                delta_value=row.likelihood_shortfall,
                notes=[
                    f"parameter status: {row.parameter_recovery_status}",
                    f"likelihood status: {row.likelihood_quality_status}",
                ],
            )
        )
    return observations


def _wrapper_status_to_suite_status(status: str) -> str:
    if status == "tolerance-match":
        return "native-advantage"
    if status == "expected-model-assumption-difference":
        return "expected-assumption-difference"
    if status == "unsupported-case":
        return "unsupported-case"
    return "native-bug"


def _suite_max_support_percent_delta(
    row: MaximumLikelihoodWrapperBenchmarkObservation,
) -> float | None:
    deltas = [
        abs(comparison.support_percent_delta)
        for comparison in row.support_comparisons
        if comparison.support_percent_delta is not None
    ]
    if not deltas:
        return None
    return max(deltas)


def _unrooted_topology_signature(tree_newick: str) -> tuple[str, ...]:
    splits = informative_unrooted_splits(loads_newick(tree_newick))
    return tuple(
        canonical_clade_id(signature)
        for signature in sorted(splits, key=split_sort_key)
    )


def _build_native_maximum_likelihood_suite_summary_rows(
    observations: list[NativeMaximumLikelihoodSuiteObservation],
) -> list[NativeMaximumLikelihoodSuiteSummaryRow]:
    rows: list[NativeMaximumLikelihoodSuiteSummaryRow] = []
    for status in _NATIVE_MAXIMUM_LIKELIHOOD_SUITE_STATUSES:
        selected = [row for row in observations if row.status == status]
        rows.append(
            NativeMaximumLikelihoodSuiteSummaryRow(
                status=status,
                observation_count=len(selected),
                dimensions=sorted({row.dimension for row in selected}),
                case_ids=sorted({row.case_id for row in selected}),
            )
        )
    return rows


def _build_summary_rows(
    observations: list[MaximumLikelihoodWrapperBenchmarkObservation],
) -> list[MaximumLikelihoodWrapperBenchmarkSummaryRow]:
    rows: list[MaximumLikelihoodWrapperBenchmarkSummaryRow] = []
    for status in _MAXIMUM_LIKELIHOOD_WRAPPER_BENCHMARK_STATUSES:
        selected = [row for row in observations if row.status == status]
        rows.append(
            MaximumLikelihoodWrapperBenchmarkSummaryRow(
                status=status,
                case_count=len(selected),
                blocking_case_count=sum(1 for row in selected if row.blocking),
                case_ids=[row.case_id for row in selected],
            )
        )
    return rows


def _alignment_path() -> Path:
    return Path(__file__).resolve().parents[6] / str(
        MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS[
            "iqtree-jc69-maximum-likelihood"
        ]["alignment_path"]
    )
