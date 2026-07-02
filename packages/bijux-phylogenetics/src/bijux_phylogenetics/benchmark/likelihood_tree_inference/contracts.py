from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class NativeMaximumLikelihoodSpeedBenchmarkObservation:
    """One native ML speed observation for a fixed model, size class, and search mode."""

    label: str
    size_class: str
    model_name: str
    search_method: str
    taxon_count: int
    site_count: int
    runtime_seconds: float
    peak_memory_bytes: int
    selected_model_name: str
    best_run_source_label: str
    search_algorithm: str
    final_log_likelihood: float


@dataclass(slots=True)
class NativeMaximumLikelihoodSpeedTrendRow:
    """One runtime trend comparison between two governed ML speed size classes."""

    model_name: str
    search_method: str
    baseline_size_class: str
    comparison_size_class: str
    baseline_taxon_count: int
    comparison_taxon_count: int
    baseline_site_count: int
    comparison_site_count: int
    baseline_runtime_seconds: float
    comparison_runtime_seconds: float
    runtime_delta_seconds: float
    runtime_ratio: float


@dataclass(slots=True)
class NativeMaximumLikelihoodSpeedBenchmarkReport:
    """Native ML runtime benchmark across governed model, size, and search lanes."""

    replicates: int
    observations: list[NativeMaximumLikelihoodSpeedBenchmarkObservation]
    trend_rows: list[NativeMaximumLikelihoodSpeedTrendRow]
    size_class_labels: list[str]
    model_names: list[str]
    search_methods: list[str]
    limitations: list[str]


@dataclass(slots=True)
class NativeMaximumLikelihoodSpeedBenchmarkBundle:
    """Written artifact bundle for one native ML speed benchmark run."""

    output_root: Path
    summary_path: Path
    observation_table_path: Path
    trend_table_path: Path
    run_json_path: Path


@dataclass(slots=True)
class NativeMaximumLikelihoodMemoryBenchmarkObservation:
    """One native ML peak-memory observation for a fixed model, size class, and search mode."""

    label: str
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


@dataclass(slots=True)
class NativeMaximumLikelihoodMemoryTrendRow:
    """One peak-memory growth comparison between two governed ML memory size classes."""

    model_name: str
    search_method: str
    baseline_size_class: str
    comparison_size_class: str
    baseline_taxon_count: int
    comparison_taxon_count: int
    baseline_site_count: int
    comparison_site_count: int
    baseline_peak_memory_bytes: int
    comparison_peak_memory_bytes: int
    peak_memory_delta_bytes: int
    peak_memory_ratio: float
    memory_growth_status: str


@dataclass(slots=True)
class NativeMaximumLikelihoodMemoryBenchmarkReport:
    """Native ML peak-memory benchmark across governed model, size, and search lanes."""

    replicates: int
    observations: list[NativeMaximumLikelihoodMemoryBenchmarkObservation]
    trend_rows: list[NativeMaximumLikelihoodMemoryTrendRow]
    size_class_labels: list[str]
    model_names: list[str]
    search_methods: list[str]
    stable_memory_case_count: int
    elevated_memory_case_count: int
    memory_blowup_case_count: int
    limitations: list[str]


@dataclass(slots=True)
class NativeMaximumLikelihoodMemoryBenchmarkBundle:
    """Written artifact bundle for one native ML peak-memory benchmark run."""

    output_root: Path
    summary_path: Path
    observation_table_path: Path
    trend_table_path: Path
    run_json_path: Path


@dataclass(slots=True)
class MaximumLikelihoodWrapperSupportComparison:
    """One canonical split-support comparison between native and cached wrapper evidence."""

    split_taxa: list[str]
    native_support_percent: float | None
    wrapper_support_percent: float | None
    support_percent_delta: float | None


@dataclass(slots=True)
class MaximumLikelihoodWrapperBenchmarkObservation:
    """One native-versus-wrapper maximum-likelihood comparison case."""

    case_id: str
    wrapper_engine: str
    native_surface: str
    wrapper_surface: str
    comparison_policy: str
    status: str
    supported: bool
    blocking: bool
    native_runtime_seconds: float | None
    native_peak_memory_bytes: int | None
    native_memory_observation_kind: str | None
    wrapper_runtime_seconds: float | None
    wrapper_peak_memory_bytes: int | None
    wrapper_memory_observation_kind: str | None
    runtime_delta_seconds: float | None
    peak_memory_delta_bytes: int | None
    native_log_likelihood: float | None
    wrapper_log_likelihood: float | None
    log_likelihood_delta: float | None
    native_tree_newick: str | None
    wrapper_tree_newick: str | None
    native_topology_signature: tuple[str, ...]
    wrapper_topology_signature: tuple[str, ...]
    robinson_foulds_distance: int | None
    normalized_robinson_foulds: float | None
    topology_equal: bool | None
    support_comparisons: list[MaximumLikelihoodWrapperSupportComparison]
    notes: list[str]


@dataclass(slots=True)
class MaximumLikelihoodWrapperBenchmarkSummaryRow:
    """One status-level summary row across maximum-likelihood wrapper benchmark cases."""

    status: str
    case_count: int
    blocking_case_count: int
    case_ids: list[str]


@dataclass(slots=True)
class MaximumLikelihoodWrapperBenchmarkReport:
    """Native-versus-wrapper maximum-likelihood benchmark over cached reference lanes."""

    observations: list[MaximumLikelihoodWrapperBenchmarkObservation]
    summary_rows: list[MaximumLikelihoodWrapperBenchmarkSummaryRow]
    case_count: int
    supported_case_count: int
    tolerance_match_case_count: int
    expected_model_assumption_difference_case_count: int
    unsupported_case_count: int
    native_bug_case_count: int
    blocking_case_count: int
    all_supported_cases_clear: bool
    limitations: list[str]


@dataclass(slots=True)
class MaximumLikelihoodWrapperBenchmarkBundle:
    """Written artifact bundle for one native-versus-wrapper ML benchmark run."""

    output_root: Path
    summary_path: Path
    observation_table_path: Path
    run_json_path: Path
