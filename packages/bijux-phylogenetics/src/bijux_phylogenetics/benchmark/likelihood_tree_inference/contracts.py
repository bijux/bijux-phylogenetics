from __future__ import annotations

from dataclasses import dataclass


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
