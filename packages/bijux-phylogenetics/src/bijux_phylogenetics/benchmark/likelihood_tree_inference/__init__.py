"""Native maximum-likelihood tree-inference benchmarks."""

from .artifact_outputs import (
    write_maximum_likelihood_wrapper_comparison_observation_table,
    write_maximum_likelihood_wrapper_comparison_run_json,
    write_maximum_likelihood_wrapper_comparison_summary_table,
    write_native_maximum_likelihood_speed_observation_table,
    write_native_maximum_likelihood_speed_run_json,
    write_native_maximum_likelihood_speed_summary_table,
    write_native_maximum_likelihood_speed_trend_table,
)
from .builder import benchmark_native_maximum_likelihood_speed
from .builder import benchmark_maximum_likelihood_wrapper_comparison
from .builder import write_native_maximum_likelihood_speed_bundle
from .builder import write_maximum_likelihood_wrapper_comparison_bundle
from .contracts import (
    MaximumLikelihoodWrapperBenchmarkBundle,
    MaximumLikelihoodWrapperBenchmarkObservation,
    MaximumLikelihoodWrapperBenchmarkReport,
    MaximumLikelihoodWrapperBenchmarkSummaryRow,
    MaximumLikelihoodWrapperSupportComparison,
    NativeMaximumLikelihoodSpeedBenchmarkBundle,
    NativeMaximumLikelihoodSpeedBenchmarkObservation,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
    NativeMaximumLikelihoodSpeedTrendRow,
)
from .wrapper_references import (
    MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS,
)

__all__ = [
    "MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS",
    "MaximumLikelihoodWrapperBenchmarkBundle",
    "MaximumLikelihoodWrapperBenchmarkObservation",
    "MaximumLikelihoodWrapperBenchmarkReport",
    "MaximumLikelihoodWrapperBenchmarkSummaryRow",
    "MaximumLikelihoodWrapperSupportComparison",
    "NativeMaximumLikelihoodSpeedBenchmarkBundle",
    "NativeMaximumLikelihoodSpeedBenchmarkObservation",
    "NativeMaximumLikelihoodSpeedBenchmarkReport",
    "NativeMaximumLikelihoodSpeedTrendRow",
    "benchmark_native_maximum_likelihood_speed",
    "benchmark_maximum_likelihood_wrapper_comparison",
    "write_native_maximum_likelihood_speed_bundle",
    "write_native_maximum_likelihood_speed_observation_table",
    "write_native_maximum_likelihood_speed_run_json",
    "write_native_maximum_likelihood_speed_summary_table",
    "write_native_maximum_likelihood_speed_trend_table",
    "write_maximum_likelihood_wrapper_comparison_bundle",
    "write_maximum_likelihood_wrapper_comparison_observation_table",
    "write_maximum_likelihood_wrapper_comparison_run_json",
    "write_maximum_likelihood_wrapper_comparison_summary_table",
]
