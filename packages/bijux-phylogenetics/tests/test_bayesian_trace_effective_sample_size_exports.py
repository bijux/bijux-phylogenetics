from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport,
    IndependentMetropolisHastingsTraceEffectiveSampleSizeReport,
    MetropolisHastingsTraceEffectiveSampleSizeReport,
    TraceEffectiveSampleSizeRow,
    compute_trace_effective_sample_size,
    compute_trace_integrated_autocorrelation_time,
    summarize_independent_metropolis_hastings_trace_effective_sample_size,
    summarize_metropolis_hastings_trace_effective_sample_size,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport as IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReportImpl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    IndependentMetropolisHastingsTraceEffectiveSampleSizeReport as IndependentMetropolisHastingsTraceEffectiveSampleSizeReportImpl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    MetropolisHastingsTraceEffectiveSampleSizeReport as MetropolisHastingsTraceEffectiveSampleSizeReportImpl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    TraceEffectiveSampleSizeRow as TraceEffectiveSampleSizeRowImpl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    compute_trace_effective_sample_size as compute_trace_effective_sample_size_impl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    compute_trace_integrated_autocorrelation_time as compute_trace_integrated_autocorrelation_time_impl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    summarize_independent_metropolis_hastings_trace_effective_sample_size as summarize_independent_metropolis_hastings_trace_effective_sample_size_impl,
)
from bijux_phylogenetics.bayesian.trace_effective_sample_size import (
    summarize_metropolis_hastings_trace_effective_sample_size as summarize_metropolis_hastings_trace_effective_sample_size_impl,
)


def test_bayesian_exports_trace_effective_sample_size_surface() -> None:
    assert TraceEffectiveSampleSizeRow is TraceEffectiveSampleSizeRowImpl
    assert (
        MetropolisHastingsTraceEffectiveSampleSizeReport
        is MetropolisHastingsTraceEffectiveSampleSizeReportImpl
    )
    assert (
        IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport
        is IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReportImpl
    )
    assert (
        IndependentMetropolisHastingsTraceEffectiveSampleSizeReport
        is IndependentMetropolisHastingsTraceEffectiveSampleSizeReportImpl
    )
    assert (
        compute_trace_effective_sample_size is compute_trace_effective_sample_size_impl
    )
    assert (
        compute_trace_integrated_autocorrelation_time
        is compute_trace_integrated_autocorrelation_time_impl
    )
    assert (
        summarize_independent_metropolis_hastings_trace_effective_sample_size
        is summarize_independent_metropolis_hastings_trace_effective_sample_size_impl
    )
    assert (
        summarize_metropolis_hastings_trace_effective_sample_size
        is summarize_metropolis_hastings_trace_effective_sample_size_impl
    )
