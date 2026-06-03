from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    IndependentMetropolisHastingsChainTraceAutocorrelationReport,
    IndependentMetropolisHastingsTraceAutocorrelationReport,
    MetropolisHastingsTraceAutocorrelationReport,
    TraceAutocorrelationLagRow,
    TraceAutocorrelationParameterReport,
    compute_trace_autocorrelation,
    summarize_independent_metropolis_hastings_trace_autocorrelation,
    summarize_metropolis_hastings_trace_autocorrelation,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    IndependentMetropolisHastingsChainTraceAutocorrelationReport as IndependentMetropolisHastingsChainTraceAutocorrelationReportImpl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    IndependentMetropolisHastingsTraceAutocorrelationReport as IndependentMetropolisHastingsTraceAutocorrelationReportImpl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    MetropolisHastingsTraceAutocorrelationReport as MetropolisHastingsTraceAutocorrelationReportImpl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    TraceAutocorrelationLagRow as TraceAutocorrelationLagRowImpl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    TraceAutocorrelationParameterReport as TraceAutocorrelationParameterReportImpl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    compute_trace_autocorrelation as compute_trace_autocorrelation_impl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    summarize_independent_metropolis_hastings_trace_autocorrelation as summarize_independent_metropolis_hastings_trace_autocorrelation_impl,
)
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    summarize_metropolis_hastings_trace_autocorrelation as summarize_metropolis_hastings_trace_autocorrelation_impl,
)


def test_bayesian_exports_trace_autocorrelation_surface() -> None:
    assert (
        IndependentMetropolisHastingsChainTraceAutocorrelationReport
        is IndependentMetropolisHastingsChainTraceAutocorrelationReportImpl
    )
    assert (
        IndependentMetropolisHastingsTraceAutocorrelationReport
        is IndependentMetropolisHastingsTraceAutocorrelationReportImpl
    )
    assert (
        MetropolisHastingsTraceAutocorrelationReport
        is MetropolisHastingsTraceAutocorrelationReportImpl
    )
    assert TraceAutocorrelationLagRow is TraceAutocorrelationLagRowImpl
    assert (
        TraceAutocorrelationParameterReport is TraceAutocorrelationParameterReportImpl
    )
    assert compute_trace_autocorrelation is compute_trace_autocorrelation_impl
    assert (
        summarize_independent_metropolis_hastings_trace_autocorrelation
        is summarize_independent_metropolis_hastings_trace_autocorrelation_impl
    )
    assert (
        summarize_metropolis_hastings_trace_autocorrelation
        is summarize_metropolis_hastings_trace_autocorrelation_impl
    )
