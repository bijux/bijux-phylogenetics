from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    HighestPosteriorDensityInterval,
    IndependentMetropolisHastingsChainTracePosteriorIntervalReport,
    IndependentMetropolisHastingsTracePosteriorIntervalReport,
    MetropolisHastingsTracePosteriorIntervalReport,
    TracePosteriorIntervalRow,
    compute_equal_tail_interval,
    compute_highest_posterior_density_interval,
    summarize_independent_metropolis_hastings_trace_posterior_intervals,
    summarize_metropolis_hastings_trace_posterior_intervals,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    HighestPosteriorDensityInterval as HighestPosteriorDensityIntervalImpl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    IndependentMetropolisHastingsChainTracePosteriorIntervalReport as IndependentMetropolisHastingsChainTracePosteriorIntervalReportImpl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    IndependentMetropolisHastingsTracePosteriorIntervalReport as IndependentMetropolisHastingsTracePosteriorIntervalReportImpl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    MetropolisHastingsTracePosteriorIntervalReport as MetropolisHastingsTracePosteriorIntervalReportImpl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    TracePosteriorIntervalRow as TracePosteriorIntervalRowImpl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    compute_equal_tail_interval as compute_equal_tail_interval_impl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    compute_highest_posterior_density_interval as compute_highest_posterior_density_interval_impl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    summarize_independent_metropolis_hastings_trace_posterior_intervals as summarize_independent_metropolis_hastings_trace_posterior_intervals_impl,
)
from bijux_phylogenetics.bayesian.trace_posterior_intervals import (
    summarize_metropolis_hastings_trace_posterior_intervals as summarize_metropolis_hastings_trace_posterior_intervals_impl,
)


def test_bayesian_exports_trace_posterior_interval_surface() -> None:
    assert HighestPosteriorDensityInterval is HighestPosteriorDensityIntervalImpl
    assert TracePosteriorIntervalRow is TracePosteriorIntervalRowImpl
    assert (
        MetropolisHastingsTracePosteriorIntervalReport
        is MetropolisHastingsTracePosteriorIntervalReportImpl
    )
    assert (
        IndependentMetropolisHastingsChainTracePosteriorIntervalReport
        is IndependentMetropolisHastingsChainTracePosteriorIntervalReportImpl
    )
    assert (
        IndependentMetropolisHastingsTracePosteriorIntervalReport
        is IndependentMetropolisHastingsTracePosteriorIntervalReportImpl
    )
    assert compute_equal_tail_interval is compute_equal_tail_interval_impl
    assert (
        compute_highest_posterior_density_interval
        is compute_highest_posterior_density_interval_impl
    )
    assert (
        summarize_independent_metropolis_hastings_trace_posterior_intervals
        is summarize_independent_metropolis_hastings_trace_posterior_intervals_impl
    )
    assert (
        summarize_metropolis_hastings_trace_posterior_intervals
        is summarize_metropolis_hastings_trace_posterior_intervals_impl
    )
