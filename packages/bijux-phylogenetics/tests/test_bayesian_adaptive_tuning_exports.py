from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    AdaptiveMetropolisHastingsRunReport,
    AdaptiveTuningController,
    AdaptiveTuningReport,
    AdaptiveTuningWindowRow,
    build_adaptive_tuning_controller,
    build_adaptive_tuning_report,
    build_adaptive_tuning_window_row,
    run_adaptive_tuned_metropolis_hastings_sampler,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveMetropolisHastingsRunReport as AdaptiveMetropolisHastingsRunReportImpl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveTuningController as AdaptiveTuningControllerImpl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveTuningReport as AdaptiveTuningReportImpl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveTuningWindowRow as AdaptiveTuningWindowRowImpl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    build_adaptive_tuning_controller as build_adaptive_tuning_controller_impl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    build_adaptive_tuning_report as build_adaptive_tuning_report_impl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    build_adaptive_tuning_window_row as build_adaptive_tuning_window_row_impl,
)
from bijux_phylogenetics.bayesian.adaptive_tuning import (
    run_adaptive_tuned_metropolis_hastings_sampler as run_adaptive_tuned_metropolis_hastings_sampler_impl,
)


def test_bayesian_exports_adaptive_tuning_surface() -> None:
    assert AdaptiveTuningController is AdaptiveTuningControllerImpl
    assert AdaptiveTuningWindowRow is AdaptiveTuningWindowRowImpl
    assert AdaptiveTuningReport is AdaptiveTuningReportImpl
    assert (
        AdaptiveMetropolisHastingsRunReport is AdaptiveMetropolisHastingsRunReportImpl
    )
    assert build_adaptive_tuning_controller is build_adaptive_tuning_controller_impl
    assert build_adaptive_tuning_window_row is build_adaptive_tuning_window_row_impl
    assert build_adaptive_tuning_report is build_adaptive_tuning_report_impl
    assert (
        run_adaptive_tuned_metropolis_hastings_sampler
        is run_adaptive_tuned_metropolis_hastings_sampler_impl
    )
