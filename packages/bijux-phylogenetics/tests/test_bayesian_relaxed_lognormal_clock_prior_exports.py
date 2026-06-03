from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CLOCK_RATE_MODEL_FAMILIES,
    RELAXED_CLOCK_RATE_POLICIES,
    RelaxedLognormalClockBranchRow,
    RelaxedLognormalClockEvaluationReport,
    RelaxedLognormalClockModel,
    build_relaxed_lognormal_clock_model,
    evaluate_relaxed_lognormal_clock_tree_log_prior,
)
from bijux_phylogenetics.bayesian.clock_models import (
    CLOCK_RATE_MODEL_FAMILIES as CLOCK_RATE_MODEL_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.clock_models import (
    RELAXED_CLOCK_RATE_POLICIES as RELAXED_CLOCK_RATE_POLICIES_IMPL,
)
from bijux_phylogenetics.bayesian.clock_models import (
    RelaxedLognormalClockBranchRow as RelaxedLognormalClockBranchRowImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    RelaxedLognormalClockEvaluationReport as RelaxedLognormalClockEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    RelaxedLognormalClockModel as RelaxedLognormalClockModelImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    build_relaxed_lognormal_clock_model as build_relaxed_lognormal_clock_model_impl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    evaluate_relaxed_lognormal_clock_tree_log_prior as evaluate_relaxed_lognormal_clock_tree_log_prior_impl,
)


def test_bayesian_exports_relaxed_lognormal_clock_prior_surface() -> None:
    assert CLOCK_RATE_MODEL_FAMILIES == CLOCK_RATE_MODEL_FAMILIES_IMPL
    assert RELAXED_CLOCK_RATE_POLICIES == RELAXED_CLOCK_RATE_POLICIES_IMPL
    assert RelaxedLognormalClockBranchRow is RelaxedLognormalClockBranchRowImpl
    assert (
        RelaxedLognormalClockEvaluationReport
        is RelaxedLognormalClockEvaluationReportImpl
    )
    assert RelaxedLognormalClockModel is RelaxedLognormalClockModelImpl
    assert (
        build_relaxed_lognormal_clock_model is build_relaxed_lognormal_clock_model_impl
    )
    assert (
        evaluate_relaxed_lognormal_clock_tree_log_prior
        is evaluate_relaxed_lognormal_clock_tree_log_prior_impl
    )
