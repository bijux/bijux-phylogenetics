from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CLOCK_RATE_MODEL_FAMILIES,
    StrictClockRateBranchRow,
    StrictClockRateModel,
    StrictClockRateModelEvaluationReport,
    build_strict_clock_rate_model,
    evaluate_strict_clock_tree_log_prior,
)
from bijux_phylogenetics.bayesian.clock_models import (
    CLOCK_RATE_MODEL_FAMILIES as CLOCK_RATE_MODEL_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.clock_models import (
    StrictClockRateBranchRow as StrictClockRateBranchRowImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    StrictClockRateModel as StrictClockRateModelImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    StrictClockRateModelEvaluationReport as StrictClockRateModelEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    build_strict_clock_rate_model as build_strict_clock_rate_model_impl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    evaluate_strict_clock_tree_log_prior as evaluate_strict_clock_tree_log_prior_impl,
)


def test_bayesian_exports_strict_clock_rate_model_surface() -> None:
    assert CLOCK_RATE_MODEL_FAMILIES == CLOCK_RATE_MODEL_FAMILIES_IMPL
    assert StrictClockRateBranchRow is StrictClockRateBranchRowImpl
    assert StrictClockRateModel is StrictClockRateModelImpl
    assert (
        StrictClockRateModelEvaluationReport is StrictClockRateModelEvaluationReportImpl
    )
    assert build_strict_clock_rate_model is build_strict_clock_rate_model_impl
    assert (
        evaluate_strict_clock_tree_log_prior
        is evaluate_strict_clock_tree_log_prior_impl
    )
