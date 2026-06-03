from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    LOCAL_CLOCK_RATE_MODEL_FAMILIES,
    LOCAL_CLOCK_TARGET_KINDS,
    LocalClockRateBranchRow,
    LocalClockRateEvaluationReport,
    LocalClockRateModel,
    LocalClockRateRegimeRow,
    LocalClockRegimeDefinition,
    build_local_clock_rate_model,
    evaluate_local_clock_tree_log_prior,
    load_local_clock_regime_definitions,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LOCAL_CLOCK_RATE_MODEL_FAMILIES as LOCAL_CLOCK_RATE_MODEL_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LOCAL_CLOCK_TARGET_KINDS as LOCAL_CLOCK_TARGET_KINDS_IMPL,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LocalClockRateBranchRow as LocalClockRateBranchRowImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LocalClockRateEvaluationReport as LocalClockRateEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LocalClockRateModel as LocalClockRateModelImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LocalClockRateRegimeRow as LocalClockRateRegimeRowImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    LocalClockRegimeDefinition as LocalClockRegimeDefinitionImpl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    build_local_clock_rate_model as build_local_clock_rate_model_impl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    evaluate_local_clock_tree_log_prior as evaluate_local_clock_tree_log_prior_impl,
)
from bijux_phylogenetics.bayesian.clock_models import (
    load_local_clock_regime_definitions as load_local_clock_regime_definitions_impl,
)


def test_bayesian_exports_local_clock_prior_surface() -> None:
    assert LOCAL_CLOCK_RATE_MODEL_FAMILIES == LOCAL_CLOCK_RATE_MODEL_FAMILIES_IMPL
    assert LOCAL_CLOCK_TARGET_KINDS == LOCAL_CLOCK_TARGET_KINDS_IMPL
    assert LocalClockRateBranchRow is LocalClockRateBranchRowImpl
    assert LocalClockRateEvaluationReport is LocalClockRateEvaluationReportImpl
    assert LocalClockRateModel is LocalClockRateModelImpl
    assert LocalClockRateRegimeRow is LocalClockRateRegimeRowImpl
    assert LocalClockRegimeDefinition is LocalClockRegimeDefinitionImpl
    assert build_local_clock_rate_model is build_local_clock_rate_model_impl
    assert (
        evaluate_local_clock_tree_log_prior is evaluate_local_clock_tree_log_prior_impl
    )
    assert (
        load_local_clock_regime_definitions is load_local_clock_regime_definitions_impl
    )
