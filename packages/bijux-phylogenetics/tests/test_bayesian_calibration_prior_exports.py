from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    CALIBRATION_PRIOR_FAMILIES,
    CalibrationPriorDefinition,
    CalibrationPriorEvaluationReport,
    CalibrationPriorRow,
    evaluate_calibration_tree_log_prior,
    load_calibration_prior_definitions,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    CALIBRATION_PRIOR_FAMILIES as CALIBRATION_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    CalibrationPriorDefinition as CalibrationPriorDefinitionImpl,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    CalibrationPriorEvaluationReport as CalibrationPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    CalibrationPriorRow as CalibrationPriorRowImpl,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    evaluate_calibration_tree_log_prior as evaluate_calibration_tree_log_prior_impl,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    load_calibration_prior_definitions as load_calibration_prior_definitions_impl,
)


def test_bayesian_exports_calibration_prior_surface() -> None:
    assert CALIBRATION_PRIOR_FAMILIES == CALIBRATION_PRIOR_FAMILIES_IMPL
    assert CalibrationPriorDefinition is CalibrationPriorDefinitionImpl
    assert CalibrationPriorRow is CalibrationPriorRowImpl
    assert CalibrationPriorEvaluationReport is CalibrationPriorEvaluationReportImpl
    assert load_calibration_prior_definitions is load_calibration_prior_definitions_impl
    assert (
        evaluate_calibration_tree_log_prior is evaluate_calibration_tree_log_prior_impl
    )
