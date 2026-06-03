from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    COALESCENT_TREE_PRIOR_FAMILIES,
    ConstantPopulationCoalescentIntervalRow,
    ConstantPopulationCoalescentPriorEvaluationReport,
    ConstantPopulationCoalescentPriorModel,
    build_constant_population_coalescent_tree_prior,
    evaluate_constant_population_coalescent_tree_log_prior,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    COALESCENT_TREE_PRIOR_FAMILIES as COALESCENT_TREE_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    ConstantPopulationCoalescentIntervalRow as ConstantPopulationCoalescentIntervalRowImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    ConstantPopulationCoalescentPriorEvaluationReport as ConstantPopulationCoalescentPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    ConstantPopulationCoalescentPriorModel as ConstantPopulationCoalescentPriorModelImpl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_constant_population_coalescent_tree_prior as build_constant_population_coalescent_tree_prior_impl,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    evaluate_constant_population_coalescent_tree_log_prior as evaluate_constant_population_coalescent_tree_log_prior_impl,
)


def test_bayesian_exports_constant_population_coalescent_prior_surface() -> None:
    assert COALESCENT_TREE_PRIOR_FAMILIES == COALESCENT_TREE_PRIOR_FAMILIES_IMPL
    assert (
        ConstantPopulationCoalescentPriorModel
        is ConstantPopulationCoalescentPriorModelImpl
    )
    assert (
        ConstantPopulationCoalescentIntervalRow
        is ConstantPopulationCoalescentIntervalRowImpl
    )
    assert (
        ConstantPopulationCoalescentPriorEvaluationReport
        is ConstantPopulationCoalescentPriorEvaluationReportImpl
    )
    assert (
        build_constant_population_coalescent_tree_prior
        is build_constant_population_coalescent_tree_prior_impl
    )
    assert (
        evaluate_constant_population_coalescent_tree_log_prior
        is evaluate_constant_population_coalescent_tree_log_prior_impl
    )
