from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BRANCH_LENGTH_PRIOR_FAMILIES,
    BranchLengthPriorBranchRow,
    BranchLengthPriorEvaluationReport,
    BranchLengthPriorModel,
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
    build_gamma_branch_length_prior,
    build_lognormal_branch_length_prior,
    evaluate_branch_length_log_prior,
    evaluate_tree_branch_length_log_prior,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    BRANCH_LENGTH_PRIOR_FAMILIES as BRANCH_LENGTH_PRIOR_FAMILIES_IMPL,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorBranchRow as BranchLengthPriorBranchRowImpl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorEvaluationReport as BranchLengthPriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    BranchLengthPriorModel as BranchLengthPriorModelImpl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior as build_exponential_branch_length_prior_impl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_fixed_branch_length_prior as build_fixed_branch_length_prior_impl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_gamma_branch_length_prior as build_gamma_branch_length_prior_impl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_lognormal_branch_length_prior as build_lognormal_branch_length_prior_impl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    evaluate_branch_length_log_prior as evaluate_branch_length_log_prior_impl,
)
from bijux_phylogenetics.bayesian.branch_length_priors import (
    evaluate_tree_branch_length_log_prior as evaluate_tree_branch_length_log_prior_impl,
)


def test_bayesian_exports_branch_length_prior_surface() -> None:
    assert BRANCH_LENGTH_PRIOR_FAMILIES == BRANCH_LENGTH_PRIOR_FAMILIES_IMPL
    assert BranchLengthPriorModel is BranchLengthPriorModelImpl
    assert BranchLengthPriorBranchRow is BranchLengthPriorBranchRowImpl
    assert BranchLengthPriorEvaluationReport is BranchLengthPriorEvaluationReportImpl
    assert (
        build_exponential_branch_length_prior
        is build_exponential_branch_length_prior_impl
    )
    assert build_gamma_branch_length_prior is build_gamma_branch_length_prior_impl
    assert (
        build_lognormal_branch_length_prior is build_lognormal_branch_length_prior_impl
    )
    assert build_fixed_branch_length_prior is build_fixed_branch_length_prior_impl
    assert evaluate_branch_length_log_prior is evaluate_branch_length_log_prior_impl
    assert (
        evaluate_tree_branch_length_log_prior
        is evaluate_tree_branch_length_log_prior_impl
    )
