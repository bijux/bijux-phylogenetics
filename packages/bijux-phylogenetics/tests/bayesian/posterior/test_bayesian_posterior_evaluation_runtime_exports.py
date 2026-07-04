from __future__ import annotations

import bijux_phylogenetics.bayesian as bayesian_api
from bijux_phylogenetics.bayesian import (
    BayesianPosteriorDecompositionRow,
    BayesianPosteriorEvaluationReport,
    evaluate_bayesian_posterior,
    summarize_bayesian_posterior_decomposition,
    write_bayesian_posterior_decomposition_table,
    write_bayesian_posterior_evaluation_artifacts,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    BayesianPosteriorDecompositionRow as BayesianPosteriorDecompositionRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    BayesianPosteriorEvaluationReport as BayesianPosteriorEvaluationReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    evaluate_bayesian_posterior as evaluate_bayesian_posterior_impl,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    summarize_bayesian_posterior_decomposition as summarize_bayesian_posterior_decomposition_impl,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    write_bayesian_posterior_decomposition_table as write_bayesian_posterior_decomposition_table_impl,
)
from bijux_phylogenetics.bayesian.posterior_evaluation import (
    write_bayesian_posterior_evaluation_artifacts as write_bayesian_posterior_evaluation_artifacts_impl,
)


def test_public_runtime_exports_bayesian_posterior_evaluation_surface() -> None:
    assert (
        bayesian_api.BayesianPosteriorDecompositionRow
        is BayesianPosteriorDecompositionRow
    )
    assert (
        bayesian_api.BayesianPosteriorDecompositionRow
        is BayesianPosteriorDecompositionRowImpl
    )
    assert (
        bayesian_api.BayesianPosteriorEvaluationReport
        is BayesianPosteriorEvaluationReport
    )
    assert (
        bayesian_api.BayesianPosteriorEvaluationReport
        is BayesianPosteriorEvaluationReportImpl
    )
    assert bayesian_api.evaluate_bayesian_posterior is evaluate_bayesian_posterior
    assert bayesian_api.evaluate_bayesian_posterior is evaluate_bayesian_posterior_impl
    assert (
        bayesian_api.summarize_bayesian_posterior_decomposition
        is summarize_bayesian_posterior_decomposition
    )
    assert (
        bayesian_api.summarize_bayesian_posterior_decomposition
        is summarize_bayesian_posterior_decomposition_impl
    )
    assert (
        bayesian_api.write_bayesian_posterior_decomposition_table
        is write_bayesian_posterior_decomposition_table
    )
    assert (
        bayesian_api.write_bayesian_posterior_decomposition_table
        is write_bayesian_posterior_decomposition_table_impl
    )
    assert (
        bayesian_api.write_bayesian_posterior_evaluation_artifacts
        is write_bayesian_posterior_evaluation_artifacts
    )
    assert (
        bayesian_api.write_bayesian_posterior_evaluation_artifacts
        is write_bayesian_posterior_evaluation_artifacts_impl
    )
