from __future__ import annotations

import bijux_phylogenetics.bayesian as bayesian_api
from bijux_phylogenetics.bayesian import (
    BayesianPriorModelProtocol,
    BayesianPriorModelSummary,
    summarize_bayesian_prior_model,
    supports_bayesian_prior_model_protocol,
)
from bijux_phylogenetics.bayesian.prior_protocol import (
    BayesianPriorModelProtocol as BayesianPriorModelProtocolImpl,
)
from bijux_phylogenetics.bayesian.prior_protocol import (
    BayesianPriorModelSummary as BayesianPriorModelSummaryImpl,
)
from bijux_phylogenetics.bayesian.prior_protocol import (
    summarize_bayesian_prior_model as summarize_bayesian_prior_model_impl,
)
from bijux_phylogenetics.bayesian.prior_protocol import (
    supports_bayesian_prior_model_protocol as supports_bayesian_prior_model_protocol_impl,
)


def test_public_runtime_exports_bayesian_prior_protocol_surface() -> None:
    assert bayesian_api.BayesianPriorModelProtocol is BayesianPriorModelProtocol
    assert bayesian_api.BayesianPriorModelProtocol is BayesianPriorModelProtocolImpl
    assert bayesian_api.BayesianPriorModelSummary is BayesianPriorModelSummary
    assert bayesian_api.BayesianPriorModelSummary is BayesianPriorModelSummaryImpl
    assert bayesian_api.summarize_bayesian_prior_model is summarize_bayesian_prior_model
    assert (
        bayesian_api.summarize_bayesian_prior_model
        is summarize_bayesian_prior_model_impl
    )
    assert (
        bayesian_api.supports_bayesian_prior_model_protocol
        is supports_bayesian_prior_model_protocol
    )
    assert (
        bayesian_api.supports_bayesian_prior_model_protocol
        is supports_bayesian_prior_model_protocol_impl
    )
