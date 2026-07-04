from __future__ import annotations

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


def test_bayesian_exports_prior_protocol_surface() -> None:
    assert BayesianPriorModelProtocol is BayesianPriorModelProtocolImpl
    assert BayesianPriorModelSummary is BayesianPriorModelSummaryImpl
    assert summarize_bayesian_prior_model is summarize_bayesian_prior_model_impl
    assert (
        supports_bayesian_prior_model_protocol
        is supports_bayesian_prior_model_protocol_impl
    )
