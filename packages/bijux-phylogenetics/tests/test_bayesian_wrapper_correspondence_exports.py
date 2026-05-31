from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES,
    BayesianWrapperCorrespondenceObservation,
    BayesianWrapperCorrespondenceReport,
    BayesianWrapperCorrespondenceSummaryRow,
    summarize_bayesian_wrapper_correspondence,
)
from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES as BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES_IMPL,
)
from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    BayesianWrapperCorrespondenceObservation as BayesianWrapperCorrespondenceObservationImpl,
)
from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    BayesianWrapperCorrespondenceReport as BayesianWrapperCorrespondenceReportImpl,
)
from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    BayesianWrapperCorrespondenceSummaryRow as BayesianWrapperCorrespondenceSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    summarize_bayesian_wrapper_correspondence as summarize_bayesian_wrapper_correspondence_impl,
)


def test_bayesian_exports_wrapper_correspondence_surface() -> None:
    assert (
        BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES
        == BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES_IMPL
    )
    assert (
        BayesianWrapperCorrespondenceObservation
        is BayesianWrapperCorrespondenceObservationImpl
    )
    assert BayesianWrapperCorrespondenceReport is BayesianWrapperCorrespondenceReportImpl
    assert (
        BayesianWrapperCorrespondenceSummaryRow
        is BayesianWrapperCorrespondenceSummaryRowImpl
    )
    assert (
        summarize_bayesian_wrapper_correspondence
        is summarize_bayesian_wrapper_correspondence_impl
    )
