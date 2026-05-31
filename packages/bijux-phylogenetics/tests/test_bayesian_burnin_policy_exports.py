from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES,
    BurninSampleRow,
    IndependentMetropolisHastingsBurninReport,
    IndependentMetropolisHastingsChainBurninReport,
    MetropolisHastingsBurninDiagnosticCandidate,
    MetropolisHastingsBurninDiagnosticReport,
    MetropolisHastingsBurninPolicy,
    MetropolisHastingsBurninReport,
    apply_independent_metropolis_hastings_burnin_policy,
    apply_metropolis_hastings_burnin_policy,
    build_metropolis_hastings_burnin_policy,
    diagnose_metropolis_hastings_burnin,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES as METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES_IMPL,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    BurninSampleRow as BurninSampleRowImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    IndependentMetropolisHastingsBurninReport as IndependentMetropolisHastingsBurninReportImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    IndependentMetropolisHastingsChainBurninReport as IndependentMetropolisHastingsChainBurninReportImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    MetropolisHastingsBurninDiagnosticCandidate as MetropolisHastingsBurninDiagnosticCandidateImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    MetropolisHastingsBurninDiagnosticReport as MetropolisHastingsBurninDiagnosticReportImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    MetropolisHastingsBurninPolicy as MetropolisHastingsBurninPolicyImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    MetropolisHastingsBurninReport as MetropolisHastingsBurninReportImpl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    apply_independent_metropolis_hastings_burnin_policy as apply_independent_metropolis_hastings_burnin_policy_impl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    apply_metropolis_hastings_burnin_policy as apply_metropolis_hastings_burnin_policy_impl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    build_metropolis_hastings_burnin_policy as build_metropolis_hastings_burnin_policy_impl,
)
from bijux_phylogenetics.bayesian.burnin_policies import (
    diagnose_metropolis_hastings_burnin as diagnose_metropolis_hastings_burnin_impl,
)


def test_bayesian_exports_burnin_policy_surface() -> None:
    assert (
        METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES
        == METROPOLIS_HASTINGS_BURNIN_POLICY_NAMES_IMPL
    )
    assert BurninSampleRow is BurninSampleRowImpl
    assert (
        IndependentMetropolisHastingsBurninReport
        is IndependentMetropolisHastingsBurninReportImpl
    )
    assert (
        IndependentMetropolisHastingsChainBurninReport
        is IndependentMetropolisHastingsChainBurninReportImpl
    )
    assert (
        MetropolisHastingsBurninDiagnosticCandidate
        is MetropolisHastingsBurninDiagnosticCandidateImpl
    )
    assert (
        MetropolisHastingsBurninDiagnosticReport
        is MetropolisHastingsBurninDiagnosticReportImpl
    )
    assert MetropolisHastingsBurninPolicy is MetropolisHastingsBurninPolicyImpl
    assert MetropolisHastingsBurninReport is MetropolisHastingsBurninReportImpl
    assert (
        apply_independent_metropolis_hastings_burnin_policy
        is apply_independent_metropolis_hastings_burnin_policy_impl
    )
    assert (
        apply_metropolis_hastings_burnin_policy
        is apply_metropolis_hastings_burnin_policy_impl
    )
    assert (
        build_metropolis_hastings_burnin_policy
        is build_metropolis_hastings_burnin_policy_impl
    )
    assert (
        diagnose_metropolis_hastings_burnin is diagnose_metropolis_hastings_burnin_impl
    )
