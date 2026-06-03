from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BAYESIAN_BURNIN_POLICY_NAMES,
    BayesianRunBurninPolicy,
    BayesianRunManifest,
    BayesianRunManifestReplayReport,
    BayesianRunPriorRow,
    build_bayesian_run_burnin_policy,
    build_bayesian_run_manifest,
    build_fixed_topology_dna_run_manifest,
    list_metropolis_hastings_retained_sample_ids,
    load_bayesian_run_manifest,
    replay_fixed_topology_dna_run_manifest,
    write_bayesian_run_manifest,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    BAYESIAN_BURNIN_POLICY_NAMES as BAYESIAN_BURNIN_POLICY_NAMES_IMPL,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    BayesianRunBurninPolicy as BayesianRunBurninPolicyImpl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    BayesianRunManifest as BayesianRunManifestImpl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    BayesianRunManifestReplayReport as BayesianRunManifestReplayReportImpl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    BayesianRunPriorRow as BayesianRunPriorRowImpl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    build_bayesian_run_burnin_policy as build_bayesian_run_burnin_policy_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    build_bayesian_run_manifest as build_bayesian_run_manifest_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    build_fixed_topology_dna_run_manifest as build_fixed_topology_dna_run_manifest_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    list_metropolis_hastings_retained_sample_ids as list_metropolis_hastings_retained_sample_ids_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    load_bayesian_run_manifest as load_bayesian_run_manifest_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    replay_fixed_topology_dna_run_manifest as replay_fixed_topology_dna_run_manifest_impl,
)
from bijux_phylogenetics.bayesian.run_manifest import (
    write_bayesian_run_manifest as write_bayesian_run_manifest_impl,
)


def test_bayesian_exports_run_manifest_surface() -> None:
    assert BAYESIAN_BURNIN_POLICY_NAMES == BAYESIAN_BURNIN_POLICY_NAMES_IMPL
    assert BayesianRunBurninPolicy is BayesianRunBurninPolicyImpl
    assert BayesianRunManifest is BayesianRunManifestImpl
    assert BayesianRunManifestReplayReport is BayesianRunManifestReplayReportImpl
    assert BayesianRunPriorRow is BayesianRunPriorRowImpl
    assert build_bayesian_run_burnin_policy is build_bayesian_run_burnin_policy_impl
    assert build_bayesian_run_manifest is build_bayesian_run_manifest_impl
    assert (
        build_fixed_topology_dna_run_manifest
        is build_fixed_topology_dna_run_manifest_impl
    )
    assert (
        list_metropolis_hastings_retained_sample_ids
        is list_metropolis_hastings_retained_sample_ids_impl
    )
    assert load_bayesian_run_manifest is load_bayesian_run_manifest_impl
    assert (
        replay_fixed_topology_dna_run_manifest
        is replay_fixed_topology_dna_run_manifest_impl
    )
    assert write_bayesian_run_manifest is write_bayesian_run_manifest_impl
