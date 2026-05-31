from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    BayesianPosteriorTreeSample,
    BayesianPosteriorTreeSampleArchive,
    build_bayesian_posterior_tree_sample,
    build_bayesian_posterior_tree_sample_archive,
    build_metropolis_hastings_posterior_tree_sample_archive,
    infer_bayesian_model_id,
    load_bayesian_posterior_tree_sample_archive,
    write_bayesian_posterior_tree_sample_archive,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    BayesianPosteriorTreeSample as BayesianPosteriorTreeSampleImpl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    BayesianPosteriorTreeSampleArchive as BayesianPosteriorTreeSampleArchiveImpl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    build_bayesian_posterior_tree_sample as build_bayesian_posterior_tree_sample_impl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    build_bayesian_posterior_tree_sample_archive as build_bayesian_posterior_tree_sample_archive_impl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    build_metropolis_hastings_posterior_tree_sample_archive as build_metropolis_hastings_posterior_tree_sample_archive_impl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    infer_bayesian_model_id as infer_bayesian_model_id_impl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    load_bayesian_posterior_tree_sample_archive as load_bayesian_posterior_tree_sample_archive_impl,
)
from bijux_phylogenetics.bayesian.posterior_tree_samples import (
    write_bayesian_posterior_tree_sample_archive as write_bayesian_posterior_tree_sample_archive_impl,
)


def test_bayesian_exports_posterior_tree_sample_archive_surface() -> None:
    assert BayesianPosteriorTreeSample is BayesianPosteriorTreeSampleImpl
    assert BayesianPosteriorTreeSampleArchive is BayesianPosteriorTreeSampleArchiveImpl
    assert (
        build_bayesian_posterior_tree_sample
        is build_bayesian_posterior_tree_sample_impl
    )
    assert (
        build_bayesian_posterior_tree_sample_archive
        is build_bayesian_posterior_tree_sample_archive_impl
    )
    assert (
        build_metropolis_hastings_posterior_tree_sample_archive
        is build_metropolis_hastings_posterior_tree_sample_archive_impl
    )
    assert infer_bayesian_model_id is infer_bayesian_model_id_impl
    assert (
        load_bayesian_posterior_tree_sample_archive
        is load_bayesian_posterior_tree_sample_archive_impl
    )
    assert (
        write_bayesian_posterior_tree_sample_archive
        is write_bayesian_posterior_tree_sample_archive_impl
    )
