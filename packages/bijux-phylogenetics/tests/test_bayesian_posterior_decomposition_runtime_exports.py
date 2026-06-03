from __future__ import annotations

import bijux_phylogenetics.bayesian as bayesian_api
from bijux_phylogenetics.bayesian import (
    BeastPosteriorDecompositionReport,
    BeastPosteriorDecompositionRow,
    MrBayesPosteriorDecompositionReport,
    MrBayesPosteriorDecompositionRow,
    summarize_beast_posterior_decomposition,
    summarize_mrbayes_posterior_decomposition,
    write_beast_posterior_decomposition_table,
    write_mrbayes_posterior_decomposition_table,
)


def test_public_runtime_exports_bayesian_posterior_decomposition_surface() -> None:
    assert (
        bayesian_api.summarize_beast_posterior_decomposition
        is summarize_beast_posterior_decomposition
    )
    assert (
        bayesian_api.write_beast_posterior_decomposition_table
        is write_beast_posterior_decomposition_table
    )
    assert bayesian_api.BeastPosteriorDecompositionReport is (
        BeastPosteriorDecompositionReport
    )
    assert bayesian_api.BeastPosteriorDecompositionRow is (
        BeastPosteriorDecompositionRow
    )
    assert (
        bayesian_api.summarize_mrbayes_posterior_decomposition
        is summarize_mrbayes_posterior_decomposition
    )
    assert (
        bayesian_api.write_mrbayes_posterior_decomposition_table
        is write_mrbayes_posterior_decomposition_table
    )
    assert bayesian_api.MrBayesPosteriorDecompositionReport is (
        MrBayesPosteriorDecompositionReport
    )
    assert bayesian_api.MrBayesPosteriorDecompositionRow is (
        MrBayesPosteriorDecompositionRow
    )
