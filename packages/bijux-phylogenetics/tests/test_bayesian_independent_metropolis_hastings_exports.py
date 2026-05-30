from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    IndependentMetropolisHastingsChainComparisonRow,
    IndependentMetropolisHastingsChainDefinition,
    IndependentMetropolisHastingsChainReport,
    IndependentMetropolisHastingsDiagnosticsReport,
    IndependentMetropolisHastingsRunReport,
    build_independent_metropolis_hastings_chain_comparison_row,
    build_independent_metropolis_hastings_chain_definition,
    build_independent_metropolis_hastings_chain_report,
    build_independent_metropolis_hastings_diagnostics_report,
    run_independent_metropolis_hastings_chains,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsChainComparisonRow as IndependentMetropolisHastingsChainComparisonRowImpl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsChainDefinition as IndependentMetropolisHastingsChainDefinitionImpl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsChainReport as IndependentMetropolisHastingsChainReportImpl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsDiagnosticsReport as IndependentMetropolisHastingsDiagnosticsReportImpl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsRunReport as IndependentMetropolisHastingsRunReportImpl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_chain_comparison_row as build_independent_metropolis_hastings_chain_comparison_row_impl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_chain_definition as build_independent_metropolis_hastings_chain_definition_impl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_chain_report as build_independent_metropolis_hastings_chain_report_impl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_diagnostics_report as build_independent_metropolis_hastings_diagnostics_report_impl,
)
from bijux_phylogenetics.bayesian.independent_chains import (
    run_independent_metropolis_hastings_chains as run_independent_metropolis_hastings_chains_impl,
)


def test_bayesian_exports_independent_metropolis_hastings_surface() -> None:
    assert (
        IndependentMetropolisHastingsChainComparisonRow
        is IndependentMetropolisHastingsChainComparisonRowImpl
    )
    assert (
        IndependentMetropolisHastingsChainDefinition
        is IndependentMetropolisHastingsChainDefinitionImpl
    )
    assert (
        IndependentMetropolisHastingsChainReport
        is IndependentMetropolisHastingsChainReportImpl
    )
    assert (
        IndependentMetropolisHastingsDiagnosticsReport
        is IndependentMetropolisHastingsDiagnosticsReportImpl
    )
    assert (
        IndependentMetropolisHastingsRunReport
        is IndependentMetropolisHastingsRunReportImpl
    )
    assert (
        build_independent_metropolis_hastings_chain_comparison_row
        is build_independent_metropolis_hastings_chain_comparison_row_impl
    )
    assert (
        build_independent_metropolis_hastings_chain_definition
        is build_independent_metropolis_hastings_chain_definition_impl
    )
    assert (
        build_independent_metropolis_hastings_chain_report
        is build_independent_metropolis_hastings_chain_report_impl
    )
    assert (
        build_independent_metropolis_hastings_diagnostics_report
        is build_independent_metropolis_hastings_diagnostics_report_impl
    )
    assert (
        run_independent_metropolis_hastings_chains
        is run_independent_metropolis_hastings_chains_impl
    )
