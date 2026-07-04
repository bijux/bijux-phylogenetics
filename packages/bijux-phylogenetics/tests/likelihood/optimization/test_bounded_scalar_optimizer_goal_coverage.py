from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy

from bijux_phylogenetics.comparative import fit_continuous_evolutionary_mode
from bijux_phylogenetics.comparative.continuous import fit_ornstein_uhlenbeck_model
from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
)
from bijux_phylogenetics.phylo.likelihood.empirical import (
    optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    optimize_jc69_branch_lengths_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    optimize_k80_kappa_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


@dataclass(slots=True)
class Goal31CoverageCase:
    case_name: str
    lower_bound: float
    upper_bound: float
    starting_parameter_policy: str
    starting_parameter_value: float | None
    converged: bool
    failure_reason: str | None


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def fixture_path(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_bounded_scalar_optimizer_goal_coverage_matches_backlog_surface() -> None:
    shared_fixture = get_shared_geiger_continuous_fixture(
        "geiger_continuous_brownian_signal_twenty_four_taxa"
    )
    jc69_branch_report = optimize_jc69_branch_lengths_from_alignment(
        fixture("trees", "jc69_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_branch_optimization_alignment_2_taxa.fasta"),
    )
    k80_report = optimize_k80_kappa_from_alignment(
        fixture("trees", "k80_kappa_optimization_tree_2_taxa.nwk"),
        fixture("alignments", "k80_kappa_optimization_alignment_2_taxa.fasta"),
    )
    ou_report = fit_ornstein_uhlenbeck_model(
        fixture_path("example_tree.nwk"),
        fixture_path("example_traits_comparative.tsv"),
        trait="response",
        taxon_column="taxon",
    )
    empirical_invariant_report = optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments",
            "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        initial_invariant_proportion=0.1,
    )
    pgls_report = run_pgls(
        fixture_path("example_tree.nwk"),
        fixture_path("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        taxon_column="taxon",
        lambda_value="estimate",
    )
    pagel_kappa_report = fit_continuous_evolutionary_mode(
        shared_fixture.tree_path,
        shared_fixture.traits_path,
        trait=shared_fixture.trait_name,
        mode="pagel-kappa",
        taxon_column=shared_fixture.taxon_column,
    )
    pagel_delta_report = fit_continuous_evolutionary_mode(
        shared_fixture.tree_path,
        shared_fixture.traits_path,
        trait=shared_fixture.trait_name,
        mode="pagel-delta",
        taxon_column=shared_fixture.taxon_column,
    )

    first_branch_diagnostics = next(
        iter(jc69_branch_report.parameter_diagnostics.values())
    )
    cases = [
        Goal31CoverageCase(
            case_name="branch length",
            lower_bound=first_branch_diagnostics.lower_bound,
            upper_bound=first_branch_diagnostics.upper_bound,
            starting_parameter_policy=first_branch_diagnostics.starting_parameter_policy,
            starting_parameter_value=first_branch_diagnostics.starting_parameter_value,
            converged=first_branch_diagnostics.converged,
            failure_reason=first_branch_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="kappa",
            lower_bound=k80_report.optimizer_diagnostics.lower_bound,
            upper_bound=k80_report.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=k80_report.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=k80_report.optimizer_diagnostics.starting_parameter_value,
            converged=k80_report.optimizer_diagnostics.converged,
            failure_reason=k80_report.optimizer_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="alpha",
            lower_bound=ou_report.optimizer_diagnostics.lower_bound,
            upper_bound=ou_report.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=ou_report.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=ou_report.optimizer_diagnostics.starting_parameter_value,
            converged=ou_report.optimizer_diagnostics.converged,
            failure_reason=ou_report.optimizer_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="invariant proportion",
            lower_bound=empirical_invariant_report.optimizer_diagnostics.lower_bound,
            upper_bound=empirical_invariant_report.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=empirical_invariant_report.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=empirical_invariant_report.optimizer_diagnostics.starting_parameter_value,
            converged=empirical_invariant_report.optimizer_diagnostics.converged,
            failure_reason=empirical_invariant_report.optimizer_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="lambda",
            lower_bound=pgls_report.lambda_fit.optimizer_diagnostics.lower_bound,
            upper_bound=pgls_report.lambda_fit.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=pgls_report.lambda_fit.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=pgls_report.lambda_fit.optimizer_diagnostics.starting_parameter_value,
            converged=pgls_report.lambda_fit.optimizer_diagnostics.converged,
            failure_reason=pgls_report.lambda_fit.optimizer_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="kappa-transform",
            lower_bound=pagel_kappa_report.optimizer_diagnostics.lower_bound,
            upper_bound=pagel_kappa_report.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=pagel_kappa_report.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=pagel_kappa_report.optimizer_diagnostics.starting_parameter_value,
            converged=pagel_kappa_report.optimizer_diagnostics.converged,
            failure_reason=pagel_kappa_report.optimizer_diagnostics.failure_reason,
        ),
        Goal31CoverageCase(
            case_name="delta-transform",
            lower_bound=pagel_delta_report.optimizer_diagnostics.lower_bound,
            upper_bound=pagel_delta_report.optimizer_diagnostics.upper_bound,
            starting_parameter_policy=pagel_delta_report.optimizer_diagnostics.starting_parameter_policy,
            starting_parameter_value=pagel_delta_report.optimizer_diagnostics.starting_parameter_value,
            converged=pagel_delta_report.optimizer_diagnostics.converged,
            failure_reason=pagel_delta_report.optimizer_diagnostics.failure_reason,
        ),
    ]

    assert [case.case_name for case in cases] == [
        "branch length",
        "kappa",
        "alpha",
        "invariant proportion",
        "lambda",
        "kappa-transform",
        "delta-transform",
    ]
    for case in cases:
        assert case.lower_bound < case.upper_bound
        assert case.starting_parameter_policy
        if case.starting_parameter_value is not None:
            assert case.lower_bound <= case.starting_parameter_value <= case.upper_bound
        assert case.converged is True
        assert case.failure_reason is None


def _compact_polar_rate_matrix() -> numpy.ndarray:
    state_order = tuple("ACDEFGHIKLMNPQRSTVWY")
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for left_state, right_state, rate in (
        ("A", "C", 0.45),
        ("C", "D", 0.35),
        ("D", "E", 0.55),
        ("A", "E", 0.20),
    ):
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())
