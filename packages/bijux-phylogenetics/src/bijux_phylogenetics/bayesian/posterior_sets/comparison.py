from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from bijux_phylogenetics.bayesian.beast.logs import assess_beast_convergence
from bijux_phylogenetics.bayesian.mrbayes import assess_mrbayes_convergence
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_tree_paths,
)

from .tree_sets import (
    BayesianRunTreeComparison,
    compare_bayesian_tree_sets,
    summarize_maximum_clade_credibility_tree,
    summarize_posterior_node_ages,
)


@dataclass(slots=True)
class BayesianParameterDifference:
    parameter: str
    left_mean: float
    right_mean: float
    mean_delta: float
    left_effective_sample_size: float
    right_effective_sample_size: float


@dataclass(slots=True)
class BayesianIndependentRunComparisonReport:
    left_tree_set_path: Path
    right_tree_set_path: Path
    left_trace_path: Path
    right_trace_path: Path
    trace_kind: str
    burnin_fraction: float
    tree_comparison: BayesianRunTreeComparison
    left_convergence: object
    right_convergence: object
    parameter_differences: list[BayesianParameterDifference]
    warnings: list[str]


@dataclass(slots=True)
class BayesianScenarioAgeDifference:
    clade: str
    left_mean_height: float
    right_mean_height: float
    delta: float


@dataclass(slots=True)
class BayesianPosteriorScenarioComparisonReport:
    comparison_axis: str
    left_label: str
    right_label: str
    left_tree_set_path: Path
    right_tree_set_path: Path
    burnin_fraction: float
    tree_comparison: BayesianRunTreeComparison
    age_differences: list[BayesianScenarioAgeDifference]
    warnings: list[str]


@dataclass(slots=True)
class BayesianMlTreeComparisonReport:
    ml_tree_path: Path
    posterior_tree_set_path: Path
    burnin_fraction: float
    mcc_tree_newick: str
    mcc_tree_index: int
    topology: object
    branch_lengths: object
    warnings: list[str]


def compare_independent_bayesian_runs(
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    *,
    left_trace_path: Path,
    right_trace_path: Path,
    trace_kind: str = "mrbayes",
    burnin_fraction: float = 0.25,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BayesianIndependentRunComparisonReport:
    """Compare two independent Bayesian runs across posterior trees and trace summaries."""
    tree_comparison = compare_bayesian_tree_sets(
        left_tree_set_path,
        right_tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    left_convergence = _assess_trace_convergence(
        left_trace_path,
        trace_kind=trace_kind,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    right_convergence = _assess_trace_convergence(
        right_trace_path,
        trace_kind=trace_kind,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    left_parameters = {
        str(row["parameter"]): row for row in left_convergence.parameter_summaries
    }
    right_parameters = {
        str(row["parameter"]): row for row in right_convergence.parameter_summaries
    }
    shared_parameters = sorted(left_parameters.keys() & right_parameters.keys())
    parameter_differences = [
        BayesianParameterDifference(
            parameter=parameter,
            left_mean=float(left_parameters[parameter]["mean"]),
            right_mean=float(right_parameters[parameter]["mean"]),
            mean_delta=round(
                float(right_parameters[parameter]["mean"])
                - float(left_parameters[parameter]["mean"]),
                6,
            ),
            left_effective_sample_size=float(
                left_parameters[parameter]["effective_sample_size"]
            ),
            right_effective_sample_size=float(
                right_parameters[parameter]["effective_sample_size"]
            ),
        )
        for parameter in shared_parameters
    ]
    warnings: list[str] = []
    if not tree_comparison.mcc_topology.topology_equal:
        warnings.append(
            "independent runs select different maximum clade credibility topologies"
        )
    if not left_convergence.converged or not right_convergence.converged:
        warnings.append(
            "one or more independent runs remain unconverged under the requested thresholds"
        )
    if any(abs(row.mean_delta) > 0.5 for row in parameter_differences):
        warnings.append(
            "one or more posterior parameter means differ materially across independent runs"
        )
    return BayesianIndependentRunComparisonReport(
        left_tree_set_path=left_tree_set_path,
        right_tree_set_path=right_tree_set_path,
        left_trace_path=left_trace_path,
        right_trace_path=right_trace_path,
        trace_kind=trace_kind,
        burnin_fraction=burnin_fraction,
        tree_comparison=tree_comparison,
        left_convergence=left_convergence,
        right_convergence=right_convergence,
        parameter_differences=parameter_differences,
        warnings=warnings,
    )


def compare_posterior_tree_sets_by_prior(
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    *,
    left_label: str = "prior_a",
    right_label: str = "prior_b",
    burnin_fraction: float = 0.25,
) -> BayesianPosteriorScenarioComparisonReport:
    """Compare posterior tree sets generated under alternative priors."""
    return _compare_posterior_tree_sets_by_scenario(
        left_tree_set_path,
        right_tree_set_path,
        comparison_axis="prior",
        left_label=left_label,
        right_label=right_label,
        burnin_fraction=burnin_fraction,
    )


def compare_posterior_tree_sets_by_clock(
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    *,
    left_label: str = "strict_clock",
    right_label: str = "relaxed_clock",
    burnin_fraction: float = 0.25,
) -> BayesianPosteriorScenarioComparisonReport:
    """Compare posterior tree sets generated under alternative clock models."""
    return _compare_posterior_tree_sets_by_scenario(
        left_tree_set_path,
        right_tree_set_path,
        comparison_axis="clock",
        left_label=left_label,
        right_label=right_label,
        burnin_fraction=burnin_fraction,
    )


def _compare_posterior_tree_sets_by_scenario(
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    *,
    comparison_axis: str,
    left_label: str,
    right_label: str,
    burnin_fraction: float,
) -> BayesianPosteriorScenarioComparisonReport:
    tree_comparison = compare_bayesian_tree_sets(
        left_tree_set_path,
        right_tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    left_ages = summarize_posterior_node_ages(
        left_tree_set_path, burnin_fraction=burnin_fraction
    )
    right_ages = summarize_posterior_node_ages(
        right_tree_set_path, burnin_fraction=burnin_fraction
    )
    left_rows = {row.clade: row for row in left_ages.rows}
    right_rows = {row.clade: row for row in right_ages.rows}
    shared_clades = sorted(left_rows.keys() & right_rows.keys())
    age_differences = [
        BayesianScenarioAgeDifference(
            clade=clade,
            left_mean_height=left_rows[clade].mean_height,
            right_mean_height=right_rows[clade].mean_height,
            delta=round(
                right_rows[clade].mean_height - left_rows[clade].mean_height, 15
            ),
        )
        for clade in shared_clades
    ]
    age_differences.sort(key=lambda row: (-abs(row.delta), row.clade))
    warnings: list[str] = []
    if not tree_comparison.mcc_topology.topology_equal:
        warnings.append(
            f"{comparison_axis} comparison changes the maximum clade credibility topology"
        )
    if any(abs(row.delta) > 0.25 for row in age_differences):
        warnings.append(
            f"{comparison_axis} comparison shifts one or more clade ages materially"
        )
    return BayesianPosteriorScenarioComparisonReport(
        comparison_axis=comparison_axis,
        left_label=left_label,
        right_label=right_label,
        left_tree_set_path=left_tree_set_path,
        right_tree_set_path=right_tree_set_path,
        burnin_fraction=burnin_fraction,
        tree_comparison=tree_comparison,
        age_differences=age_differences,
        warnings=warnings,
    )


def compare_ml_tree_to_bayesian_posterior(
    ml_tree_path: Path,
    posterior_tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> BayesianMlTreeComparisonReport:
    """Compare one ML summary tree against the MCC tree from a Bayesian posterior tree set."""
    _mcc_tree, mcc = summarize_maximum_clade_credibility_tree(
        posterior_tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    mcc_tree_path = Path(
        tempfile.mkstemp(prefix=f"{posterior_tree_set_path.stem}.mcc-", suffix=".nwk")[
            1
        ]
    )
    mcc_tree_path.write_text(mcc.mcc_newick + "\n", encoding="utf-8")
    try:
        topology = compare_tree_paths(ml_tree_path, mcc_tree_path)
        branch_lengths = compare_branch_lengths(ml_tree_path, mcc_tree_path)
        warnings: list[str] = []
        if not topology.topology_equal:
            warnings.append(
                "maximum-likelihood and Bayesian maximum clade credibility trees disagree in rooted topology"
            )
        if any(pair.delta not in {None, 0.0} for pair in branch_lengths.shared_splits):
            warnings.append(
                "shared clades differ in branch lengths between ML and Bayesian summaries"
            )
        if topology.same_taxa_different_rooting:
            warnings.append(
                "ML and Bayesian summaries share taxa and unrooted topology but disagree in rooting"
            )
        return BayesianMlTreeComparisonReport(
            ml_tree_path=ml_tree_path,
            posterior_tree_set_path=posterior_tree_set_path,
            burnin_fraction=burnin_fraction,
            mcc_tree_newick=mcc.mcc_newick,
            mcc_tree_index=mcc.selected_tree_index,
            topology=topology,
            branch_lengths=branch_lengths,
            warnings=warnings,
        )
    finally:
        mcc_tree_path.unlink(missing_ok=True)


def _assess_trace_convergence(
    path: Path,
    *,
    trace_kind: str,
    ess_threshold: float,
    mean_shift_threshold: float,
) -> object:
    normalized = trace_kind.lower()
    if normalized == "mrbayes":
        return assess_mrbayes_convergence(
            path,
            ess_threshold=ess_threshold,
            mean_shift_threshold=mean_shift_threshold,
        )
    if normalized == "beast":
        return assess_beast_convergence(
            path,
            ess_threshold=ess_threshold,
            mean_shift_threshold=mean_shift_threshold,
        )
    raise ValueError(f"unsupported trace_kind: {trace_kind}")
