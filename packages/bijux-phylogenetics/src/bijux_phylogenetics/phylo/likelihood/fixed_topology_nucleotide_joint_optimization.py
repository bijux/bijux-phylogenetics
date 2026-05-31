from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    optimize_fixed_topology_nucleotide_branch_lengths,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    evaluate_gtr_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    evaluate_hky85_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    BranchLengthOptimizationRow,
    FixedTopologyNucleotideJointOptimizationReport,
    JointNucleotideOptimizationUpdateRow,
    SubstitutionParameterOptimizationRow,
)
from bijux_phylogenetics.phylo.likelihood.substitution_parameters import (
    optimize_nucleotide_substitution_parameters,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_FIXED_TOPOLOGY_JOINT_MODELS = frozenset({"k80", "hky85", "gtr"})
_BOUNDARY_TOLERANCE = 1e-9


def validate_fixed_topology_nucleotide_joint_optimization_model(model_name: str) -> str:
    """Validate one nucleotide model with a real joint branch-and-model optimization surface."""
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _SUPPORTED_FIXED_TOPOLOGY_JOINT_MODELS:
        raise ValueError(
            "fixed-topology nucleotide joint optimization model must be one of "
            + ", ".join(sorted(_SUPPORTED_FIXED_TOPOLOGY_JOINT_MODELS))
        )
    return normalized_model_name


def optimize_fixed_topology_nucleotide_branches_and_model(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_kappa: float | None = None,
    lower_kappa_bound: float | None = None,
    upper_kappa_bound: float | None = None,
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    lower_exchangeability_bound: float | None = None,
    upper_exchangeability_bound: float | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_joint_passes: int = 8,
    max_branch_coordinate_passes: int = 12,
    max_model_coordinate_passes: int = 24,
) -> FixedTopologyNucleotideJointOptimizationReport:
    """Jointly optimize branch lengths and substitution parameters on one fixed nucleotide topology."""
    normalized_model_name = validate_fixed_topology_nucleotide_joint_optimization_model(
        model_name
    )
    if improvement_tolerance < 0.0:
        raise ValueError(
            "fixed-topology nucleotide joint optimization improvement_tolerance must be nonnegative"
        )
    if max_joint_passes < 1:
        raise ValueError(
            "fixed-topology nucleotide joint optimization requires at least one joint pass"
        )
    if max_branch_coordinate_passes < 1:
        raise ValueError(
            "fixed-topology nucleotide joint optimization requires at least one branch coordinate pass"
        )
    if max_model_coordinate_passes < 1:
        raise ValueError(
            "fixed-topology nucleotide joint optimization requires at least one model coordinate pass"
        )

    initial_tree = tree.copy().refresh()
    initial_tree_newick = dumps_newick(initial_tree)
    current_tree = initial_tree.copy().refresh()
    provided_base_frequencies = base_frequencies
    current_base_frequencies = base_frequencies
    current_kappa = initial_kappa
    current_exchangeabilities = (
        {
            "AC": 1.0,
            "AG": 1.0,
            "AT": 1.0,
            "CG": 1.0,
            "CT": 1.0,
            "GT": 1.0,
        }
        if normalized_model_name == "gtr" and initial_exchangeabilities is None
        else initial_exchangeabilities
    )
    initial_log_likelihood = _evaluate_joint_state_log_likelihood(
        current_tree,
        records,
        model_name=normalized_model_name,
        base_frequencies=current_base_frequencies,
        kappa=current_kappa,
        exchangeabilities=current_exchangeabilities,
    )
    current_log_likelihood = initial_log_likelihood

    update_rows: list[JointNucleotideOptimizationUpdateRow] = []
    warnings: list[str] = []
    total_function_evaluation_count = 0
    executed_joint_pass_count = 0
    schedule_converged = False
    inner_converged = True
    latest_branch_report = None
    latest_model_report = None

    for joint_pass_index in range(1, max_joint_passes + 1):
        executed_joint_pass_count = joint_pass_index
        branch_report = optimize_fixed_topology_nucleotide_branch_lengths(
            current_tree,
            records,
            model_name=normalized_model_name,
            kappa=current_kappa,
            base_frequencies=current_base_frequencies,
            exchangeabilities=current_exchangeabilities,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_branch_coordinate_passes,
        )
        latest_branch_report = branch_report
        total_function_evaluation_count += branch_report.function_evaluation_count
        branch_delta = (
            branch_report.optimized_log_likelihood - current_log_likelihood
        )
        current_tree = loads_newick(branch_report.optimized_tree_newick).refresh()
        current_log_likelihood = branch_report.optimized_log_likelihood
        update_rows.append(
            JointNucleotideOptimizationUpdateRow(
                joint_pass_index=joint_pass_index,
                update_kind="branch-lengths",
                starting_log_likelihood=branch_report.initial_log_likelihood,
                optimized_log_likelihood=branch_report.optimized_log_likelihood,
                log_likelihood_delta=branch_report.optimized_log_likelihood
                - branch_report.initial_log_likelihood,
                function_evaluation_count=branch_report.function_evaluation_count,
                optimization_pass_count=branch_report.optimization_pass_count,
                converged=branch_report.converged,
                optimized_branch_count=len(branch_report.branches),
                optimized_branch_ids=[
                    row.branch_id
                    for row in branch_report.branches
                ],
                updated_parameter_names=[],
            )
        )
        if not branch_report.converged:
            inner_converged = False

        model_report = optimize_nucleotide_substitution_parameters(
            current_tree,
            records,
            model_name=normalized_model_name,
            base_frequencies=provided_base_frequencies,
            initial_kappa=current_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
            initial_exchangeabilities=current_exchangeabilities,
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
            max_coordinate_passes=max_model_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
        latest_model_report = model_report
        total_function_evaluation_count += model_report.function_evaluation_count
        model_delta = model_report.optimized_log_likelihood - current_log_likelihood
        current_log_likelihood = model_report.optimized_log_likelihood
        current_base_frequencies = _resolved_base_frequencies(model_report)
        current_kappa = _resolved_kappa(model_report)
        current_exchangeabilities = _resolved_exchangeabilities(model_report)
        update_rows.append(
            JointNucleotideOptimizationUpdateRow(
                joint_pass_index=joint_pass_index,
                update_kind="substitution-parameters",
                starting_log_likelihood=model_report.initial_log_likelihood,
                optimized_log_likelihood=model_report.optimized_log_likelihood,
                log_likelihood_delta=model_report.optimized_log_likelihood
                - model_report.initial_log_likelihood,
                function_evaluation_count=model_report.function_evaluation_count,
                optimization_pass_count=model_report.optimization_pass_count,
                converged=model_report.converged,
                optimized_branch_count=0,
                optimized_branch_ids=[],
                updated_parameter_names=[
                    row.parameter_name
                    for row in model_report.parameter_rows
                ],
            )
        )
        if not model_report.converged:
            inner_converged = False
        for warning in model_report.warnings:
            if warning not in warnings:
                warnings.append(warning)

        if (
            branch_delta <= improvement_tolerance
            and model_delta <= improvement_tolerance
        ):
            schedule_converged = True
            break

    if latest_branch_report is None or latest_model_report is None:
        raise AssertionError("joint optimization did not execute any branch-and-model pass")

    if not schedule_converged:
        convergence_reason = "max-joint-passes-exhausted"
    elif not inner_converged:
        convergence_reason = "inner-optimizer-did-not-converge"
    else:
        convergence_reason = "joint-schedule-converged"

    final_parameter_rows = _build_final_parameter_rows(
        model_name=normalized_model_name,
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
        final_kappa=current_kappa,
        initial_exchangeabilities=initial_exchangeabilities,
        lower_exchangeability_bound=lower_exchangeability_bound,
        upper_exchangeability_bound=upper_exchangeability_bound,
        final_exchangeabilities=current_exchangeabilities,
    )
    return FixedTopologyNucleotideJointOptimizationReport(
        model_name=latest_model_report.model_name,
        taxa=latest_branch_report.taxa,
        site_count=latest_branch_report.site_count,
        pattern_count=latest_branch_report.pattern_count,
        branch_count=latest_branch_report.branch_count,
        initial_tree_newick=initial_tree_newick,
        optimized_tree_newick=dumps_newick(current_tree),
        state_count=latest_branch_report.state_count,
        observation_policy=latest_branch_report.observation_policy,
        root_prior_source=latest_branch_report.root_prior_source,
        parameter_count=latest_model_report.parameter_count,
        base_frequency_source=latest_model_report.base_frequency_source,
        base_frequency_a=latest_model_report.base_frequency_a,
        base_frequency_c=latest_model_report.base_frequency_c,
        base_frequency_g=latest_model_report.base_frequency_g,
        base_frequency_t=latest_model_report.base_frequency_t,
        fixed_parameter_values=dict(latest_model_report.fixed_parameter_values),
        parameter_rows=final_parameter_rows,
        branch_rows=_build_final_branch_rows(initial_tree, current_tree),
        initial_log_likelihood=initial_log_likelihood,
        optimized_log_likelihood=current_log_likelihood,
        function_evaluation_count=total_function_evaluation_count,
        joint_optimization_pass_count=executed_joint_pass_count,
        converged=(schedule_converged and inner_converged),
        convergence_reason=convergence_reason,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        update_rows=update_rows,
        warnings=warnings,
    )


def optimize_fixed_topology_nucleotide_branches_and_model_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_kappa: float | None = None,
    lower_kappa_bound: float | None = None,
    upper_kappa_bound: float | None = None,
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    lower_exchangeability_bound: float | None = None,
    upper_exchangeability_bound: float | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_joint_passes: int = 8,
    max_branch_coordinate_passes: int = 12,
    max_model_coordinate_passes: int = 24,
) -> FixedTopologyNucleotideJointOptimizationReport:
    """Jointly optimize one fixed-topology nucleotide model from tree and alignment paths."""
    return optimize_fixed_topology_nucleotide_branches_and_model(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        base_frequencies=base_frequencies,
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
        initial_exchangeabilities=initial_exchangeabilities,
        lower_exchangeability_bound=lower_exchangeability_bound,
        upper_exchangeability_bound=upper_exchangeability_bound,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_joint_passes=max_joint_passes,
        max_branch_coordinate_passes=max_branch_coordinate_passes,
        max_model_coordinate_passes=max_model_coordinate_passes,
    )


def _evaluate_joint_state_log_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    kappa: float | None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ),
) -> float:
    if model_name == "k80":
        report = evaluate_k80_tree_likelihood(
            tree,
            records,
            kappa=(1.0 if kappa is None else kappa),
        )
        return report.log_likelihood
    if model_name == "hky85":
        report = evaluate_hky85_tree_likelihood(
            tree,
            records,
            base_frequencies=base_frequencies,
            kappa=(1.0 if kappa is None else kappa),
        )
        return report.log_likelihood
    report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        base_frequencies=base_frequencies,
        exchangeabilities=(
            {
                "AC": 1.0,
                "AG": 1.0,
                "AT": 1.0,
                "CG": 1.0,
                "CT": 1.0,
                "GT": 1.0,
            }
            if exchangeabilities is None
            else exchangeabilities
        ),
    )
    return report.log_likelihood


def _resolved_base_frequencies(
    optimization_report,
) -> dict[str, float] | None:
    if optimization_report.base_frequency_a is None:
        return None
    return {
        "A": optimization_report.base_frequency_a,
        "C": optimization_report.base_frequency_c or 0.0,
        "G": optimization_report.base_frequency_g or 0.0,
        "T": optimization_report.base_frequency_t or 0.0,
    }


def _resolved_kappa(
    optimization_report,
) -> float | None:
    for row in optimization_report.parameter_rows:
        if row.parameter_name == "kappa":
            return row.optimized_value
    return None


def _resolved_exchangeabilities(
    optimization_report,
) -> dict[str, float] | None:
    if optimization_report.model_name != "GTR":
        return None
    optimized_by_name = {
        row.parameter_name: row.optimized_value
        for row in optimization_report.parameter_rows
    }
    return {
        "AC": optimization_report.fixed_parameter_values.get("AC", 1.0),
        "AG": optimized_by_name["AG"],
        "AT": optimized_by_name["AT"],
        "CG": optimized_by_name["CG"],
        "CT": optimized_by_name["CT"],
        "GT": optimized_by_name["GT"],
    }


def _build_final_branch_rows(
    initial_tree: PhyloTree,
    optimized_tree: PhyloTree,
) -> list[BranchLengthOptimizationRow]:
    initial_tree = initial_tree.copy().refresh()
    optimized_tree = optimized_tree.copy().refresh()
    optimized_lengths = {
        node.node_id or "": float(node.branch_length or 0.0)
        for _parent, node in optimized_tree.iter_edges()
    }
    return [
        BranchLengthOptimizationRow(
            branch_id=node.node_id or "",
            child_name=node.name,
            descendant_taxa=node.descendant_taxa,
            initial_branch_length=float(node.branch_length or 0.0),
            optimized_branch_length=optimized_lengths[node.node_id or ""],
        )
        for _parent, node in initial_tree.iter_edges()
    ]


def _build_final_parameter_rows(
    *,
    model_name: str,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
    final_kappa: float | None,
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ),
    lower_exchangeability_bound: float | None,
    upper_exchangeability_bound: float | None,
    final_exchangeabilities: dict[str, float] | None,
) -> list[SubstitutionParameterOptimizationRow]:
    if model_name in {"k80", "hky85"}:
        resolved_initial_kappa = 1.0 if initial_kappa is None else float(initial_kappa)
        resolved_lower_bound = 0.05 if lower_kappa_bound is None else float(lower_kappa_bound)
        resolved_upper_bound = 20.0 if upper_kappa_bound is None else float(upper_kappa_bound)
        if final_kappa is None:
            raise AssertionError("joint optimization lost final kappa")
        return [
            _build_parameter_row(
                parameter_name="kappa",
                initial_value=resolved_initial_kappa,
                optimized_value=final_kappa,
                lower_bound=resolved_lower_bound,
                upper_bound=resolved_upper_bound,
            )
        ]
    if final_exchangeabilities is None:
        raise AssertionError("joint optimization lost final exchangeabilities")
    resolved_initial_exchangeabilities = _normalize_initial_exchangeabilities(
        initial_exchangeabilities
    )
    resolved_lower_bound = (
        0.05 if lower_exchangeability_bound is None else float(lower_exchangeability_bound)
    )
    resolved_upper_bound = (
        20.0 if upper_exchangeability_bound is None else float(upper_exchangeability_bound)
    )
    return [
        _build_parameter_row(
            parameter_name=parameter_name,
            initial_value=resolved_initial_exchangeabilities[parameter_name],
            optimized_value=final_exchangeabilities[parameter_name],
            lower_bound=resolved_lower_bound,
            upper_bound=resolved_upper_bound,
        )
        for parameter_name in ("AG", "AT", "CG", "CT", "GT")
    ]


def _normalize_initial_exchangeabilities(
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ),
) -> dict[str, float]:
    if initial_exchangeabilities is None:
        return {
            "AG": 1.0,
            "AT": 1.0,
            "CG": 1.0,
            "CT": 1.0,
            "GT": 1.0,
        }
    if isinstance(initial_exchangeabilities, dict):
        if ("A", "C") in initial_exchangeabilities:
            named = {
                "".join(key): float(value)
                for key, value in initial_exchangeabilities.items()
            }
        else:
            named = {
                str(key): float(value)
                for key, value in initial_exchangeabilities.items()
            }
        ac_value = named.get("AC", 1.0)
        return {
            parameter_name: named[parameter_name] / ac_value
            for parameter_name in ("AG", "AT", "CG", "CT", "GT")
        }
    values = [float(value) for value in initial_exchangeabilities]
    if len(values) != 6:
        raise ValueError("initial_exchangeabilities must contain six ordered values")
    anchor = values[0]
    return {
        "AG": values[1] / anchor,
        "AT": values[2] / anchor,
        "CG": values[3] / anchor,
        "CT": values[4] / anchor,
        "GT": values[5] / anchor,
    }


def _build_parameter_row(
    *,
    parameter_name: str,
    initial_value: float,
    optimized_value: float,
    lower_bound: float,
    upper_bound: float,
) -> SubstitutionParameterOptimizationRow:
    return SubstitutionParameterOptimizationRow(
        parameter_name=parameter_name,
        initial_value=initial_value,
        optimized_value=optimized_value,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        hit_lower_bound=math.isclose(
            optimized_value,
            lower_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ),
        hit_upper_bound=math.isclose(
            optimized_value,
            upper_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ),
    )
