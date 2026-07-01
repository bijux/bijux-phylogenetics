from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    normalize_dna_exchangeabilities_by_anchor,
)
from bijux_phylogenetics.phylo.likelihood.f81 import (
    evaluate_f81_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    optimize_gtr_exchangeabilities,
)
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    evaluate_hky85_tree_likelihood,
    optimize_hky85_kappa,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    evaluate_jc69_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    optimize_k80_kappa,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideSubstitutionParameterOptimizationReport,
    SubstitutionParameterOptimizationRow,
)
from bijux_phylogenetics.phylo.likelihood.optimization_boundary_warnings import (
    boundary_warning_messages,
    build_base_frequency_boundary_warnings,
    build_substitution_parameter_boundary_warnings,
)
from bijux_phylogenetics.phylo.likelihood.parameter_bounds import (
    validate_parameter_within_bounds,
    validate_positive_parameter_bounds,
)
from bijux_phylogenetics.phylo.likelihood.simplex_search import (
    run_simplex_likelihood_search,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_NUCLEOTIDE_SUBSTITUTION_OPTIMIZATION_MODELS = frozenset(
    {"jc69", "k80", "f81", "hky85", "gtr"}
)
_DEFAULT_KAPPA_BOUNDS = (0.05, 20.0)
_DEFAULT_EXCHANGEABILITY_BOUNDS = (0.05, 20.0)
_DEFAULT_MINIMUM_BASE_FREQUENCY = 1e-6
_BOUNDARY_TOLERANCE = 1e-9


def validate_nucleotide_substitution_optimization_model(model_name: str) -> str:
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _NUCLEOTIDE_SUBSTITUTION_OPTIMIZATION_MODELS:
        raise ValueError(
            "nucleotide substitution optimization model must be one of "
            + ", ".join(sorted(_NUCLEOTIDE_SUBSTITUTION_OPTIMIZATION_MODELS))
        )
    return normalized_model_name


def optimize_nucleotide_substitution_parameters(
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
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> NucleotideSubstitutionParameterOptimizationReport:
    """Optimize one fixed-topology nucleotide substitution surface on explicit branches."""
    normalized_model_name = validate_nucleotide_substitution_optimization_model(
        model_name
    )
    if normalized_model_name == "jc69":
        return _optimize_jc69_substitution_parameters(
            tree,
            records,
            base_frequencies=base_frequencies,
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
            initial_exchangeabilities=initial_exchangeabilities,
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
        )
    if normalized_model_name == "k80":
        return _optimize_k80_substitution_parameters(
            tree,
            records,
            base_frequencies=base_frequencies,
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
            initial_exchangeabilities=initial_exchangeabilities,
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
        )
    if normalized_model_name == "f81":
        return _optimize_f81_substitution_parameters(
            tree,
            records,
            base_frequencies=base_frequencies,
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
            initial_exchangeabilities=initial_exchangeabilities,
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
            max_coordinate_passes=max_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
    if normalized_model_name == "hky85":
        return _optimize_hky85_substitution_parameters(
            tree,
            records,
            base_frequencies=base_frequencies,
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
            initial_exchangeabilities=initial_exchangeabilities,
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
            max_coordinate_passes=max_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
    return _optimize_gtr_substitution_parameters(
        tree,
        records,
        base_frequencies=base_frequencies,
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
        initial_exchangeabilities=initial_exchangeabilities,
        lower_exchangeability_bound=lower_exchangeability_bound,
        upper_exchangeability_bound=upper_exchangeability_bound,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )


def optimize_nucleotide_substitution_parameters_from_alignment(
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
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> NucleotideSubstitutionParameterOptimizationReport:
    """Optimize nucleotide substitution parameters from one tree path and alignment."""
    return optimize_nucleotide_substitution_parameters(
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
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )


def _optimize_jc69_substitution_parameters(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
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
) -> NucleotideSubstitutionParameterOptimizationReport:
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "base_frequencies",
        base_frequencies,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "initial_kappa",
        initial_kappa,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "lower_kappa_bound",
        lower_kappa_bound,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "upper_kappa_bound",
        upper_kappa_bound,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "initial_exchangeabilities",
        initial_exchangeabilities,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "lower_exchangeability_bound",
        lower_exchangeability_bound,
    )
    _reject_irrelevant_parameter(
        "JC69 substitution optimization",
        "upper_exchangeability_bound",
        upper_exchangeability_bound,
    )

    evaluation_report = evaluate_jc69_tree_likelihood(tree, records)
    log_likelihood = evaluation_report.log_likelihood
    aic = -2.0 * log_likelihood
    return NucleotideSubstitutionParameterOptimizationReport(
        model_name="JC69",
        taxa=evaluation_report.taxa,
        site_count=evaluation_report.site_count,
        pattern_count=evaluation_report.pattern_count,
        tree_newick=evaluation_report.tree_newick,
        parameter_count=0,
        base_frequency_source=None,
        base_frequency_a=None,
        base_frequency_c=None,
        base_frequency_g=None,
        base_frequency_t=None,
        fixed_parameter_values={},
        parameter_rows=[],
        initial_log_likelihood=log_likelihood,
        optimized_log_likelihood=log_likelihood,
        initial_aic=aic,
        optimized_aic=aic,
        function_evaluation_count=1,
        optimization_pass_count=0,
        converged=True,
        boundary_warnings=[],
        warnings=[
            "JC69 has no free substitution parameters; skipping parameter search"
        ],
    )


def _optimize_k80_substitution_parameters(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
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
) -> NucleotideSubstitutionParameterOptimizationReport:
    _reject_irrelevant_parameter(
        "K80 substitution optimization",
        "base_frequencies",
        base_frequencies,
    )
    _reject_irrelevant_parameter(
        "K80 substitution optimization",
        "initial_exchangeabilities",
        initial_exchangeabilities,
    )
    _reject_irrelevant_parameter(
        "K80 substitution optimization",
        "lower_exchangeability_bound",
        lower_exchangeability_bound,
    )
    _reject_irrelevant_parameter(
        "K80 substitution optimization",
        "upper_exchangeability_bound",
        upper_exchangeability_bound,
    )

    validated_initial_kappa, validated_lower_bound, validated_upper_bound = (
        _coalesce_kappa_search_parameters(
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
        )
    )
    optimization_report = optimize_k80_kappa(
        tree,
        records,
        initial_kappa=validated_initial_kappa,
        lower_kappa_bound=validated_lower_bound,
        upper_kappa_bound=validated_upper_bound,
    )
    parameter_row, warnings = _parameter_row_with_boundary_warnings(
        parameter_name="kappa",
        initial_value=validated_initial_kappa,
        optimized_value=optimization_report.optimized_kappa,
        lower_bound=optimization_report.lower_kappa_bound,
        upper_bound=optimization_report.upper_kappa_bound,
    )
    if not optimization_report.converged:
        warnings.append("K80 substitution-parameter search did not converge")
    boundary_warnings = build_substitution_parameter_boundary_warnings([parameter_row])
    return NucleotideSubstitutionParameterOptimizationReport(
        model_name="K80",
        taxa=optimization_report.taxa,
        site_count=optimization_report.site_count,
        pattern_count=optimization_report.pattern_count,
        tree_newick=optimization_report.tree_newick,
        parameter_count=1,
        base_frequency_source=None,
        base_frequency_a=None,
        base_frequency_c=None,
        base_frequency_g=None,
        base_frequency_t=None,
        fixed_parameter_values={},
        parameter_rows=[parameter_row],
        initial_log_likelihood=optimization_report.initial_log_likelihood,
        optimized_log_likelihood=optimization_report.optimized_log_likelihood,
        initial_aic=optimization_report.initial_aic,
        optimized_aic=optimization_report.optimized_aic,
        function_evaluation_count=optimization_report.function_evaluation_count,
        optimization_pass_count=1,
        converged=optimization_report.converged,
        boundary_warnings=boundary_warnings,
        warnings=warnings,
    )


def _optimize_f81_substitution_parameters(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
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
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> NucleotideSubstitutionParameterOptimizationReport:
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "initial_kappa",
        initial_kappa,
    )
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "lower_kappa_bound",
        lower_kappa_bound,
    )
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "upper_kappa_bound",
        upper_kappa_bound,
    )
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "initial_exchangeabilities",
        initial_exchangeabilities,
    )
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "lower_exchangeability_bound",
        lower_exchangeability_bound,
    )
    _reject_irrelevant_parameter(
        "F81 substitution optimization",
        "upper_exchangeability_bound",
        upper_exchangeability_bound,
    )

    if base_frequencies is not None:
        evaluation_report = evaluate_f81_tree_likelihood(
            tree,
            records,
            base_frequencies=base_frequencies,
        )
        initial_log_likelihood = evaluation_report.log_likelihood
        optimized_log_likelihood = evaluation_report.log_likelihood
        initial_aic = evaluation_report.aic
        optimized_aic = evaluation_report.aic
        optimization_pass_count = 0
        function_evaluation_count = 1
        converged = True
        warnings: list[str] = []
    else:
        initial_report = evaluate_f81_tree_likelihood(tree, records)
        search_result = run_simplex_likelihood_search(
            initial_values_by_name={
                "A": initial_report.base_frequency_a,
                "C": initial_report.base_frequency_c,
                "G": initial_report.base_frequency_g,
                "T": initial_report.base_frequency_t,
            },
            expected_component_names=DNA_STATE_ORDER,
            owner_name="F81 substitution optimization",
            evaluate=lambda candidate: _evaluate_f81_frequency_candidate(
                tree,
                records,
                candidate,
            ),
            minimum_component_value=_DEFAULT_MINIMUM_BASE_FREQUENCY,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
        evaluation_report = replace(
            search_result.payload,
            base_frequency_source="estimated",
        )
        initial_log_likelihood = initial_report.log_likelihood
        optimized_log_likelihood = evaluation_report.log_likelihood
        initial_aic = initial_report.aic
        optimized_aic = evaluation_report.aic
        optimization_pass_count = search_result.optimization_pass_count
        function_evaluation_count = search_result.function_evaluation_count
        converged = search_result.converged
        warnings = []
        if not search_result.converged:
            warnings.append("F81 frequency search did not converge")
    boundary_warnings = build_base_frequency_boundary_warnings(
        base_frequency_source=evaluation_report.base_frequency_source,
        base_frequency_a=evaluation_report.base_frequency_a,
        base_frequency_c=evaluation_report.base_frequency_c,
        base_frequency_g=evaluation_report.base_frequency_g,
        base_frequency_t=evaluation_report.base_frequency_t,
        lower_frequency_bound=(
            _DEFAULT_MINIMUM_BASE_FREQUENCY
            if base_frequencies is None
            else 0.0
        ),
    )
    warnings.extend(
        warning
        for warning in boundary_warning_messages(boundary_warnings)
        if warning not in warnings
    )
    return NucleotideSubstitutionParameterOptimizationReport(
        model_name="F81",
        taxa=evaluation_report.taxa,
        site_count=evaluation_report.site_count,
        pattern_count=evaluation_report.pattern_count,
        tree_newick=evaluation_report.tree_newick,
        parameter_count=evaluation_report.parameter_count,
        base_frequency_source=evaluation_report.base_frequency_source,
        base_frequency_a=evaluation_report.base_frequency_a,
        base_frequency_c=evaluation_report.base_frequency_c,
        base_frequency_g=evaluation_report.base_frequency_g,
        base_frequency_t=evaluation_report.base_frequency_t,
        fixed_parameter_values={},
        parameter_rows=[],
        initial_log_likelihood=initial_log_likelihood,
        optimized_log_likelihood=optimized_log_likelihood,
        initial_aic=initial_aic,
        optimized_aic=optimized_aic,
        function_evaluation_count=function_evaluation_count,
        optimization_pass_count=optimization_pass_count,
        converged=converged,
        boundary_warnings=boundary_warnings,
        warnings=warnings,
    )


def _optimize_hky85_substitution_parameters(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
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
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> NucleotideSubstitutionParameterOptimizationReport:
    _reject_irrelevant_parameter(
        "HKY85 substitution optimization",
        "initial_exchangeabilities",
        initial_exchangeabilities,
    )
    _reject_irrelevant_parameter(
        "HKY85 substitution optimization",
        "lower_exchangeability_bound",
        lower_exchangeability_bound,
    )
    _reject_irrelevant_parameter(
        "HKY85 substitution optimization",
        "upper_exchangeability_bound",
        upper_exchangeability_bound,
    )

    validated_initial_kappa, validated_lower_bound, validated_upper_bound = (
        _coalesce_kappa_search_parameters(
            initial_kappa=initial_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
        )
    )
    if base_frequencies is not None:
        optimization_report = optimize_hky85_kappa(
            tree,
            records,
            base_frequencies=base_frequencies,
            initial_kappa=validated_initial_kappa,
            lower_kappa_bound=validated_lower_bound,
            upper_kappa_bound=validated_upper_bound,
        )
        optimization_pass_count = 1
        function_evaluation_count = optimization_report.function_evaluation_count
        converged = optimization_report.converged
    else:
        (
            optimization_report,
            function_evaluation_count,
            optimization_pass_count,
            converged,
        ) = _optimize_estimated_hky85_kappa_and_base_frequencies(
            tree,
            records,
            initial_kappa=validated_initial_kappa,
            lower_kappa_bound=validated_lower_bound,
            upper_kappa_bound=validated_upper_bound,
            max_joint_passes=max_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
    parameter_row, warnings = _parameter_row_with_boundary_warnings(
        parameter_name="kappa",
        initial_value=validated_initial_kappa,
        optimized_value=optimization_report.optimized_kappa,
        lower_bound=optimization_report.lower_kappa_bound,
        upper_bound=optimization_report.upper_kappa_bound,
    )
    if not converged:
        warnings.append("HKY85 substitution-parameter search did not converge")
    boundary_warnings = build_substitution_parameter_boundary_warnings([parameter_row])
    boundary_warnings.extend(
        build_base_frequency_boundary_warnings(
            base_frequency_source=optimization_report.base_frequency_source,
            base_frequency_a=optimization_report.base_frequency_a,
            base_frequency_c=optimization_report.base_frequency_c,
            base_frequency_g=optimization_report.base_frequency_g,
            base_frequency_t=optimization_report.base_frequency_t,
            lower_frequency_bound=(
                _DEFAULT_MINIMUM_BASE_FREQUENCY
                if base_frequencies is None
                else 0.0
            ),
        )
    )
    warnings.extend(
        warning
        for warning in boundary_warning_messages(boundary_warnings)
        if warning not in warnings
    )
    return NucleotideSubstitutionParameterOptimizationReport(
        model_name="HKY85",
        taxa=optimization_report.taxa,
        site_count=optimization_report.site_count,
        pattern_count=optimization_report.pattern_count,
        tree_newick=optimization_report.tree_newick,
        parameter_count=optimization_report.parameter_count,
        base_frequency_source=optimization_report.base_frequency_source,
        base_frequency_a=optimization_report.base_frequency_a,
        base_frequency_c=optimization_report.base_frequency_c,
        base_frequency_g=optimization_report.base_frequency_g,
        base_frequency_t=optimization_report.base_frequency_t,
        fixed_parameter_values={},
        parameter_rows=[parameter_row],
        initial_log_likelihood=optimization_report.initial_log_likelihood,
        optimized_log_likelihood=optimization_report.optimized_log_likelihood,
        initial_aic=optimization_report.initial_aic,
        optimized_aic=optimization_report.optimized_aic,
        function_evaluation_count=function_evaluation_count,
        optimization_pass_count=optimization_pass_count,
        converged=converged,
        boundary_warnings=boundary_warnings,
        warnings=warnings,
    )


def _optimize_gtr_substitution_parameters(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
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
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> NucleotideSubstitutionParameterOptimizationReport:
    _reject_irrelevant_parameter(
        "GTR substitution optimization",
        "initial_kappa",
        initial_kappa,
    )
    _reject_irrelevant_parameter(
        "GTR substitution optimization",
        "lower_kappa_bound",
        lower_kappa_bound,
    )
    _reject_irrelevant_parameter(
        "GTR substitution optimization",
        "upper_kappa_bound",
        upper_kappa_bound,
    )

    validated_lower_bound, validated_upper_bound = (
        _coalesce_exchangeability_search_bounds(
            lower_exchangeability_bound=lower_exchangeability_bound,
            upper_exchangeability_bound=upper_exchangeability_bound,
        )
    )
    if initial_exchangeabilities is None:
        normalized_initial_exchangeabilities = numpy.ones(6, dtype=float)
    else:
        normalized_initial_exchangeabilities = (
            normalize_dna_exchangeabilities_by_anchor(
                initial_exchangeabilities,
                model_name="GTR substitution optimization",
            )
        )
    optimization_report = optimize_gtr_exchangeabilities(
        tree,
        records,
        base_frequencies=base_frequencies,
        initial_exchangeabilities=initial_exchangeabilities,
        lower_exchangeability_bound=validated_lower_bound,
        upper_exchangeability_bound=validated_upper_bound,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )
    parameter_rows: list[SubstitutionParameterOptimizationRow] = []
    warnings: list[str] = []
    optimized_values = {
        "AG": optimization_report.exchangeability_ag,
        "AT": optimization_report.exchangeability_at,
        "CG": optimization_report.exchangeability_cg,
        "CT": optimization_report.exchangeability_ct,
        "GT": optimization_report.exchangeability_gt,
    }
    for parameter_index, parameter_name in enumerate(
        ("AG", "AT", "CG", "CT", "GT"), start=1
    ):
        parameter_row, parameter_warnings = _parameter_row_with_boundary_warnings(
            parameter_name=parameter_name,
            initial_value=float(normalized_initial_exchangeabilities[parameter_index]),
            optimized_value=optimized_values[parameter_name],
            lower_bound=optimization_report.lower_exchangeability_bound,
            upper_bound=optimization_report.upper_exchangeability_bound,
        )
        parameter_rows.append(parameter_row)
        warnings.extend(parameter_warnings)
    if not optimization_report.converged:
        warnings.append("GTR substitution-parameter search did not converge")
    boundary_warnings = build_substitution_parameter_boundary_warnings(parameter_rows)
    boundary_warnings.extend(
        build_base_frequency_boundary_warnings(
            base_frequency_source=optimization_report.base_frequency_source,
            base_frequency_a=optimization_report.base_frequency_a,
            base_frequency_c=optimization_report.base_frequency_c,
            base_frequency_g=optimization_report.base_frequency_g,
            base_frequency_t=optimization_report.base_frequency_t,
        )
    )
    warnings.extend(
        warning
        for warning in boundary_warning_messages(boundary_warnings)
        if warning not in warnings
    )
    return NucleotideSubstitutionParameterOptimizationReport(
        model_name="GTR",
        taxa=optimization_report.taxa,
        site_count=optimization_report.site_count,
        pattern_count=optimization_report.pattern_count,
        tree_newick=optimization_report.tree_newick,
        parameter_count=optimization_report.parameter_count,
        base_frequency_source=optimization_report.base_frequency_source,
        base_frequency_a=optimization_report.base_frequency_a,
        base_frequency_c=optimization_report.base_frequency_c,
        base_frequency_g=optimization_report.base_frequency_g,
        base_frequency_t=optimization_report.base_frequency_t,
        fixed_parameter_values={"AC": optimization_report.exchangeability_ac},
        parameter_rows=parameter_rows,
        initial_log_likelihood=optimization_report.initial_log_likelihood,
        optimized_log_likelihood=optimization_report.optimized_log_likelihood,
        initial_aic=optimization_report.initial_aic,
        optimized_aic=optimization_report.optimized_aic,
        function_evaluation_count=optimization_report.function_evaluation_count,
        optimization_pass_count=optimization_report.optimization_pass_count,
        converged=optimization_report.converged,
        boundary_warnings=boundary_warnings,
        warnings=warnings,
    )


def _evaluate_f81_frequency_candidate(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    candidate_base_frequencies: dict[str, float],
) -> tuple[object, float]:
    report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=candidate_base_frequencies,
    )
    return report, report.log_likelihood


def _optimize_estimated_hky85_kappa_and_base_frequencies(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    initial_kappa: float,
    lower_kappa_bound: float,
    upper_kappa_bound: float,
    max_joint_passes: int,
    improvement_tolerance: float,
):
    current_optimization_report = optimize_hky85_kappa(
        tree,
        records,
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
    )
    current_base_frequencies = _resolved_report_base_frequencies(
        current_optimization_report
    )
    current_kappa = current_optimization_report.optimized_kappa
    current_objective = current_optimization_report.optimized_log_likelihood
    total_function_evaluation_count = (
        current_optimization_report.function_evaluation_count
    )
    scalar_searches_converged = current_optimization_report.converged
    executed_joint_pass_count = 0
    schedule_converged = False

    for joint_pass_index in range(1, max_joint_passes + 1):
        executed_joint_pass_count = joint_pass_index
        frequency_search = run_simplex_likelihood_search(
            initial_values_by_name=current_base_frequencies,
            expected_component_names=DNA_STATE_ORDER,
            owner_name="HKY85 substitution optimization",
            evaluate=lambda candidate: _evaluate_hky85_frequency_candidate(
                tree,
                records,
                candidate,
                current_kappa=current_kappa,
            ),
            minimum_component_value=_DEFAULT_MINIMUM_BASE_FREQUENCY,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=1,
        )
        total_function_evaluation_count += frequency_search.function_evaluation_count
        current_base_frequencies = _resolved_report_base_frequencies(
            frequency_search.payload
        )

        kappa_report = optimize_hky85_kappa(
            tree,
            records,
            base_frequencies=current_base_frequencies,
            initial_kappa=current_kappa,
            lower_kappa_bound=lower_kappa_bound,
            upper_kappa_bound=upper_kappa_bound,
        )
        total_function_evaluation_count += kappa_report.function_evaluation_count
        scalar_searches_converged = scalar_searches_converged and kappa_report.converged
        if (
            frequency_search.payload.log_likelihood
            <= (current_objective + improvement_tolerance)
            and kappa_report.optimized_log_likelihood
            <= (current_objective + improvement_tolerance)
        ):
            schedule_converged = True
            break
        current_kappa = kappa_report.optimized_kappa
        current_objective = kappa_report.optimized_log_likelihood
        current_optimization_report = replace(
            kappa_report,
            base_frequency_source="estimated",
            base_frequency_a=current_base_frequencies["A"],
            base_frequency_c=current_base_frequencies["C"],
            base_frequency_g=current_base_frequencies["G"],
            base_frequency_t=current_base_frequencies["T"],
        )

    return (
        current_optimization_report,
        total_function_evaluation_count,
        max(1, executed_joint_pass_count),
        schedule_converged and scalar_searches_converged,
    )


def _evaluate_hky85_frequency_candidate(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    candidate_base_frequencies: dict[str, float],
    *,
    current_kappa: float,
) -> tuple[object, float]:
    report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        base_frequencies=candidate_base_frequencies,
        kappa=current_kappa,
    )
    return report, report.log_likelihood


def _resolved_report_base_frequencies(report) -> dict[str, float]:
    return {
        "A": float(report.base_frequency_a),
        "C": float(report.base_frequency_c),
        "G": float(report.base_frequency_g),
        "T": float(report.base_frequency_t),
    }


def _parameter_row_with_boundary_warnings(
    *,
    parameter_name: str,
    initial_value: float,
    optimized_value: float,
    lower_bound: float,
    upper_bound: float,
) -> tuple[SubstitutionParameterOptimizationRow, list[str]]:
    hit_lower_bound = math.isclose(
        optimized_value,
        lower_bound,
        rel_tol=0.0,
        abs_tol=_BOUNDARY_TOLERANCE,
    )
    hit_upper_bound = math.isclose(
        optimized_value,
        upper_bound,
        rel_tol=0.0,
        abs_tol=_BOUNDARY_TOLERANCE,
    )
    warnings: list[str] = []
    if hit_lower_bound:
        warnings.append(f"{parameter_name} hit lower search boundary")
    if hit_upper_bound:
        warnings.append(f"{parameter_name} hit upper search boundary")
    return (
        SubstitutionParameterOptimizationRow(
            parameter_name=parameter_name,
            initial_value=initial_value,
            optimized_value=optimized_value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            hit_lower_bound=hit_lower_bound,
            hit_upper_bound=hit_upper_bound,
        ),
        warnings,
    )


def _coalesce_kappa_search_parameters(
    *,
    initial_kappa: float | None,
    lower_kappa_bound: float | None,
    upper_kappa_bound: float | None,
) -> tuple[float, float, float]:
    resolved_lower_bound, resolved_upper_bound = validate_positive_parameter_bounds(
        parameter_name="kappa",
        lower_bound=(
            _DEFAULT_KAPPA_BOUNDS[0]
            if lower_kappa_bound is None
            else float(lower_kappa_bound)
        ),
        upper_bound=(
            _DEFAULT_KAPPA_BOUNDS[1]
            if upper_kappa_bound is None
            else float(upper_kappa_bound)
        ),
        owner_name="nucleotide substitution optimization",
    )
    resolved_initial_kappa = validate_parameter_within_bounds(
        parameter_name="kappa",
        value=(1.0 if initial_kappa is None else float(initial_kappa)),
        lower_bound=resolved_lower_bound,
        upper_bound=resolved_upper_bound,
        owner_name="nucleotide substitution optimization",
    )
    return (
        resolved_initial_kappa,
        resolved_lower_bound,
        resolved_upper_bound,
    )


def _coalesce_exchangeability_search_bounds(
    *,
    lower_exchangeability_bound: float | None,
    upper_exchangeability_bound: float | None,
) -> tuple[float, float]:
    return validate_positive_parameter_bounds(
        parameter_name="exchangeability",
        lower_bound=(
            _DEFAULT_EXCHANGEABILITY_BOUNDS[0]
            if lower_exchangeability_bound is None
            else float(lower_exchangeability_bound)
        ),
        upper_bound=(
            _DEFAULT_EXCHANGEABILITY_BOUNDS[1]
            if upper_exchangeability_bound is None
            else float(upper_exchangeability_bound)
        ),
        owner_name="GTR substitution optimization",
    )


def _reject_irrelevant_parameter(
    owner_name: str,
    parameter_name: str,
    value: object,
) -> None:
    if value is not None:
        raise ValueError(
            f"{owner_name} does not accept '{parameter_name}' because that model does not use it"
        )
