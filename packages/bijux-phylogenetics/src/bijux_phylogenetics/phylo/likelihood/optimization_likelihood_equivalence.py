from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.f81 import evaluate_f81_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.gtr import evaluate_gtr_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.hky85 import evaluate_hky85_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.jc69 import evaluate_jc69_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.k80 import evaluate_k80_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.models import (
    FixedTopologyNucleotideBranchLengthOptimizationReport,
    FixedTopologyNucleotideJointOptimizationReport,
    FixedTopologyNucleotideJointOptimizationRestartReport,
    FixedTopologyNucleotideSingleBranchOptimizationReport,
    NucleotideLikelihoodOptimizationEquivalenceReport,
    NucleotideSubstitutionParameterOptimizationReport,
)


def validate_nucleotide_likelihood_optimization_equivalence_tolerances(
    *,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> tuple[float, float]:
    """Validate one declared pair of optimization equivalence tolerances."""
    if not math.isfinite(absolute_tolerance) or absolute_tolerance < 0.0:
        raise ValueError(
            "nucleotide likelihood optimization equivalence absolute_tolerance must be finite and nonnegative"
        )
    if not math.isfinite(relative_tolerance) or relative_tolerance < 0.0:
        raise ValueError(
            "nucleotide likelihood optimization equivalence relative_tolerance must be finite and nonnegative"
        )
    return absolute_tolerance, relative_tolerance


def check_nucleotide_likelihood_optimization_equivalence(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
        | FixedTopologyNucleotideJointOptimizationRestartReport
    ),
    records: list[AlignmentRecord],
    *,
    absolute_tolerance: float = 1e-12,
    relative_tolerance: float = 1e-12,
) -> NucleotideLikelihoodOptimizationEquivalenceReport:
    """Independently rescore one optimized nucleotide likelihood result."""
    validated_absolute_tolerance, validated_relative_tolerance = (
        validate_nucleotide_likelihood_optimization_equivalence_tolerances(
            absolute_tolerance=absolute_tolerance,
            relative_tolerance=relative_tolerance,
        )
    )
    optimization_surface, resolved_report = _resolve_optimization_surface_and_report(
        optimization_report
    )
    parameter_values = _resolve_optimized_parameter_values(resolved_report)
    independently_rescored_log_likelihood = _independently_rescore_optimized_report(
        resolved_report,
        records,
        parameter_values=parameter_values,
    )
    stored_log_likelihood = _resolve_stored_log_likelihood(resolved_report)
    absolute_difference = abs(
        independently_rescored_log_likelihood - stored_log_likelihood
    )
    relative_difference = absolute_difference / max(
        1.0,
        abs(independently_rescored_log_likelihood),
        abs(stored_log_likelihood),
    )
    return NucleotideLikelihoodOptimizationEquivalenceReport(
        optimization_surface=optimization_surface,
        model_name=resolved_report.model_name,
        taxa=list(resolved_report.taxa),
        site_count=resolved_report.site_count,
        pattern_count=resolved_report.pattern_count,
        optimized_tree_newick=_resolve_optimized_tree_newick(resolved_report),
        parameter_values=parameter_values,
        root_prior_source=getattr(resolved_report, "root_prior_source", None),
        stored_log_likelihood=stored_log_likelihood,
        independently_rescored_log_likelihood=independently_rescored_log_likelihood,
        absolute_difference=absolute_difference,
        relative_difference=relative_difference,
        absolute_tolerance=validated_absolute_tolerance,
        relative_tolerance=validated_relative_tolerance,
        equivalent=math.isclose(
            independently_rescored_log_likelihood,
            stored_log_likelihood,
            rel_tol=validated_relative_tolerance,
            abs_tol=validated_absolute_tolerance,
        ),
    )


def check_nucleotide_likelihood_optimization_equivalence_from_alignment(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
        | FixedTopologyNucleotideJointOptimizationRestartReport
    ),
    alignment_path: Path,
    *,
    absolute_tolerance: float = 1e-12,
    relative_tolerance: float = 1e-12,
) -> NucleotideLikelihoodOptimizationEquivalenceReport:
    """Independently rescore one optimized nucleotide result from one alignment path."""
    return check_nucleotide_likelihood_optimization_equivalence(
        optimization_report,
        load_fasta_alignment(alignment_path),
        absolute_tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance,
    )


def _resolve_optimization_surface_and_report(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
        | FixedTopologyNucleotideJointOptimizationRestartReport
    ),
) -> tuple[
    str,
    (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
]:
    if isinstance(
        optimization_report,
        FixedTopologyNucleotideBranchLengthOptimizationReport,
    ):
        return (
            "fixed-topology-nucleotide-branch-length-optimization",
            optimization_report,
        )
    if isinstance(
        optimization_report,
        FixedTopologyNucleotideSingleBranchOptimizationReport,
    ):
        return (
            "fixed-topology-nucleotide-single-branch-optimization",
            optimization_report,
        )
    if isinstance(
        optimization_report,
        NucleotideSubstitutionParameterOptimizationReport,
    ):
        return (
            "nucleotide-substitution-parameter-optimization",
            optimization_report,
        )
    if isinstance(
        optimization_report,
        FixedTopologyNucleotideJointOptimizationReport,
    ):
        return (
            "fixed-topology-nucleotide-joint-optimization",
            optimization_report,
        )
    if isinstance(
        optimization_report,
        FixedTopologyNucleotideJointOptimizationRestartReport,
    ):
        return (
            "fixed-topology-nucleotide-joint-optimization-restarts",
            optimization_report.selected_report,
        )
    raise TypeError(
        "optimization_report must be a supported fixed-topology nucleotide optimization report"
    )


def _resolve_optimized_parameter_values(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
) -> dict[str, float]:
    if isinstance(
        optimization_report,
        (
            FixedTopologyNucleotideBranchLengthOptimizationReport,
            FixedTopologyNucleotideSingleBranchOptimizationReport,
        ),
    ):
        return _normalize_parameter_names(
            optimization_report.fixed_parameter_values
        )
    parameter_values = dict(optimization_report.fixed_parameter_values)
    if optimization_report.base_frequency_source is not None:
        parameter_values.update(
            {
                "A": float(optimization_report.base_frequency_a or 0.0),
                "C": float(optimization_report.base_frequency_c or 0.0),
                "G": float(optimization_report.base_frequency_g or 0.0),
                "T": float(optimization_report.base_frequency_t or 0.0),
            }
        )
    for row in optimization_report.parameter_rows:
        parameter_values[row.parameter_name] = row.optimized_value
    return parameter_values


def _normalize_parameter_names(
    parameter_values: dict[str, float],
) -> dict[str, float]:
    normalized_parameter_values: dict[str, float] = {}
    for parameter_name, value in parameter_values.items():
        if parameter_name.startswith("base_frequency_"):
            normalized_parameter_values[parameter_name.removeprefix("base_frequency_").upper()] = value
            continue
        if parameter_name.startswith("exchangeability_"):
            normalized_parameter_values[
                parameter_name.removeprefix("exchangeability_").upper()
            ] = value
            continue
        normalized_parameter_values[parameter_name] = value
    return normalized_parameter_values


def _independently_rescore_optimized_report(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
    records: list[AlignmentRecord],
    *,
    parameter_values: dict[str, float],
) -> float:
    evaluator = _resolve_model_evaluator(optimization_report.model_name)
    reevaluated_report = evaluator(
        loads_newick(_resolve_optimized_tree_newick(optimization_report)),
        records,
        **_resolve_evaluator_kwargs(
            optimization_report,
            parameter_values=parameter_values,
        ),
    )
    return float(reevaluated_report.log_likelihood)


def _resolve_model_evaluator(model_name: str):
    normalized_model_name = model_name.strip().upper()
    if normalized_model_name == "JC69":
        return evaluate_jc69_tree_likelihood
    if normalized_model_name == "K80":
        return evaluate_k80_tree_likelihood
    if normalized_model_name == "F81":
        return evaluate_f81_tree_likelihood
    if normalized_model_name == "HKY85":
        return evaluate_hky85_tree_likelihood
    if normalized_model_name == "GTR":
        return evaluate_gtr_tree_likelihood
    raise ValueError(f"unsupported optimized nucleotide model '{model_name}'")


def _resolve_evaluator_kwargs(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
    *,
    parameter_values: dict[str, float],
) -> dict[str, object]:
    model_name = optimization_report.model_name.strip().upper()
    kwargs: dict[str, object] = {}
    if model_name in {"K80", "HKY85"}:
        kwargs["kappa"] = parameter_values["kappa"]
    if model_name in {"F81", "HKY85", "GTR"}:
        kwargs["base_frequencies"] = {
            state: parameter_values[state]
            for state in ("A", "C", "G", "T")
        }
    if model_name == "GTR":
        kwargs["exchangeabilities"] = {
            label: parameter_values[label]
            for label in ("AC", "AG", "AT", "CG", "CT", "GT")
        }
    if isinstance(
        optimization_report,
        (
            FixedTopologyNucleotideBranchLengthOptimizationReport,
            FixedTopologyNucleotideSingleBranchOptimizationReport,
            FixedTopologyNucleotideJointOptimizationReport,
        ),
    ):
        kwargs["observation_policy"] = optimization_report.observation_policy
        kwargs["root_prior_policy"] = "provided"
        kwargs["root_prior"] = list(optimization_report.root_prior_values)
    return kwargs


def _resolve_optimized_tree_newick(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
) -> str:
    if isinstance(
        optimization_report,
        NucleotideSubstitutionParameterOptimizationReport,
    ):
        return optimization_report.tree_newick
    return optimization_report.optimized_tree_newick


def _resolve_stored_log_likelihood(
    optimization_report: (
        FixedTopologyNucleotideBranchLengthOptimizationReport
        | FixedTopologyNucleotideSingleBranchOptimizationReport
        | NucleotideSubstitutionParameterOptimizationReport
        | FixedTopologyNucleotideJointOptimizationReport
    ),
) -> float:
    return float(optimization_report.optimized_log_likelihood)
