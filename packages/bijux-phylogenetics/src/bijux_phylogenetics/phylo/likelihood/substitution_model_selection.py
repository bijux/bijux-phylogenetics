from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_INDEX,
    UNIFORM_DNA_ROOT_PRIOR,
    estimate_empirical_dna_base_frequencies,
    evaluate_fixed_topology_dna_site_log_likelihood,
    normalize_dna_exchangeabilities_by_anchor,
    normalize_unambiguous_dna_records,
    validate_positive_kappa,
)
from bijux_phylogenetics.phylo.likelihood.f81 import (
    f81_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.gamma import (
    build_discrete_gamma_rate_categories,
)
from bijux_phylogenetics.phylo.likelihood.gtr import (
    gtr_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.hky85 import (
    hky85_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.k80 import (
    k80_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    SubstitutionModelSelectionReport,
    SubstitutionModelSelectionRow,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.substitution_parameters import (
    optimize_nucleotide_substitution_parameters,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
    InvalidBranchLengthError,
)

_BASE_MODEL_ORDER = ("JC69", "K80", "F81", "HKY85", "GTR")
_RATE_HETEROGENEITY_ORDER = (
    ("fixed-rate", ""),
    ("discrete-gamma", "+G"),
    ("invariant", "+I"),
    ("discrete-gamma-invariant", "+G+I"),
)
_DEFAULT_ALPHA = 0.8
_DEFAULT_INVARIANT_PROPORTION = 0.1
_DEFAULT_ALPHA_BOUNDS = (0.05, 20.0)
_DEFAULT_INVARIANT_BOUNDS = (0.0, 0.95)
_DEFAULT_KAPPA_BOUNDS = (0.05, 20.0)
_DEFAULT_EXCHANGEABILITY_BOUNDS = (0.05, 20.0)
_BOUNDARY_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class _CandidateSpecification:
    model_name: str
    base_model_name: str
    rate_heterogeneity_model: str


def default_substitution_model_selection_candidates() -> tuple[str, ...]:
    """Return the default fixed-topology nucleotide model-selection candidate set."""
    return tuple(
        f"{base_model_name}{suffix}"
        for base_model_name in _BASE_MODEL_ORDER
        for _rate_model, suffix in _RATE_HETEROGENEITY_ORDER
    )


def compare_nucleotide_substitution_models(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    candidate_models: tuple[str, ...] | None = None,
    category_count: int = 4,
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> SubstitutionModelSelectionReport:
    """Compare fixed-topology nucleotide substitution candidates by information criteria."""
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name="substitution model selection",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    validate_explicit_branch_lengths(tree, model_name="substitution model selection")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="substitution model selection",
    )
    empirical_base_frequencies = estimate_empirical_dna_base_frequencies(
        normalized_records
    )
    selected_candidates = (
        default_substitution_model_selection_candidates()
        if candidate_models is None
        else candidate_models
    )
    rows = [
        _fit_substitution_model_candidate(
            tree,
            normalized_records,
            compressed_patterns,
            empirical_base_frequencies=empirical_base_frequencies,
            candidate_model_name=candidate_model_name,
            category_count=category_count,
            max_coordinate_passes=max_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
        for candidate_model_name in selected_candidates
    ]
    warnings = _rank_substitution_model_selection_rows(
        rows, sample_size=compressed_patterns.alignment_length
    )
    return SubstitutionModelSelectionReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        tree_newick=dumps_newick(tree),
        rows=rows,
        best_model_aic=next(
            (row.model_name for row in rows if row.selected_by_aic), None
        ),
        best_model_aicc=next(
            (row.model_name for row in rows if row.selected_by_aicc), None
        ),
        best_model_bic=next(
            (row.model_name for row in rows if row.selected_by_bic), None
        ),
        warnings=warnings,
    )


def compare_nucleotide_substitution_models_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    candidate_models: tuple[str, ...] | None = None,
    category_count: int = 4,
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> SubstitutionModelSelectionReport:
    """Compare fixed-topology nucleotide substitution candidates from file paths."""
    return compare_nucleotide_substitution_models(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        candidate_models=candidate_models,
        category_count=category_count,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )


def _fit_substitution_model_candidate(
    tree: PhyloTree,
    normalized_records: list[AlignmentRecord],
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    empirical_base_frequencies: numpy.ndarray,
    candidate_model_name: str,
    category_count: int,
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> SubstitutionModelSelectionRow:
    try:
        specification = _parse_candidate_model_name(candidate_model_name)
        if specification.rate_heterogeneity_model == "fixed-rate":
            return _fit_fixed_rate_candidate(
                tree,
                normalized_records,
                compressed_patterns,
                _empirical_base_frequencies=empirical_base_frequencies,
                specification=specification,
            )
        return _fit_rate_heterogeneous_candidate(
            tree,
            compressed_patterns,
            empirical_base_frequencies=empirical_base_frequencies,
            specification=specification,
            category_count=category_count,
            max_coordinate_passes=max_coordinate_passes,
            improvement_tolerance=improvement_tolerance,
        )
    except (
        ValueError,
        InvalidAlignmentError,
        AlignmentTaxonMismatchError,
        InvalidBranchLengthError,
    ) as error:
        specification = _coerce_failed_candidate(candidate_model_name)
        return SubstitutionModelSelectionRow(
            model_name=specification.model_name,
            base_model_name=specification.base_model_name,
            rate_heterogeneity_model=specification.rate_heterogeneity_model,
            fit_succeeded=False,
            parameter_count=None,
            log_likelihood=None,
            aic=None,
            aicc=None,
            bic=None,
            delta_aic=None,
            akaike_weight=None,
            rank=None,
            comparable_on_aic=False,
            comparable_on_aicc=False,
            comparable_on_bic=False,
            selected_by_aic=False,
            selected_by_aicc=False,
            selected_by_bic=False,
            parameter_values={},
            warnings=[str(error)],
        )


def _fit_fixed_rate_candidate(
    tree: PhyloTree,
    normalized_records: list[AlignmentRecord],
    _compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    _empirical_base_frequencies: numpy.ndarray,
    specification: _CandidateSpecification,
) -> SubstitutionModelSelectionRow:
    fit_report = optimize_nucleotide_substitution_parameters(
        tree,
        normalized_records,
        model_name=specification.base_model_name.lower(),
    )
    return _successful_selection_row(
        specification=specification,
        parameter_count=fit_report.parameter_count,
        log_likelihood=fit_report.optimized_log_likelihood,
        parameter_values={
            **(
                {}
                if fit_report.base_frequency_a is None
                else {
                    "base_frequency_a": fit_report.base_frequency_a,
                    "base_frequency_c": fit_report.base_frequency_c
                    if fit_report.base_frequency_c is not None
                    else 0.0,
                    "base_frequency_g": fit_report.base_frequency_g
                    if fit_report.base_frequency_g is not None
                    else 0.0,
                    "base_frequency_t": fit_report.base_frequency_t
                    if fit_report.base_frequency_t is not None
                    else 0.0,
                }
            ),
            **fit_report.fixed_parameter_values,
            **{
                row.parameter_name: row.optimized_value
                for row in fit_report.parameter_rows
            },
        },
        warnings=list(fit_report.warnings),
    )


def _fit_rate_heterogeneous_candidate(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    empirical_base_frequencies: numpy.ndarray,
    specification: _CandidateSpecification,
    category_count: int,
    max_coordinate_passes: int,
    improvement_tolerance: float,
) -> SubstitutionModelSelectionRow:
    initial_values: dict[str, float] = {}
    bounds_by_name: dict[str, tuple[float, float]] = {}
    parameter_count = _base_model_parameter_count(specification.base_model_name)
    parameter_values: dict[str, float] = {}

    if specification.base_model_name in {"F81", "HKY85", "GTR"}:
        parameter_values.update(
            _base_frequency_parameter_values(empirical_base_frequencies)
        )
    if specification.base_model_name in {"K80", "HKY85"}:
        initial_values["kappa"] = 1.0
        bounds_by_name["kappa"] = _DEFAULT_KAPPA_BOUNDS
    if specification.base_model_name == "GTR":
        for parameter_name in ("AG", "AT", "CG", "CT", "GT"):
            initial_values[parameter_name] = 1.0
            bounds_by_name[parameter_name] = _DEFAULT_EXCHANGEABILITY_BOUNDS
    if specification.rate_heterogeneity_model in {
        "discrete-gamma",
        "discrete-gamma-invariant",
    }:
        initial_values["alpha"] = _DEFAULT_ALPHA
        bounds_by_name["alpha"] = _DEFAULT_ALPHA_BOUNDS
        parameter_count += 1
    if specification.rate_heterogeneity_model in {
        "invariant",
        "discrete-gamma-invariant",
    }:
        initial_values["invariant_proportion"] = _DEFAULT_INVARIANT_PROPORTION
        bounds_by_name["invariant_proportion"] = _DEFAULT_INVARIANT_BOUNDS
        parameter_count += 1

    def evaluate_candidate(
        candidate_values: dict[str, float],
    ) -> tuple[dict[str, float], float]:
        resolved_parameters = dict(parameter_values)
        resolved_parameters.update(candidate_values)
        log_likelihood = _evaluate_candidate_log_likelihood(
            tree,
            compressed_patterns,
            empirical_base_frequencies=empirical_base_frequencies,
            specification=specification,
            parameter_values=resolved_parameters,
            category_count=category_count,
        )
        return resolved_parameters, log_likelihood

    if initial_values:
        search_result = run_bounded_coordinate_likelihood_search(
            initial_values=initial_values,
            bounds_by_name=bounds_by_name,
            evaluate=evaluate_candidate,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
        resolved_parameters = dict(search_result.payload)
        log_likelihood = float(search_result.objective_value)
        converged = search_result.converged
        function_evaluation_count = search_result.function_evaluation_count
        optimization_pass_count = search_result.optimization_pass_count
    else:
        resolved_parameters, raw_log_likelihood = evaluate_candidate({})
        log_likelihood = float(raw_log_likelihood)
        converged = True
        function_evaluation_count = 1
        optimization_pass_count = 0
    warnings = _selection_parameter_warnings(
        specification=specification,
        parameter_values=resolved_parameters,
        bounds_by_name=bounds_by_name,
        converged=converged,
        optimization_pass_count=optimization_pass_count,
        function_evaluation_count=function_evaluation_count,
    )
    return _successful_selection_row(
        specification=specification,
        parameter_count=parameter_count,
        log_likelihood=log_likelihood,
        parameter_values=resolved_parameters,
        warnings=warnings,
    )


def _parse_candidate_model_name(candidate_model_name: str) -> _CandidateSpecification:
    normalized_name = candidate_model_name.strip().upper()
    if not normalized_name:
        raise ValueError("candidate substitution model name must not be empty")
    for _rate_heterogeneity_model, suffix in reversed(_RATE_HETEROGENEITY_ORDER):
        if suffix and normalized_name.endswith(suffix):
            rate_heterogeneity_model = _rate_heterogeneity_model
            base_model_name = normalized_name[: -len(suffix)]
            break
    else:
        rate_heterogeneity_model = "fixed-rate"
        base_model_name = normalized_name
    if base_model_name == "HKY":
        base_model_name = "HKY85"
    if base_model_name not in _BASE_MODEL_ORDER:
        raise ValueError(
            "candidate substitution model must use one of "
            + ", ".join(_BASE_MODEL_ORDER)
        )
    canonical_model_name = base_model_name + next(
        suffix
        for model, suffix in _RATE_HETEROGENEITY_ORDER
        if model == rate_heterogeneity_model
    )
    return _CandidateSpecification(
        model_name=canonical_model_name,
        base_model_name=base_model_name,
        rate_heterogeneity_model=rate_heterogeneity_model,
    )


def _coerce_failed_candidate(candidate_model_name: str) -> _CandidateSpecification:
    try:
        return _parse_candidate_model_name(candidate_model_name)
    except ValueError:
        return _CandidateSpecification(
            model_name=candidate_model_name.strip() or "<empty>",
            base_model_name="unknown",
            rate_heterogeneity_model="unknown",
        )


def _evaluate_candidate_log_likelihood(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    empirical_base_frequencies: numpy.ndarray,
    specification: _CandidateSpecification,
    parameter_values: dict[str, float],
    category_count: int,
) -> float:
    root_prior, transition_matrix_for_scaled_branch_length = (
        _resolve_candidate_transition_surface(
            empirical_base_frequencies=empirical_base_frequencies,
            specification=specification,
            parameter_values=parameter_values,
        )
    )
    category_rates = None
    if specification.rate_heterogeneity_model in {
        "discrete-gamma",
        "discrete-gamma-invariant",
    }:
        category_rates = build_discrete_gamma_rate_categories(
            alpha=float(parameter_values["alpha"]),
            category_count=category_count,
        )
    invariant_proportion = (
        None
        if specification.rate_heterogeneity_model == "fixed-rate"
        or specification.rate_heterogeneity_model == "discrete-gamma"
        else float(parameter_values["invariant_proportion"])
    )
    total_log_likelihood = 0.0
    for pattern in compressed_patterns.patterns:
        site_likelihood = _candidate_pattern_likelihood(
            tree,
            pattern.states,
            taxon_order=compressed_patterns.taxon_order,
            specification=specification,
            root_prior=root_prior,
            transition_matrix_for_scaled_branch_length=transition_matrix_for_scaled_branch_length,
            category_rates=category_rates,
            invariant_proportion=invariant_proportion,
        )
        if site_likelihood <= 0.0 or not math.isfinite(site_likelihood):
            raise ValueError(
                f"{specification.model_name} retained a nonpositive site likelihood"
            )
        total_log_likelihood += pattern.weight * math.log(site_likelihood)
    return total_log_likelihood


def _candidate_pattern_likelihood(
    tree: PhyloTree,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
    specification: _CandidateSpecification,
    root_prior: numpy.ndarray,
    transition_matrix_for_scaled_branch_length,
    category_rates,
    invariant_proportion: float | None,
) -> float:
    if category_rates is None:
        variable_component_likelihood = math.exp(
            evaluate_fixed_topology_dna_site_log_likelihood(
                tree,
                states,
                taxon_order=taxon_order,
                model_name=specification.model_name,
                observation_policy="reject",
                root_prior=root_prior,
                transition_matrix_for_child=lambda child: (
                    transition_matrix_for_scaled_branch_length(
                        float(child.branch_length or 0.0),
                        1.0,
                    )
                ),
            )
        )
    else:
        variable_component_likelihood = sum(
            category.weight
            * math.exp(
                evaluate_fixed_topology_dna_site_log_likelihood(
                    tree,
                    states,
                    taxon_order=taxon_order,
                    model_name=specification.model_name,
                    observation_policy="reject",
                    root_prior=root_prior,
                    transition_matrix_for_child=lambda child, rate=category.rate: (
                        transition_matrix_for_scaled_branch_length(
                            float(child.branch_length or 0.0),
                            rate,
                        )
                    ),
                )
            )
            for category in category_rates
        )
    if invariant_proportion is None:
        return variable_component_likelihood
    invariant_component_likelihood = _dna_invariant_component_likelihood(
        states,
        root_prior=root_prior,
    )
    return invariant_proportion * invariant_component_likelihood + (
        (1.0 - invariant_proportion) * variable_component_likelihood
    )


def _dna_invariant_component_likelihood(
    states: tuple[str, ...],
    *,
    root_prior: numpy.ndarray,
) -> float:
    observed_states = set(states)
    if len(observed_states) != 1:
        return 0.0
    state = next(iter(observed_states))
    return float(root_prior[DNA_STATE_INDEX[state]])


def _resolve_candidate_transition_surface(
    *,
    empirical_base_frequencies: numpy.ndarray,
    specification: _CandidateSpecification,
    parameter_values: dict[str, float],
):
    if specification.base_model_name == "JC69":
        return (
            UNIFORM_DNA_ROOT_PRIOR,
            lambda branch_length, rate: jc69_transition_probability_matrix(
                branch_length * rate
            ),
        )
    if specification.base_model_name == "K80":
        kappa = validate_positive_kappa(
            float(parameter_values["kappa"]),
            model_name=specification.model_name,
        )
        return (
            UNIFORM_DNA_ROOT_PRIOR,
            lambda branch_length, rate: k80_transition_probability_matrix(
                branch_length * rate, kappa=kappa
            ),
        )
    if specification.base_model_name == "F81":
        stationary = empirical_base_frequencies
        return (
            stationary,
            lambda branch_length, rate: f81_transition_probability_matrix(
                branch_length * rate,
                base_frequencies=stationary,
            ),
        )
    if specification.base_model_name == "HKY85":
        stationary = empirical_base_frequencies
        kappa = validate_positive_kappa(
            float(parameter_values["kappa"]),
            model_name=specification.model_name,
        )
        return (
            stationary,
            lambda branch_length, rate: hky85_transition_probability_matrix(
                branch_length * rate,
                base_frequencies=stationary,
                kappa=kappa,
            ),
        )
    stationary = empirical_base_frequencies
    exchangeabilities = normalize_dna_exchangeabilities_by_anchor(
        {
            "AC": 1.0,
            "AG": float(parameter_values["AG"]),
            "AT": float(parameter_values["AT"]),
            "CG": float(parameter_values["CG"]),
            "CT": float(parameter_values["CT"]),
            "GT": float(parameter_values["GT"]),
        },
        model_name=specification.model_name,
    )
    return stationary, lambda branch_length, rate: gtr_transition_probability_matrix(
        branch_length * rate,
        exchangeabilities=exchangeabilities,
        base_frequencies=stationary,
    )


def _selection_parameter_warnings(
    *,
    specification: _CandidateSpecification,
    parameter_values: dict[str, float],
    bounds_by_name: dict[str, tuple[float, float]],
    converged: bool,
    optimization_pass_count: int,
    function_evaluation_count: int,
) -> list[str]:
    warnings: list[str] = []
    for parameter_name, (lower_bound, upper_bound) in bounds_by_name.items():
        parameter_value = float(parameter_values[parameter_name])
        if math.isclose(
            parameter_value,
            lower_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ):
            warnings.append(f"{parameter_name} hit lower search boundary")
        if math.isclose(
            parameter_value,
            upper_bound,
            rel_tol=0.0,
            abs_tol=_BOUNDARY_TOLERANCE,
        ):
            warnings.append(f"{parameter_name} hit upper search boundary")
    if not converged:
        warnings.append(
            f"{specification.model_name} model search did not converge after {optimization_pass_count} coordinate passes"
        )
    if function_evaluation_count < 1:
        warnings.append(
            f"{specification.model_name} model search did not evaluate any candidate surface"
        )
    return warnings


def _successful_selection_row(
    *,
    specification: _CandidateSpecification,
    parameter_count: int,
    log_likelihood: float,
    parameter_values: dict[str, float],
    warnings: list[str],
) -> SubstitutionModelSelectionRow:
    aic = _compute_aic(log_likelihood, parameter_count=parameter_count)
    return SubstitutionModelSelectionRow(
        model_name=specification.model_name,
        base_model_name=specification.base_model_name,
        rate_heterogeneity_model=specification.rate_heterogeneity_model,
        fit_succeeded=True,
        parameter_count=parameter_count,
        log_likelihood=log_likelihood,
        aic=aic,
        aicc=None,
        bic=None,
        delta_aic=None,
        akaike_weight=None,
        rank=None,
        comparable_on_aic=True,
        comparable_on_aicc=True,
        comparable_on_bic=True,
        selected_by_aic=False,
        selected_by_aicc=False,
        selected_by_bic=False,
        parameter_values=parameter_values,
        warnings=warnings,
    )


def _rank_substitution_model_selection_rows(
    rows: list[SubstitutionModelSelectionRow],
    *,
    sample_size: int,
) -> list[str]:
    warnings: list[str] = []
    successful_rows = [row for row in rows if row.fit_succeeded]
    for row in successful_rows:
        row.aicc = _compute_aicc(
            row.aic if row.aic is not None else math.inf,
            sample_size=sample_size,
            parameter_count=row.parameter_count or 0,
        )
        row.bic = _compute_bic(
            row.log_likelihood if row.log_likelihood is not None else float("nan"),
            sample_size=sample_size,
            parameter_count=row.parameter_count or 0,
        )
        if row.aicc is None:
            row.comparable_on_aicc = False
            row.warnings.append(
                "sample size is too small to compute finite AICc for this parameter count"
            )
    comparable_aic_rows = [row for row in successful_rows if row.aic is not None]
    if comparable_aic_rows:
        best_aic = min(row.aic for row in comparable_aic_rows if row.aic is not None)
        ranked_aic_rows = sorted(
            comparable_aic_rows,
            key=lambda row: (
                row.aic if row.aic is not None else math.inf,
                row.model_name,
            ),
        )
        raw_weights = [
            math.exp(-0.5 * ((row.aic or math.inf) - best_aic))
            for row in ranked_aic_rows
        ]
        weight_total = sum(raw_weights)
        for rank, (row, raw_weight) in enumerate(
            zip(ranked_aic_rows, raw_weights, strict=True),
            start=1,
        ):
            row.rank = rank
            row.delta_aic = (row.aic or best_aic) - best_aic
            row.akaike_weight = raw_weight / weight_total if weight_total else 0.0
            row.selected_by_aic = math.isclose(
                row.delta_aic,
                0.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
    comparable_aicc_rows = [row for row in successful_rows if row.aicc is not None]
    if comparable_aicc_rows:
        best_aicc = min(
            row.aicc for row in comparable_aicc_rows if row.aicc is not None
        )
        for row in comparable_aicc_rows:
            row.selected_by_aicc = math.isclose(
                (row.aicc or best_aicc) - best_aicc,
                0.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
    else:
        warnings.append(
            "no finite AICc model remained available for substitution-model selection"
        )
    comparable_bic_rows = [row for row in successful_rows if row.bic is not None]
    if comparable_bic_rows:
        best_bic = min(row.bic for row in comparable_bic_rows if row.bic is not None)
        for row in comparable_bic_rows:
            row.selected_by_bic = math.isclose(
                (row.bic or best_bic) - best_bic,
                0.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
    failed_rows = [row.model_name for row in rows if not row.fit_succeeded]
    if failed_rows:
        warnings.append(
            "one or more substitution candidates failed but remain visible in the table: "
            + ", ".join(failed_rows)
        )
    rows.sort(
        key=lambda row: (
            row.rank is None,
            math.inf if row.rank is None else row.rank,
            row.model_name,
        )
    )
    return warnings


def _compute_aic(log_likelihood: float, *, parameter_count: int) -> float:
    return (2.0 * parameter_count) - (2.0 * log_likelihood)


def _compute_aicc(
    aic: float,
    *,
    sample_size: int,
    parameter_count: int,
) -> float | None:
    denominator = sample_size - parameter_count - 1
    if denominator <= 0:
        return None
    return aic + ((2.0 * parameter_count * (parameter_count + 1)) / denominator)


def _compute_bic(
    log_likelihood: float,
    *,
    sample_size: int,
    parameter_count: int,
) -> float:
    return (math.log(sample_size) * parameter_count) - (2.0 * log_likelihood)


def _base_frequency_parameter_values(
    stationary: numpy.ndarray,
) -> dict[str, float]:
    return {
        "base_frequency_a": float(stationary[0]),
        "base_frequency_c": float(stationary[1]),
        "base_frequency_g": float(stationary[2]),
        "base_frequency_t": float(stationary[3]),
    }


def _base_model_parameter_count(base_model_name: str) -> int:
    return {
        "JC69": 0,
        "K80": 1,
        "F81": 3,
        "HKY85": 4,
        "GTR": 8,
    }[base_model_name]
