from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
    empirical_protein_rate_matrix,
    evaluate_codon_ctmc_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma,
    evaluate_protein_poisson_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    compress_codon_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    build_codon_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_named_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .contracts import (
    GammaInvariantMixtureFamilyCoverageEntry,
    GammaInvariantMixtureFamilyCoverageReport,
    GammaInvariantMixtureModelCoverageCase,
)

_ALPHA = 0.8
_CATEGORY_COUNT = 4
_INVARIANT_PROPORTION = 0.35
_UPPER_INVARIANT_PROPORTION_BOUND = 0.5
_LIKELIHOOD_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class _GammaInvariantModelCase:
    tree: PhyloTree
    compressed_patterns: CompressedAlignmentSitePatterns
    surface_model_name: str
    gamma_report: object
    combined_report: (
        FiniteStateDiscreteGammaInvariantTreeLikelihoodReport
        | ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport
    )


@dataclass(frozen=True, slots=True)
class _GammaInvariantModelSpecification:
    model_name: str
    notes: list[str]
    case_builder: Callable[[], _GammaInvariantModelCase]


@dataclass(frozen=True, slots=True)
class _GammaInvariantFamilySpecification:
    family_name: str
    notes: list[str]
    models: tuple[_GammaInvariantModelSpecification, ...]


def build_gamma_invariant_family_coverage() -> (
    GammaInvariantMixtureFamilyCoverageReport
):
    """Prove that eligible model families keep gamma and invariant components active together."""

    families = [_evaluate_family_case(specification) for specification in _catalog()]
    passed_family_count = sum(1 for family in families if family.passed)
    model_count = sum(family.model_count for family in families)
    passed_model_count = sum(family.passed_model_count for family in families)
    return GammaInvariantMixtureFamilyCoverageReport(
        goal_id=30,
        passed=all(entry.passed for entry in families),
        family_count=len(families),
        model_count=model_count,
        passed_family_count=passed_family_count,
        passed_model_count=passed_model_count,
        families=families,
        limitations=[
            "the coverage report proves combined gamma and invariant activation on one governed fixture for each eligible finite-state model family; it is not a large-dataset optimizer study",
            "the report exercises the public family-level likelihood surfaces so durable coverage is tied to the owned API, not just the shared internal finite-state engine",
            "the fixtures are intentionally small and deterministic because the review target is combined +G+I activation, not empirical model-selection behavior on production alignments",
        ],
    )


def _catalog() -> tuple[_GammaInvariantFamilySpecification, ...]:
    return (
        _GammaInvariantFamilySpecification(
            family_name="nucleotide",
            notes=[
                "covers the owned nucleotide model catalogue through the selected-nucleotide likelihood surface",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "JC69",
                    [
                        "uses a constant DNA pattern so both gamma and invariant mass stay identifiable"
                    ],
                    _build_jc69_case,
                ),
                _GammaInvariantModelSpecification(
                    "K80",
                    [
                        "keeps transition bias active while proving combined +G+I handling"
                    ],
                    _build_k80_case,
                ),
                _GammaInvariantModelSpecification(
                    "F81",
                    [
                        "uses unequal base frequencies while testing combined rate heterogeneity"
                    ],
                    _build_f81_case,
                ),
                _GammaInvariantModelSpecification(
                    "HKY85",
                    ["covers simultaneous kappa, unequal frequencies, and +G+I"],
                    _build_hky85_case,
                ),
                _GammaInvariantModelSpecification(
                    "GTR",
                    [
                        "covers the most general native nucleotide surface under combined +G+I"
                    ],
                    _build_gtr_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="protein-poisson",
            notes=[
                "checks the native protein Poisson likelihood surface separately from empirical matrix models",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "protein Poisson",
                    [
                        "uses one constant amino-acid pattern so both mixture components stay measurable"
                    ],
                    _build_protein_poisson_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="protein-empirical",
            notes=[
                "covers the named empirical matrix catalogue through the empirical protein +G+I wrapper",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "JTT",
                    ["covers the JTT empirical matrix under combined +G+I"],
                    _build_protein_empirical_jtt_case,
                ),
                _GammaInvariantModelSpecification(
                    "WAG",
                    ["covers the WAG empirical matrix under combined +G+I"],
                    _build_protein_empirical_wag_case,
                ),
                _GammaInvariantModelSpecification(
                    "LG",
                    ["covers the LG empirical matrix under combined +G+I"],
                    _build_protein_empirical_lg_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="codon-ctmc",
            notes=[
                "covers the native 61-state codon surface with combined rate heterogeneity and invariant mass",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "codon CTMC",
                    [
                        "uses a constant coding pattern so codon-state handling remains inside the combined likelihood"
                    ],
                    _build_codon_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="discrete-mk-unordered",
            notes=[
                "proves combined +G+I coverage for unordered ER, SYM, and ARD discrete Mk surfaces",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "equal-rates",
                    ["covers unordered ER under combined +G+I"],
                    _build_discrete_unordered_equal_rates_case,
                ),
                _GammaInvariantModelSpecification(
                    "symmetric",
                    ["covers unordered SYM under combined +G+I"],
                    _build_discrete_unordered_symmetric_case,
                ),
                _GammaInvariantModelSpecification(
                    "all-rates-different",
                    ["covers unordered ARD under combined +G+I"],
                    _build_discrete_unordered_ard_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="discrete-mk-ordered",
            notes=[
                "proves combined +G+I coverage for ordered ER, SYM, and ARD discrete Mk surfaces",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "equal-rates",
                    ["covers ordered ER with adjacent transitions only"],
                    _build_discrete_ordered_equal_rates_case,
                ),
                _GammaInvariantModelSpecification(
                    "symmetric",
                    ["covers ordered SYM with adjacent-state symmetry"],
                    _build_discrete_ordered_symmetric_case,
                ),
                _GammaInvariantModelSpecification(
                    "all-rates-different",
                    ["covers ordered ARD with per-direction adjacent rates"],
                    _build_discrete_ordered_ard_case,
                ),
            ),
        ),
        _GammaInvariantFamilySpecification(
            family_name="discrete-mk-irreversible",
            notes=[
                "covers the one-way irreversible discrete surface under combined +G+I",
            ],
            models=(
                _GammaInvariantModelSpecification(
                    "irreversible",
                    [
                        "uses one one-way loss surface so reverse transitions remain structurally blocked"
                    ],
                    _build_discrete_irreversible_case,
                ),
            ),
        ),
    )


def _evaluate_family_case(
    specification: _GammaInvariantFamilySpecification,
) -> GammaInvariantMixtureFamilyCoverageEntry:
    model_cases = [
        _evaluate_model_case(model_specification)
        for model_specification in specification.models
    ]
    passed_model_count = sum(1 for case in model_cases if case.passed)
    return GammaInvariantMixtureFamilyCoverageEntry(
        family_name=specification.family_name,
        supported_model_names=[
            model_specification.model_name
            for model_specification in specification.models
        ],
        covered_model_names=[case.model_name for case in model_cases if case.passed],
        model_count=len(model_cases),
        passed_model_count=passed_model_count,
        passed=(passed_model_count == len(model_cases)),
        model_cases=model_cases,
        notes=specification.notes,
    )


def _evaluate_model_case(
    specification: _GammaInvariantModelSpecification,
) -> GammaInvariantMixtureModelCoverageCase:
    case = specification.case_builder()
    first_site = case.combined_report.site_likelihoods[0]
    likelihood_delta = (
        case.combined_report.log_likelihood - case.gamma_report.log_likelihood
    )
    likelihood_changed = not math.isclose(
        case.combined_report.log_likelihood,
        case.gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=_LIKELIHOOD_TOLERANCE,
    )
    invariant_component_active = first_site.invariant_component_likelihood > 0.0
    variable_component_active = first_site.variable_component_likelihood > 0.0
    category_likelihood_count = len(first_site.category_likelihoods)
    boundary_warning_count = len(case.combined_report.boundary_warnings)
    boundary_warning_emitted = (
        case.combined_report.hit_lower_invariant_proportion_boundary
        and boundary_warning_count > 0
    )
    passed = (
        likelihood_changed
        and invariant_component_active
        and variable_component_active
        and category_likelihood_count == _CATEGORY_COUNT
        and boundary_warning_emitted
    )
    return GammaInvariantMixtureModelCoverageCase(
        model_name=specification.model_name,
        surface_model_name=case.surface_model_name,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
        invariant_proportion=case.combined_report.invariant_proportion,
        site_count=case.combined_report.site_count,
        pattern_count=case.combined_report.pattern_count,
        gamma_log_likelihood=case.gamma_report.log_likelihood,
        gamma_invariant_log_likelihood=case.combined_report.log_likelihood,
        likelihood_delta=likelihood_delta,
        invariant_component_likelihood=first_site.invariant_component_likelihood,
        variable_component_likelihood=first_site.variable_component_likelihood,
        category_likelihood_count=category_likelihood_count,
        likelihood_changed=likelihood_changed,
        invariant_component_active=invariant_component_active,
        variable_component_active=variable_component_active,
        boundary_warning_count=boundary_warning_count,
        hit_lower_invariant_proportion_boundary=(
            case.combined_report.hit_lower_invariant_proportion_boundary
        ),
        hit_upper_invariant_proportion_boundary=(
            case.combined_report.hit_upper_invariant_proportion_boundary
        ),
        boundary_warnings=list(case.combined_report.boundary_warnings),
        passed=passed,
        notes=specification.notes,
    )


def _shared_tree() -> PhyloTree:
    return loads_newick("(A:0.1,B:0.2);")


def _dna_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="A", sequence="AAAA"),
        AlignmentRecord(identifier="B", sequence="AAAA"),
    ]


def _protein_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="A", sequence="AAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAA"),
    ]


def _codon_records() -> list[AlignmentRecord]:
    return [
        AlignmentRecord(identifier="A", sequence="AAAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAAA"),
    ]


def _build_jc69_case() -> _GammaInvariantModelCase:
    records = _dna_records()
    return _build_selected_nucleotide_case(records, model_name="jc69")


def _build_k80_case() -> _GammaInvariantModelCase:
    records = _dna_records()
    return _build_selected_nucleotide_case(records, model_name="k80", kappa=3.2)


def _build_f81_case() -> _GammaInvariantModelCase:
    records = _dna_records()
    return _build_selected_nucleotide_case(
        records,
        model_name="f81",
        base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
    )


def _build_hky85_case() -> _GammaInvariantModelCase:
    records = _dna_records()
    return _build_selected_nucleotide_case(
        records,
        model_name="hky85",
        kappa=4.0,
        base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
    )


def _build_gtr_case() -> _GammaInvariantModelCase:
    records = _dna_records()
    return _build_selected_nucleotide_case(
        records,
        model_name="gtr",
        base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
        exchangeabilities={
            ("A", "C"): 0.4,
            ("A", "G"): 2.0,
            ("A", "T"): 0.6,
            ("C", "G"): 0.5,
            ("C", "T"): 1.6,
            ("G", "T"): 0.7,
        },
    )


def _build_selected_nucleotide_case(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | None = None,
    exchangeabilities: dict[tuple[str, str], float] | None = None,
) -> _GammaInvariantModelCase:
    tree = _shared_tree()
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name="gamma invariant family coverage",
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
    )
    combined_report = evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        model_name=model_name,
        owner_name="gamma invariant family coverage",
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
        invariant_proportion=_INVARIANT_PROPORTION,
        lower_invariant_proportion_bound=_INVARIANT_PROPORTION,
        upper_invariant_proportion_bound=_UPPER_INVARIANT_PROPORTION_BOUND,
        emit_boundary_warnings=True,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )
    return _GammaInvariantModelCase(
        tree=tree,
        compressed_patterns=compressed_patterns,
        surface_model_name=surface.model_name,
        gamma_report=gamma_report,
        combined_report=combined_report,
    )


def _build_protein_poisson_case() -> _GammaInvariantModelCase:
    tree = _shared_tree()
    records = _protein_records()
    surface = build_protein_poisson_finite_state_surface(records)
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
    )
    combined_report = evaluate_protein_poisson_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
        invariant_proportion=_INVARIANT_PROPORTION,
        lower_invariant_proportion_bound=_INVARIANT_PROPORTION,
        upper_invariant_proportion_bound=_UPPER_INVARIANT_PROPORTION_BOUND,
        emit_boundary_warnings=True,
    )
    return _GammaInvariantModelCase(
        tree=tree,
        compressed_patterns=compressed_patterns,
        surface_model_name=surface.model_name,
        gamma_report=gamma_report,
        combined_report=combined_report,
    )


def _build_protein_empirical_jtt_case() -> _GammaInvariantModelCase:
    return _build_named_empirical_protein_case("JTT")


def _build_protein_empirical_wag_case() -> _GammaInvariantModelCase:
    return _build_named_empirical_protein_case("WAG")


def _build_protein_empirical_lg_case() -> _GammaInvariantModelCase:
    return _build_named_empirical_protein_case("LG")


def _build_named_empirical_protein_case(model_name: str) -> _GammaInvariantModelCase:
    tree = _shared_tree()
    records = _protein_records()
    surface = build_named_empirical_protein_finite_state_surface(
        records, model_name=model_name
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
    )
    combined_report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        rate_matrix=empirical_protein_rate_matrix(model_name),
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
        invariant_proportion=_INVARIANT_PROPORTION,
        lower_invariant_proportion_bound=_INVARIANT_PROPORTION,
        upper_invariant_proportion_bound=_UPPER_INVARIANT_PROPORTION_BOUND,
        emit_boundary_warnings=True,
        matrix_label=model_name,
    )
    return _GammaInvariantModelCase(
        tree=tree,
        compressed_patterns=compressed_patterns,
        surface_model_name=surface.model_name,
        gamma_report=gamma_report,
        combined_report=combined_report,
    )


def _build_codon_case() -> _GammaInvariantModelCase:
    tree = _shared_tree()
    records = _codon_records()
    surface = build_codon_finite_state_surface(records)
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
    )
    combined_report = (
        evaluate_codon_ctmc_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            alpha=_ALPHA,
            category_count=_CATEGORY_COUNT,
            invariant_proportion=_INVARIANT_PROPORTION,
            lower_invariant_proportion_bound=_INVARIANT_PROPORTION,
            upper_invariant_proportion_bound=_UPPER_INVARIANT_PROPORTION_BOUND,
            emit_boundary_warnings=True,
        )
    )
    return _GammaInvariantModelCase(
        tree=tree,
        compressed_patterns=compressed_patterns,
        surface_model_name=surface.model_name,
        gamma_report=gamma_report,
        combined_report=combined_report,
    )


def _build_discrete_unordered_equal_rates_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="unordered discrete Mk",
        state_order=("north", "south", "west"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8], dtype=float)),
            state_order=["north", "south", "west"],
            model="equal-rates",
            state_ordering="unordered",
            allowed_transition_pairs={
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            },
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_unordered_symmetric_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="unordered discrete Mk",
        state_order=("north", "south", "west"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8, 0.35, 1.1], dtype=float)),
            state_order=["north", "south", "west"],
            model="symmetric",
            state_ordering="unordered",
            allowed_transition_pairs={
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            },
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_unordered_ard_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="unordered discrete Mk",
        state_order=("north", "south", "west"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8, 0.35, 1.1, 0.7, 0.55, 0.9], dtype=float)),
            state_order=["north", "south", "west"],
            model="all-rates-different",
            state_ordering="unordered",
            allowed_transition_pairs={
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            },
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_ordered_equal_rates_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="ordered discrete Mk",
        state_order=("north", "south", "island"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8], dtype=float)),
            state_order=["north", "south", "island"],
            model="equal-rates",
            state_ordering="ordered",
            allowed_transition_pairs={(0, 1), (1, 0), (1, 2), (2, 1)},
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_ordered_symmetric_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="ordered discrete Mk",
        state_order=("north", "south", "island"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8, 0.35], dtype=float)),
            state_order=["north", "south", "island"],
            model="symmetric",
            state_ordering="ordered",
            allowed_transition_pairs={(0, 1), (1, 0), (1, 2), (2, 1)},
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_ordered_ard_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="ordered discrete Mk",
        state_order=("north", "south", "island"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8, 0.35, 1.1, 0.7], dtype=float)),
            state_order=["north", "south", "island"],
            model="all-rates-different",
            state_ordering="ordered",
            allowed_transition_pairs={(0, 1), (1, 0), (1, 2), (2, 1)},
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="north"),
            AlignmentRecord(identifier="B", sequence="north"),
        ],
    )


def _build_discrete_irreversible_case() -> _GammaInvariantModelCase:
    return _build_discrete_case(
        model_name="irreversible discrete Mk",
        state_order=("present", "lost"),
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8], dtype=float)),
            state_order=["present", "lost"],
            model="all-rates-different",
            state_ordering="ordered",
            allowed_transition_pairs={(0, 1)},
        ),
        records=[
            AlignmentRecord(identifier="A", sequence="present"),
            AlignmentRecord(identifier="B", sequence="present"),
        ],
    )


def _build_discrete_case(
    *,
    model_name: str,
    state_order: tuple[str, ...],
    rate_matrix: numpy.ndarray,
    records: list[AlignmentRecord],
) -> _GammaInvariantModelCase:
    tree = _shared_tree()
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=1,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=(records[0].sequence, records[1].sequence),
                weight=1,
                site_positions=[1],
            )
        ],
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        model_name=model_name,
    )
    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=_ALPHA,
        category_count=_CATEGORY_COUNT,
    )
    combined_report = (
        evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            state_order=state_order,
            rate_matrix=rate_matrix,
            model_name=model_name,
            alpha=_ALPHA,
            category_count=_CATEGORY_COUNT,
            invariant_proportion=_INVARIANT_PROPORTION,
            lower_invariant_proportion_bound=_INVARIANT_PROPORTION,
            upper_invariant_proportion_bound=_UPPER_INVARIANT_PROPORTION_BOUND,
            emit_boundary_warnings=True,
        )
    )
    return _GammaInvariantModelCase(
        tree=tree,
        compressed_patterns=compressed_patterns,
        surface_model_name=surface.model_name,
        gamma_report=gamma_report,
        combined_report=combined_report,
    )
