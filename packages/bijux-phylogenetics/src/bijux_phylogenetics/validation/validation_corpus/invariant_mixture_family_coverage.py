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
    evaluate_finite_state_tree_likelihood_with_invariant_mixture,
    evaluate_finite_state_tree_log_likelihood,
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
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_named_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .contracts import (
    InvariantMixtureFamilyCoverageEntry,
    InvariantMixtureFamilyCoverageReport,
    InvariantMixtureModelCoverageCase,
)

_INVARIANT_PROPORTION = 0.35
_LIKELIHOOD_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class _InvariantModelSpecification:
    model_name: str
    notes: list[str]
    case_builder: Callable[
        [], tuple[PhyloTree, CompressedAlignmentSitePatterns, object]
    ]


@dataclass(frozen=True, slots=True)
class _InvariantFamilySpecification:
    family_name: str
    notes: list[str]
    models: tuple[_InvariantModelSpecification, ...]


def build_invariant_mixture_family_coverage() -> InvariantMixtureFamilyCoverageReport:
    """Prove that eligible finite-state models keep invariant mass inside scoring."""

    families = [_evaluate_family_case(specification) for specification in _catalog()]
    passed_family_count = sum(1 for family in families if family.passed)
    model_count = sum(family.model_count for family in families)
    passed_model_count = sum(family.passed_model_count for family in families)
    return InvariantMixtureFamilyCoverageReport(
        goal_id=29,
        passed=all(entry.passed for entry in families),
        family_count=len(families),
        model_count=model_count,
        passed_family_count=passed_family_count,
        passed_model_count=passed_model_count,
        families=families,
        limitations=[
            "the coverage report proves invariant mass changes the likelihood on one invariant-rich fixture for each eligible finite-state model; it is not a whole-repository parity or optimizer report",
            "the report exercises the owned shared finite-state likelihood surfaces, so wrapper-specific ergonomics outside those surfaces are out of scope here",
            "the fixtures are intentionally small and deterministic because the review target is invariant-mixture activation, not large-dataset model selection behavior",
        ],
    )


def _catalog() -> tuple[_InvariantFamilySpecification, ...]:
    return (
        _InvariantFamilySpecification(
            family_name="nucleotide",
            notes=[
                "covers the owned nucleotide likelihood catalogue instead of using JC69 as a stand-in for every DNA model",
                "keeps kappa, base-frequency, and exchangeability-sensitive models inside the same shared invariant-mixture proof ledger",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="JC69",
                    notes=[
                        "uses one constant DNA pattern on a two-tip tree so invariant mass must remain available inside the likelihood",
                    ],
                    case_builder=_build_jc69_case,
                ),
                _InvariantModelSpecification(
                    model_name="K80",
                    notes=[
                        "exercises the transition-biased DNA surface with an explicit kappa value instead of falling back to equal-rate JC69 behavior",
                    ],
                    case_builder=_build_k80_case,
                ),
                _InvariantModelSpecification(
                    model_name="F81",
                    notes=[
                        "uses unequal base frequencies so invariant activation is proven on a composition-biased nucleotide surface",
                    ],
                    case_builder=_build_f81_case,
                ),
                _InvariantModelSpecification(
                    model_name="HKY85",
                    notes=[
                        "keeps both kappa and unequal base frequencies active while testing invariant mass on the shared DNA finite-state engine",
                    ],
                    case_builder=_build_hky85_case,
                ),
                _InvariantModelSpecification(
                    model_name="GTR",
                    notes=[
                        "uses explicit exchangeabilities and base frequencies so the most general native nucleotide surface remains covered by the invariant-mixture proof",
                    ],
                    case_builder=_build_gtr_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="protein-poisson",
            notes=[
                "checks the native protein Poisson surface separately from empirical matrix models",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="protein Poisson",
                    notes=[
                        "uses one constant amino-acid pattern so the invariant component remains identifiable on the shared protein Poisson surface",
                    ],
                    case_builder=_build_protein_poisson_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="protein-empirical",
            notes=[
                "proves invariant activation on every named empirical matrix the repository currently exposes",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="JTT",
                    notes=[
                        "covers the JTT empirical matrix on the shared finite-state engine rather than through wrapper-only reporting surfaces",
                    ],
                    case_builder=_build_protein_empirical_jtt_case,
                ),
                _InvariantModelSpecification(
                    model_name="WAG",
                    notes=[
                        "retains the existing WAG proof while moving it into the full empirical-model catalogue",
                    ],
                    case_builder=_build_protein_empirical_wag_case,
                ),
                _InvariantModelSpecification(
                    model_name="LG",
                    notes=[
                        "extends empirical-matrix coverage to LG so the invariant proof is not limited to one named protein matrix",
                    ],
                    case_builder=_build_protein_empirical_lg_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="codon-ctmc",
            notes=[
                "covers the 61-state codon surface that enforces codon observation and root-prior rules",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="codon CTMC",
                    notes=[
                        "uses a constant coding pattern so the codon surface must preserve invariant mass without collapsing to nucleotide scoring",
                    ],
                    case_builder=_build_codon_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="discrete-mk-unordered",
            notes=[
                "proves invariant activation for the unordered ER, SYM, and ARD discrete Mk surfaces",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="equal-rates",
                    notes=[
                        "keeps the unordered ER rate matrix explicit so the invariant signal is not confounded with a hidden default",
                    ],
                    case_builder=_build_discrete_unordered_equal_rates_case,
                ),
                _InvariantModelSpecification(
                    model_name="symmetric",
                    notes=[
                        "covers the unordered SYM surface with pairwise rate heterogeneity while preserving directional symmetry",
                    ],
                    case_builder=_build_discrete_unordered_symmetric_case,
                ),
                _InvariantModelSpecification(
                    model_name="all-rates-different",
                    notes=[
                        "covers the unordered ARD surface so invariant mass remains active under fully directional multistate rates",
                    ],
                    case_builder=_build_discrete_unordered_ard_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="discrete-mk-ordered",
            notes=[
                "covers the ordered ER, SYM, and ARD surfaces where state adjacency restricts legal transitions",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="equal-rates",
                    notes=[
                        "uses one ordered ER surface with only adjacent transitions allowed",
                    ],
                    case_builder=_build_discrete_ordered_equal_rates_case,
                ),
                _InvariantModelSpecification(
                    model_name="symmetric",
                    notes=[
                        "uses ordered SYM so the invariant proof remains active under adjacent-state symmetry constraints",
                    ],
                    case_builder=_build_discrete_ordered_symmetric_case,
                ),
                _InvariantModelSpecification(
                    model_name="all-rates-different",
                    notes=[
                        "uses ordered ARD so invariant activation is covered even when each adjacent direction has its own rate",
                    ],
                    case_builder=_build_discrete_ordered_ard_case,
                ),
            ),
        ),
        _InvariantFamilySpecification(
            family_name="discrete-mk-irreversible",
            notes=[
                "covers the governed one-way irreversible discrete surface instead of assuming the unordered ARD proof is sufficient",
            ],
            models=(
                _InvariantModelSpecification(
                    model_name="irreversible",
                    notes=[
                        "uses one one-way loss surface so the invariant proof remains active when reverse transitions are structurally forbidden",
                    ],
                    case_builder=_build_discrete_irreversible_case,
                ),
            ),
        ),
    )


def _evaluate_family_case(
    specification: _InvariantFamilySpecification,
) -> InvariantMixtureFamilyCoverageEntry:
    model_cases = [
        _evaluate_model_case(model_specification)
        for model_specification in specification.models
    ]
    passed_model_count = sum(1 for case in model_cases if case.passed)
    return InvariantMixtureFamilyCoverageEntry(
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
    specification: _InvariantModelSpecification,
) -> InvariantMixtureModelCoverageCase:
    tree, compressed_patterns, surface = specification.case_builder()
    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compressed_patterns,
        surface=surface,
    )
    invariant_report = evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=_INVARIANT_PROPORTION,
    )
    first_site = invariant_report.site_likelihoods[0]
    likelihood_delta = invariant_report.log_likelihood - fixed_rate_log_likelihood
    likelihood_changed = not math.isclose(
        invariant_report.log_likelihood,
        fixed_rate_log_likelihood,
        rel_tol=0.0,
        abs_tol=_LIKELIHOOD_TOLERANCE,
    )
    invariant_component_active = first_site.invariant_component_likelihood > 0.0
    return InvariantMixtureModelCoverageCase(
        model_name=specification.model_name,
        surface_model_name=surface.model_name,
        invariant_proportion=_INVARIANT_PROPORTION,
        site_count=invariant_report.site_count,
        pattern_count=invariant_report.pattern_count,
        fixed_rate_log_likelihood=fixed_rate_log_likelihood,
        invariant_log_likelihood=invariant_report.log_likelihood,
        likelihood_delta=likelihood_delta,
        invariant_component_likelihood=first_site.invariant_component_likelihood,
        variable_component_likelihood=first_site.variable_component_likelihood,
        likelihood_changed=likelihood_changed,
        invariant_component_active=invariant_component_active,
        passed=likelihood_changed and invariant_component_active,
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


def _build_jc69_case():
    records = _dna_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_selected_nucleotide_finite_state_surface(
            records,
            model_name="jc69",
            owner_name="invariant mixture family coverage",
        ),
    )


def _build_k80_case():
    records = _dna_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_selected_nucleotide_finite_state_surface(
            records,
            model_name="k80",
            owner_name="invariant mixture family coverage",
            kappa=3.2,
        ),
    )


def _build_f81_case():
    records = _dna_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_selected_nucleotide_finite_state_surface(
            records,
            model_name="f81",
            owner_name="invariant mixture family coverage",
            base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
        ),
    )


def _build_hky85_case():
    records = _dna_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_selected_nucleotide_finite_state_surface(
            records,
            model_name="hky85",
            owner_name="invariant mixture family coverage",
            kappa=4.0,
            base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
        ),
    )


def _build_gtr_case():
    records = _dna_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_selected_nucleotide_finite_state_surface(
            records,
            model_name="gtr",
            owner_name="invariant mixture family coverage",
            base_frequencies={"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
            exchangeabilities={
                ("A", "C"): 0.4,
                ("A", "G"): 2.0,
                ("A", "T"): 0.6,
                ("C", "G"): 0.5,
                ("C", "T"): 1.6,
                ("G", "T"): 0.7,
            },
        ),
    )


def _build_protein_poisson_case():
    records = _protein_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_protein_poisson_finite_state_surface(records),
    )


def _build_protein_empirical_jtt_case():
    return _build_named_empirical_protein_case("JTT")


def _build_protein_empirical_wag_case():
    return _build_named_empirical_protein_case("WAG")


def _build_protein_empirical_lg_case():
    return _build_named_empirical_protein_case("LG")


def _build_named_empirical_protein_case(model_name: str):
    records = _protein_records()
    return (
        _shared_tree(),
        compress_alignment_site_patterns_from_records(records),
        build_named_empirical_protein_finite_state_surface(
            records,
            model_name=model_name,
        ),
    )


def _build_codon_case():
    records = _codon_records()
    return (
        _shared_tree(),
        compress_codon_site_patterns_from_records(records),
        build_codon_finite_state_surface(records),
    )


def _build_discrete_unordered_equal_rates_case():
    return _build_discrete_case(
        model_name="equal-rates",
        state_ordering="unordered",
        state_order=("north", "south", "west"),
        log_parameters=numpy.log(numpy.array([0.8], dtype=float)),
        allowed_transition_pairs=frozenset(
            {
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            }
        ),
        invariant_state="north",
    )


def _build_discrete_unordered_symmetric_case():
    return _build_discrete_case(
        model_name="symmetric",
        state_ordering="unordered",
        state_order=("north", "south", "west"),
        log_parameters=numpy.log(numpy.array([0.8, 0.45, 1.1], dtype=float)),
        allowed_transition_pairs=frozenset(
            {
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            }
        ),
        invariant_state="north",
    )


def _build_discrete_unordered_ard_case():
    return _build_discrete_case(
        model_name="all-rates-different",
        state_ordering="unordered",
        state_order=("north", "south", "west"),
        log_parameters=numpy.log(
            numpy.array([0.8, 0.35, 1.1, 0.7, 0.5, 0.9], dtype=float)
        ),
        allowed_transition_pairs=frozenset(
            {
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            }
        ),
        invariant_state="north",
    )


def _build_discrete_ordered_equal_rates_case():
    return _build_discrete_case(
        model_name="equal-rates",
        state_ordering="ordered",
        state_order=("north", "south", "island"),
        log_parameters=numpy.log(numpy.array([0.8], dtype=float)),
        allowed_transition_pairs=frozenset({(0, 1), (1, 0), (1, 2), (2, 1)}),
        invariant_state="north",
    )


def _build_discrete_ordered_symmetric_case():
    return _build_discrete_case(
        model_name="symmetric",
        state_ordering="ordered",
        state_order=("north", "south", "island"),
        log_parameters=numpy.log(numpy.array([0.8, 0.45], dtype=float)),
        allowed_transition_pairs=frozenset({(0, 1), (1, 0), (1, 2), (2, 1)}),
        invariant_state="north",
    )


def _build_discrete_ordered_ard_case():
    return _build_discrete_case(
        model_name="all-rates-different",
        state_ordering="ordered",
        state_order=("north", "south", "island"),
        log_parameters=numpy.log(numpy.array([0.8, 0.35, 1.1, 0.7], dtype=float)),
        allowed_transition_pairs=frozenset({(0, 1), (1, 0), (1, 2), (2, 1)}),
        invariant_state="north",
    )


def _build_discrete_irreversible_case():
    return _build_discrete_case(
        model_name="all-rates-different",
        state_ordering="ordered",
        state_order=("present", "lost"),
        log_parameters=numpy.log(numpy.array([0.8], dtype=float)),
        allowed_transition_pairs=frozenset({(0, 1)}),
        invariant_state="present",
        surface_model_name="irreversible discrete Mk",
    )


def _build_discrete_case(
    *,
    model_name: str,
    state_ordering: str,
    state_order: tuple[str, ...],
    log_parameters: numpy.ndarray,
    allowed_transition_pairs: frozenset[tuple[int, int]],
    invariant_state: str,
    surface_model_name: str | None = None,
):
    tree = _shared_tree()
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=3,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=(invariant_state, invariant_state),
                weight=3,
                site_positions=[1, 2, 3],
            )
        ],
    )
    rate_matrix = rate_matrix_from_log_parameters(
        log_parameters,
        state_order=list(state_order),
        model=model_name,
        state_ordering=state_ordering,
        allowed_transition_pairs=set(allowed_transition_pairs),
    )
    model_label = (
        surface_model_name
        if surface_model_name is not None
        else f"{state_ordering} discrete Mk ({model_name})"
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=list(state_order),
        rate_matrix=rate_matrix,
        root_prior=numpy.full(len(state_order), 1.0 / len(state_order), dtype=float),
        model_name=model_label,
    )
    return tree, compressed_patterns, surface
