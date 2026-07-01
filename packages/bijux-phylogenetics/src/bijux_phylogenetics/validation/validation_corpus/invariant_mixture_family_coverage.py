from __future__ import annotations

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

from .contracts import (
    InvariantMixtureFamilyCoverageEntry,
    InvariantMixtureFamilyCoverageReport,
)

_INVARIANT_PROPORTION = 0.35
_LIKELIHOOD_TOLERANCE = 1e-12


def build_invariant_mixture_family_coverage() -> (
    InvariantMixtureFamilyCoverageReport
):
    """Prove that eligible finite-state families keep invariant mass inside scoring."""

    families = [
        _evaluate_family_case(
            family_name="nucleotide-jc69",
            notes=[
                "uses one constant DNA pattern on a two-tip tree so invariant mass must remain available inside the likelihood",
                "routes through the shared nucleotide finite-state surface instead of a family-specific invariant shortcut",
            ],
            case_builder=_build_nucleotide_case,
        ),
        _evaluate_family_case(
            family_name="protein-poisson",
            notes=[
                "uses one constant amino-acid pattern so the invariant component remains identifiable on the shared protein Poisson surface",
                "checks the native protein Poisson family separately from the empirical matrix wrapper surface",
            ],
            case_builder=_build_protein_poisson_case,
        ),
        _evaluate_family_case(
            family_name="protein-empirical-wag",
            notes=[
                "uses the named empirical protein surface so invariant mass is proven on the same finite-state engine that underlies WAG-style scoring",
                "keeps the test independent from the empirical protein optimizer surfaces by evaluating the shared compiled likelihood directly",
            ],
            case_builder=_build_protein_empirical_case,
        ),
        _evaluate_family_case(
            family_name="codon-ctmc",
            notes=[
                "uses a constant coding pattern so the 61-state codon surface must preserve invariant mass without collapsing to nucleotide scoring",
                "covers the codon finite-state surface that enforces codon observation and root-prior rules",
            ],
            case_builder=_build_codon_case,
        ),
        _evaluate_family_case(
            family_name="discrete-mk-unordered",
            notes=[
                "uses one constant unordered Mk pattern so the discrete invariant component is tested on the same finite-state pruning engine as continuous parity fixtures depend on",
                "keeps the discrete rate matrix explicit so the invariant signal is not confounded with a hidden ER default",
            ],
            case_builder=_build_discrete_unordered_case,
        ),
    ]
    return InvariantMixtureFamilyCoverageReport(
        goal_id=29,
        passed=all(entry.passed for entry in families),
        family_count=len(families),
        families=families,
        limitations=[
            "the coverage report proves invariant mass changes the likelihood on one invariant-rich fixture for each eligible family; it is not a whole-repository parity or optimizer report",
            "the report exercises the owned shared finite-state likelihood surfaces, so wrapper-specific ergonomics outside those surfaces are out of scope here",
            "the fixtures are intentionally small and deterministic because the review target is invariant-mixture activation, not large-dataset model selection behavior",
        ],
    )


def _evaluate_family_case(
    *,
    family_name: str,
    notes: list[str],
    case_builder,
) -> InvariantMixtureFamilyCoverageEntry:
    tree, compressed_patterns, surface = case_builder()
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
    likelihood_delta = (
        invariant_report.log_likelihood - fixed_rate_log_likelihood
    )
    likelihood_changed = not math.isclose(
        invariant_report.log_likelihood,
        fixed_rate_log_likelihood,
        rel_tol=0.0,
        abs_tol=_LIKELIHOOD_TOLERANCE,
    )
    invariant_component_active = first_site.invariant_component_likelihood > 0.0
    return InvariantMixtureFamilyCoverageEntry(
        family_name=family_name,
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
        notes=notes,
    )


def _build_nucleotide_case():
    tree = loads_newick("(A:0.1,B:0.2);")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAA"),
        AlignmentRecord(identifier="B", sequence="AAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="jc69",
        owner_name="invariant mixture family coverage",
    )
    return tree, compressed_patterns, surface


def _build_protein_poisson_case():
    tree = loads_newick("(A:0.1,B:0.2);")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)
    return tree, compressed_patterns, surface


def _build_protein_empirical_case():
    tree = loads_newick("(A:0.1,B:0.2);")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_named_empirical_protein_finite_state_surface(
        records,
        model_name="WAG",
    )
    return tree, compressed_patterns, surface


def _build_codon_case():
    tree = loads_newick("(A:0.1,B:0.2);")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAAA"),
    ]
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    surface = build_codon_finite_state_surface(records)
    return tree, compressed_patterns, surface


def _build_discrete_unordered_case():
    tree = loads_newick("(A:0.1,B:0.2);")
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=3,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("north", "north"),
                weight=3,
                site_positions=[1, 2, 3],
            )
        ],
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=["north", "south", "west"],
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
        root_prior=numpy.full(3, 1.0 / 3.0, dtype=float),
        model_name="unordered discrete Mk",
    )
    return tree, compressed_patterns, surface
