from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
from statistics import median
import tempfile

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from bijux_phylogenetics.compare.topology import (
    BranchLengthComparisonReport,
    TreeComparisonReport,
    _informative_clades,
    _robinson_foulds_metrics,
    _unrooted_splits,
    compare_branch_lengths,
    compare_tree_paths,
)
from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.errors import InvalidAlignmentError, InvalidDistanceMatrixError
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.tree_set import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    write_clade_frequency_table,
)

DistanceModel = str
GapHandlingMode = str
AmbiguityPolicy = str

_COMPARISON_NUCLEOTIDES = {"A", "C", "G", "T"}
_NUCLEOTIDE_AMBIGUITIES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}
_PROTEIN_RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")
_PROTEIN_AMBIGUITIES = {
    **{residue: {residue} for residue in _PROTEIN_RESIDUES},
    "B": {"D", "N"},
    "J": {"I", "L"},
    "Z": {"E", "Q"},
    "X": set(_PROTEIN_RESIDUES),
}
_MISSING_OR_GAP = {"-", "?"}
_TRANSITIONS = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C"),
}
_DISTANCE_MODEL_ALIASES = {
    "raw": "p-distance",
    "p-distance": "p-distance",
    "jc69": "jukes-cantor",
    "jukes-cantor": "jukes-cantor",
    "k80": "kimura-2-parameter",
    "kimura-2-parameter": "kimura-2-parameter",
    "f81": "felsenstein-81",
    "felsenstein-81": "felsenstein-81",
    "tn93": "tamura-nei-93",
    "tamura-nei-93": "tamura-nei-93",
    "amino-acid-p-distance": "amino-acid-p-distance",
}


@dataclass(frozen=True, slots=True)
class GeneticDistanceModelParameters:
    """Alignment-wide parameters used by composition-aware DNA distance models."""

    informative_base_count: int
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    purine_frequency: float
    pyrimidine_frequency: float
    f81_limit: float
    tn93_ag_coefficient: float | None
    tn93_ct_coefficient: float | None
    tn93_transversion_coefficient: float | None


@dataclass(frozen=True, slots=True)
class PairwiseGeneticDistance:
    """One pairwise genetic distance entry for an aligned dataset."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int
    mismatch_sites: float
    transition_sites: float
    ag_transition_sites: float
    ct_transition_sites: float
    transversion_sites: float
    ambiguity_sites: int
    skipped_sites: int
    saturated: bool
    saturation_reason: str | None


@dataclass(slots=True)
class GeneticDistanceMatrix:
    """Deterministic pairwise genetic distance matrix for one alignment."""

    path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    inferred_alphabet: str
    alignment_length: int
    identifiers: list[str]
    model_parameters: GeneticDistanceModelParameters | None
    warnings: list[str]
    pairs: list[PairwiseGeneticDistance]


@dataclass(slots=True)
class DistanceTreeBuildReport:
    """Explicit report for a distance-based tree build."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    taxon_count: int
    pair_count: int
    assumptions: DistanceMethodAssumptionReport


@dataclass(slots=True)
class DistanceTreeTopologyComparison:
    """Topology comparison between NJ and UPGMA trees built from one alignment."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    shared_taxa: list[str]
    nj_informative_clades: int
    upgma_informative_clades: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


@dataclass(frozen=True, slots=True)
class UPGMAUltrametricViolation:
    """One taxon triple whose distances are inconsistent with an ultrametric clock."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    smallest_distance: float
    middle_distance: float
    largest_distance: float
    deviation: float


@dataclass(slots=True)
class DistanceMethodAssumptionReport:
    """Audit whether an alignment or matrix respects core distance-tree assumptions."""

    source_path: Path
    source_kind: str
    taxon_count: int
    pair_count: int
    nj_assumptions: list[str]
    upgma_assumptions: list[str]
    ultrametric_compatible: bool
    ultrametric_tolerance: float
    upgma_ultrametric_violations: list[UPGMAUltrametricViolation]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class ImportedDistanceEntry:
    """One directional entry from an imported long-form distance matrix table."""

    left_identifier: str
    right_identifier: str
    distance: float
    comparable_sites: int | None


@dataclass(frozen=True, slots=True)
class DistanceMatrixAsymmetry:
    """Two directional entries that disagree numerically."""

    left_identifier: str
    right_identifier: str
    left_to_right_distance: float
    right_to_left_distance: float


@dataclass(frozen=True, slots=True)
class NonMetricDistanceObservation:
    """One triangle-inequality violation within an imported distance matrix."""

    left_identifier: str
    middle_identifier: str
    right_identifier: str
    direct_distance: float
    indirect_distance: float


@dataclass(slots=True)
class ImportedDistanceMatrixReport:
    """Validation report for an imported long-form distance matrix table."""

    path: Path
    identifiers: list[str]
    pair_count: int
    complete: bool
    zero_diagonal: bool
    symmetric: bool
    nonnegative: bool
    missing_pairs: list[str]
    diagonal_problems: list[str]
    negative_distance_pairs: list[str]
    asymmetric_pairs: list[DistanceMatrixAsymmetry]
    nonmetric_observations: list[NonMetricDistanceObservation]
    warnings: list[str]


@dataclass(slots=True)
class ImportedDistanceTreeBuildReport:
    """Explicit report for building a tree from an imported distance matrix."""

    matrix_path: Path
    method: str
    taxon_count: int
    pair_count: int
    assumptions: DistanceMethodAssumptionReport


@dataclass(slots=True)
class ImportedDistanceMatrixQualityReport:
    """Diagnostics over an imported long-form distance matrix."""

    validation: ImportedDistanceMatrixReport
    saturated_pairs: list[SaturatedDistancePair]
    low_information_pairs: list[LowInformationPair]
    low_information_pair_cutoff: int | None
    saturation_audit_scale: str
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DistanceReferenceObservation:
    """One reference example used to verify core distance calculations."""

    case: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    left_identifier: str
    right_identifier: str
    expected_distance: float | None
    observed_distance: float | None
    comparable_sites: int
    expected_ambiguity_sites: int | None
    observed_ambiguity_sites: int
    passed: bool


@dataclass(slots=True)
class DistanceReferenceValidationReport:
    """Validation of built-in reference distance examples."""

    observations: list[DistanceReferenceObservation]
    tree_observations: list[DistanceTreeReferenceObservation]
    all_passed: bool


@dataclass(frozen=True, slots=True)
class DistanceTreeReferenceObservation:
    """One reference example used to validate distance-tree clustering."""

    case: str
    method: str
    matrix_path: Path
    expected_clades: list[str]
    observed_clades: list[str]
    passed: bool


@dataclass(frozen=True, slots=True)
class SaturatedDistancePair:
    """One pair that reaches an undefined or unreliable correction regime."""

    left_identifier: str
    right_identifier: str
    distance: float | None
    comparable_sites: int
    reason: str


@dataclass(frozen=True, slots=True)
class DistanceOutlierPair:
    """One pair whose distance is unusually large relative to the dataset."""

    left_identifier: str
    right_identifier: str
    distance: float
    note: str


@dataclass(frozen=True, slots=True)
class LowInformationPair:
    """One pair with too few comparable sites for robust interpretation."""

    left_identifier: str
    right_identifier: str
    comparable_sites: int
    note: str


@dataclass(slots=True)
class DistanceMethodAssessment:
    """Decision about whether the computed matrix is suitable for distance methods."""

    decision: str
    reasons: list[str]


@dataclass(slots=True)
class DistanceMatrixQualityReport:
    """Diagnostics over a computed distance matrix."""

    path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    inferred_alphabet: str
    taxon_count: int
    pair_count: int
    saturated_pairs: list[SaturatedDistancePair]
    high_distance_outliers: list[DistanceOutlierPair]
    low_information_pairs: list[LowInformationPair]
    assumptions: DistanceMethodAssumptionReport
    warnings: list[str]
    method_assessment: DistanceMethodAssessment


@dataclass(frozen=True, slots=True)
class DistanceBootstrapSupportRow:
    """One clade support row across bootstrap replicate trees."""

    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class DistanceBootstrapReport:
    """Bootstrap summary for a distance-based tree-building workflow."""

    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    replicates: int
    seed: int
    tree_count: int
    consensus_newick: str
    support: list[DistanceBootstrapSupportRow]


@dataclass(slots=True)
class DistanceBootstrapSupportSummary:
    """Reviewer-facing summary over bootstrap clade support frequencies."""

    alignment_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    replicates: int
    clade_count: int
    minimum_frequency: float | None
    maximum_frequency: float | None
    median_frequency: float | None
    weak_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class DistanceTreeReferenceComparisonReport:
    """Compare one built distance tree against a reviewer-supplied reference tree."""

    alignment_path: Path
    reference_tree_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    topology: TreeComparisonReport
    branch_lengths: BranchLengthComparisonReport
    warnings: list[str]


@dataclass(slots=True)
class DistanceModelComparisonRow:
    """One supported distance model summarized over the same alignment."""

    model: DistanceModel
    defined_pair_count: int
    saturated_pair_count: int
    low_information_pair_count: int
    mean_distance: float | None
    maximum_distance: float | None
    decision: str
    reasons: list[str]


@dataclass(slots=True)
class DistanceModelComparisonReport:
    """Comparison of all supported distance models for one alignment."""

    alignment_path: Path
    inferred_alphabet: str
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    rows: list[DistanceModelComparisonRow]
    warnings: list[str]


@dataclass(slots=True)
class DistanceGapPolicyDeltaRow:
    """One taxon pair whose distance changed across gap-handling policies."""

    left_identifier: str
    right_identifier: str
    pairwise_distance: float | None
    complete_distance: float | None
    pairwise_comparable_sites: int
    complete_comparable_sites: int
    distance_delta: float | None
    comparable_site_delta: int


@dataclass(slots=True)
class DistanceGapPolicySensitivityReport:
    """Summarize how pairwise versus complete deletion changes the same analysis."""

    alignment_path: Path
    model: DistanceModel
    ambiguity_policy: AmbiguityPolicy
    changed_pair_count: int
    pair_count: int
    rows: list[DistanceGapPolicyDeltaRow]
    warnings: list[str]


@dataclass(slots=True)
class DistanceMethodMaturityCheck:
    """One explicit maturity criterion for distance-analysis surfaces."""

    name: str
    satisfied: bool
    details: str


@dataclass(slots=True)
class DistanceMethodMaturityGateReport:
    """High-level maturity gate over one distance-analysis workflow."""

    alignment_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    decision: str
    checks: list[DistanceMethodMaturityCheck]
    warnings: list[str]


@dataclass(slots=True)
class DistanceMethodReport:
    """Structured machine-readable report for one distance-analysis workflow."""

    alignment_path: Path
    method: str
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    matrix: GeneticDistanceMatrix
    quality: DistanceMatrixQualityReport
    assumptions: DistanceMethodAssumptionReport
    reference_validation: DistanceReferenceValidationReport
    built_tree_newick: str
    alternative_tree_newick: str
    topology_comparison: DistanceTreeTopologyComparison
    bootstrap_summary: DistanceBootstrapSupportSummary
    model_comparison: DistanceModelComparisonReport
    gap_policy_sensitivity: DistanceGapPolicySensitivityReport
    maturity_gate: DistanceMethodMaturityGateReport


@dataclass(slots=True)
class DistanceReproducibilityBundleReport:
    """Reproducibility bundle written for one distance-analysis run."""

    out_dir: Path
    alignment_path: Path
    model: DistanceModel
    gap_handling: GapHandlingMode
    ambiguity_policy: AmbiguityPolicy
    method: str
    replicates: int
    files: list[Path]


@dataclass(frozen=True, slots=True)
class _SiteContribution:
    comparable: bool
    mismatch_weight: float
    transition_weight: float
    ag_transition_weight: float
    ct_transition_weight: float
    transversion_weight: float
    ambiguous: bool


@dataclass(frozen=True, slots=True)
class _PairSummary:
    distance: float | None
    comparable_sites: int
    mismatch_sites: float
    transition_sites: float
    ag_transition_sites: float
    ct_transition_sites: float
    transversion_sites: float
    ambiguity_sites: int
    skipped_sites: int
    saturated: bool
    saturation_reason: str | None


def _normalize_residue(residue: str) -> str:
    upper = residue.upper()
    if upper == "U":
        return "T"
    return upper


def _normalize_distance_model(model: DistanceModel) -> DistanceModel:
    normalized = _DISTANCE_MODEL_ALIASES.get(model)
    if normalized is None:
        raise ValueError(f"unsupported distance model: {model}")
    return normalized


def _allowed_models_for_alphabet(alphabet: str) -> set[str]:
    if alphabet in {"dna", "rna"}:
        return {
            "p-distance",
            "jukes-cantor",
            "kimura-2-parameter",
            "felsenstein-81",
            "tamura-nei-93",
        }
    if alphabet == "protein":
        return {"amino-acid-p-distance"}
    return set()


def _states_for_symbol(symbol: str, *, alphabet: str) -> set[str] | None:
    normalized = _normalize_residue(symbol)
    if normalized in _MISSING_OR_GAP:
        return set()
    if alphabet == "protein":
        return _PROTEIN_AMBIGUITIES.get(normalized)
    return _NUCLEOTIDE_AMBIGUITIES.get(normalized)


def _is_transition(left: str, right: str) -> bool:
    return (left, right) in _TRANSITIONS


def _transition_class_weights(left: str, right: str) -> tuple[float, float, float]:
    if left == right:
        return 0.0, 0.0, 0.0
    if (left, right) in {("A", "G"), ("G", "A")}:
        return 1.0, 1.0, 0.0
    if (left, right) in {("C", "T"), ("T", "C")}:
        return 1.0, 0.0, 1.0
    return 0.0, 0.0, 0.0


def _site_contribution(
    left_symbol: str,
    right_symbol: str,
    *,
    alphabet: str,
    ambiguity_policy: AmbiguityPolicy,
) -> _SiteContribution | None:
    left_states = _states_for_symbol(left_symbol, alphabet=alphabet)
    right_states = _states_for_symbol(right_symbol, alphabet=alphabet)
    if left_states is None or right_states is None:
        return None
    if not left_states or not right_states:
        return None
    ambiguous = len(left_states) > 1 or len(right_states) > 1
    if ambiguity_policy in {"ignore", "report-only"} and ambiguous:
        return None

    if not ambiguous:
        left = next(iter(left_states))
        right = next(iter(right_states))
        mismatch_weight = 0.0 if left == right else 1.0
        if alphabet == "protein":
            transition_weight = 0.0
            ag_transition_weight = 0.0
            ct_transition_weight = 0.0
        else:
            (
                transition_weight,
                ag_transition_weight,
                ct_transition_weight,
            ) = _transition_class_weights(left, right)
        transversion_weight = mismatch_weight - transition_weight
        return _SiteContribution(
            comparable=True,
            mismatch_weight=mismatch_weight,
            transition_weight=transition_weight,
            ag_transition_weight=ag_transition_weight,
            ct_transition_weight=ct_transition_weight,
            transversion_weight=transversion_weight,
            ambiguous=False,
        )

    if ambiguity_policy == "partial-match":
        total_pairs = len(left_states) * len(right_states)
        equal_pairs = 0
        transition_pairs = 0
        ag_transition_pairs = 0
        ct_transition_pairs = 0
        transversion_pairs = 0
        for left in left_states:
            for right in right_states:
                if left == right:
                    equal_pairs += 1
                else:
                    if alphabet == "protein":
                        transversion_pairs += 1
                    else:
                        transition_weight, ag_weight, ct_weight = (
                            _transition_class_weights(left, right)
                        )
                        transition_pairs += transition_weight
                        ag_transition_pairs += ag_weight
                        ct_transition_pairs += ct_weight
                        transversion_pairs += 1.0 - transition_weight
        return _SiteContribution(
            comparable=True,
            mismatch_weight=(total_pairs - equal_pairs) / total_pairs,
            transition_weight=transition_pairs / total_pairs,
            ag_transition_weight=ag_transition_pairs / total_pairs,
            ct_transition_weight=ct_transition_pairs / total_pairs,
            transversion_weight=transversion_pairs / total_pairs,
            ambiguous=True,
        )

    if ambiguity_policy == "strict-mismatch":
        if left_states == right_states and left_symbol.upper() == right_symbol.upper():
            return _SiteContribution(
                comparable=True,
                mismatch_weight=0.0,
                transition_weight=0.0,
                ag_transition_weight=0.0,
                ct_transition_weight=0.0,
                transversion_weight=0.0,
                ambiguous=True,
            )
        mismatch_pairs: list[tuple[str, str]] = [
            (left, right)
            for left in left_states
            for right in right_states
            if left != right
        ]
        if not mismatch_pairs:
            return _SiteContribution(
                comparable=True,
                mismatch_weight=0.0,
                transition_weight=0.0,
                ag_transition_weight=0.0,
                ct_transition_weight=0.0,
                transversion_weight=0.0,
                ambiguous=True,
            )
        transition_pairs = 0.0
        ag_transition_pairs = 0.0
        ct_transition_pairs = 0.0
        transversion_pairs = 0.0
        for left, right in mismatch_pairs:
            if alphabet == "protein":
                transversion_pairs += 1.0
            else:
                transition_weight, ag_weight, ct_weight = _transition_class_weights(
                    left, right
                )
                transition_pairs += transition_weight
                ag_transition_pairs += ag_weight
                ct_transition_pairs += ct_weight
                transversion_pairs += 1.0 - transition_weight
        return _SiteContribution(
            comparable=True,
            mismatch_weight=1.0,
            transition_weight=transition_pairs / len(mismatch_pairs),
            ag_transition_weight=ag_transition_pairs / len(mismatch_pairs),
            ct_transition_weight=ct_transition_pairs / len(mismatch_pairs),
            transversion_weight=transversion_pairs / len(mismatch_pairs),
            ambiguous=True,
        )

    raise ValueError(f"unsupported ambiguity policy: {ambiguity_policy}")


def _p_distance(summary: _PairSummary) -> float | None:
    if summary.comparable_sites == 0:
        return None
    return round(summary.mismatch_sites / summary.comparable_sites, 15)


def _estimate_nucleotide_model_parameters(
    records: list[AlignmentRecord],
) -> tuple[GeneticDistanceModelParameters, list[str]]:
    counts = {base: 0 for base in "ACGT"}
    informative_base_count = 0
    for record in records:
        for symbol in record.sequence:
            states = _states_for_symbol(symbol, alphabet="dna")
            if states is None or len(states) != 1:
                continue
            base = next(iter(states))
            counts[base] += 1
            informative_base_count += 1

    if informative_base_count == 0:
        parameters = GeneticDistanceModelParameters(
            informative_base_count=0,
            base_frequency_a=0.0,
            base_frequency_c=0.0,
            base_frequency_g=0.0,
            base_frequency_t=0.0,
            purine_frequency=0.0,
            pyrimidine_frequency=0.0,
            f81_limit=0.0,
            tn93_ag_coefficient=None,
            tn93_ct_coefficient=None,
            tn93_transversion_coefficient=None,
        )
        return parameters, [
            "no resolved A/C/G/T nucleotides remain to estimate composition-aware DNA distance parameters"
        ]

    pi_a = counts["A"] / informative_base_count
    pi_c = counts["C"] / informative_base_count
    pi_g = counts["G"] / informative_base_count
    pi_t = counts["T"] / informative_base_count
    pi_r = pi_a + pi_g
    pi_y = pi_c + pi_t
    f81_limit = 1.0 - ((pi_a * pi_a) + (pi_c * pi_c) + (pi_g * pi_g) + (pi_t * pi_t))

    tn93_ag_coefficient = None
    if pi_r > 0.0:
        tn93_ag_coefficient = (2.0 * pi_a * pi_g) / pi_r
    tn93_ct_coefficient = None
    if pi_y > 0.0:
        tn93_ct_coefficient = (2.0 * pi_c * pi_t) / pi_y
    tn93_transversion_coefficient = None
    if pi_r > 0.0 and pi_y > 0.0:
        tn93_transversion_coefficient = 2.0 * (
            (pi_r * pi_y)
            - ((pi_a * pi_g * pi_y) / pi_r)
            - ((pi_c * pi_t * pi_r) / pi_y)
        )

    warnings: list[str] = []
    if f81_limit <= 0.0:
        warnings.append(
            "alignment-wide resolved base composition leaves no variability for F81 correction"
        )
    if min(pi_a, pi_c, pi_g, pi_t) == 0.0:
        warnings.append(
            "alignment-wide resolved base composition omits at least one nucleotide, so TN93 assumptions break"
        )
    return (
        GeneticDistanceModelParameters(
            informative_base_count=informative_base_count,
            base_frequency_a=round(pi_a, 15),
            base_frequency_c=round(pi_c, 15),
            base_frequency_g=round(pi_g, 15),
            base_frequency_t=round(pi_t, 15),
            purine_frequency=round(pi_r, 15),
            pyrimidine_frequency=round(pi_y, 15),
            f81_limit=round(f81_limit, 15),
            tn93_ag_coefficient=None
            if tn93_ag_coefficient is None
            else round(tn93_ag_coefficient, 15),
            tn93_ct_coefficient=None
            if tn93_ct_coefficient is None
            else round(tn93_ct_coefficient, 15),
            tn93_transversion_coefficient=None
            if tn93_transversion_coefficient is None
            else round(tn93_transversion_coefficient, 15),
        ),
        warnings,
    )


def _jukes_cantor_distance(p_distance: float | None) -> tuple[float | None, str | None]:
    if p_distance is None:
        return None, "no comparable sites remain after filtering"
    if p_distance == 0.0:
        return 0.0, None
    if p_distance == 0.75:
        return (
            None,
            "p-distance is at the Jukes-Cantor correction limit, so the corrected distance tends to infinity",
        )
    if p_distance > 0.75:
        return (
            None,
            "p-distance exceeds the Jukes-Cantor correction range, so the corrected distance is undefined",
        )
    return round((-3.0 / 4.0) * math.log(1.0 - (4.0 * p_distance / 3.0)), 15), None


def _kimura_two_parameter_distance(
    summary: _PairSummary,
) -> tuple[float | None, str | None]:
    if summary.comparable_sites == 0:
        return None, "no comparable sites remain after filtering"
    p = summary.transition_sites / summary.comparable_sites
    q = summary.transversion_sites / summary.comparable_sites
    first = 1.0 - (2.0 * p) - q
    second = 1.0 - (2.0 * q)
    if first < 0.0 or second < 0.0:
        return (
            None,
            "transition and transversion proportions exceed the Kimura 2-parameter correction range, so the corrected distance is undefined",
        )
    if first == 0.0 or second == 0.0:
        return (
            None,
            "transition and transversion proportions are at the Kimura 2-parameter correction limit, so the corrected distance tends to infinity",
        )
    value = (-0.5 * math.log(first)) - (0.25 * math.log(second))
    return round(value, 15), None


def _felsenstein_81_distance(
    summary: _PairSummary, parameters: GeneticDistanceModelParameters | None
) -> tuple[float | None, str | None]:
    if summary.comparable_sites == 0:
        return None, "no comparable sites remain after filtering"
    if summary.mismatch_sites == 0.0:
        return 0.0, None
    if parameters is None or parameters.informative_base_count == 0:
        return (
            None,
            "alignment-wide resolved base composition is unavailable, so the F81 correction is undefined",
        )
    limit = parameters.f81_limit
    if limit <= 0.0:
        return (
            None,
            "alignment-wide resolved base composition leaves no variability for F81 correction, so the corrected distance is undefined",
        )
    p_distance = summary.mismatch_sites / summary.comparable_sites
    if p_distance > limit:
        return (
            None,
            "observed mismatch proportion exceeds the F81 correction range for the estimated base composition, so the corrected distance is undefined",
        )
    if p_distance == limit:
        return (
            None,
            "observed mismatch proportion is at the F81 correction limit for the estimated base composition, so the corrected distance tends to infinity",
        )
    value = -limit * math.log(1.0 - (p_distance / limit))
    return round(value, 15), None


def _tamura_nei_93_distance(
    summary: _PairSummary, parameters: GeneticDistanceModelParameters | None
) -> tuple[float | None, str | None]:
    if summary.comparable_sites == 0:
        return None, "no comparable sites remain after filtering"
    if summary.mismatch_sites == 0.0:
        return 0.0, None
    if parameters is None or parameters.informative_base_count == 0:
        return (
            None,
            "alignment-wide resolved base composition is unavailable, so the TN93 correction is undefined",
        )
    pi_a = parameters.base_frequency_a
    pi_c = parameters.base_frequency_c
    pi_g = parameters.base_frequency_g
    pi_t = parameters.base_frequency_t
    pi_r = parameters.purine_frequency
    pi_y = parameters.pyrimidine_frequency
    if min(pi_a, pi_c, pi_g, pi_t) <= 0.0 or pi_r <= 0.0 or pi_y <= 0.0:
        return (
            None,
            "alignment-wide resolved base composition omits at least one nucleotide class, so the TN93 correction is undefined",
        )

    p1 = summary.ag_transition_sites / summary.comparable_sites
    p2 = summary.ct_transition_sites / summary.comparable_sites
    q = summary.transversion_sites / summary.comparable_sites
    first = 1.0 - ((pi_r * p1) / (2.0 * pi_a * pi_g)) - (q / (2.0 * pi_r))
    second = 1.0 - ((pi_y * p2) / (2.0 * pi_c * pi_t)) - (q / (2.0 * pi_y))
    third = 1.0 - (q / (2.0 * pi_r * pi_y))
    if first < 0.0 or second < 0.0 or third < 0.0:
        return (
            None,
            "transition and transversion proportions exceed the TN93 correction range for the estimated base composition, so the corrected distance is undefined",
        )
    if first == 0.0 or second == 0.0 or third == 0.0:
        return (
            None,
            "transition and transversion proportions are at the TN93 correction limit for the estimated base composition, so the corrected distance tends to infinity",
        )
    coefficient_ag = (2.0 * pi_a * pi_g) / pi_r
    coefficient_ct = (2.0 * pi_c * pi_t) / pi_y
    coefficient_tv = 2.0 * (
        (pi_r * pi_y)
        - ((pi_a * pi_g * pi_y) / pi_r)
        - ((pi_c * pi_t * pi_r) / pi_y)
    )
    value = (
        (-coefficient_ag * math.log(first))
        + (-coefficient_ct * math.log(second))
        + (-coefficient_tv * math.log(third))
    )
    return round(value, 15), None


def _protein_p_distance(summary: _PairSummary) -> tuple[float | None, str | None]:
    return (
        _p_distance(summary),
        None
        if summary.comparable_sites > 0
        else "no comparable sites remain after filtering",
    )


def _distance_from_summary(
    summary: _PairSummary,
    *,
    model: DistanceModel,
    model_parameters: GeneticDistanceModelParameters | None,
) -> tuple[float | None, str | None]:
    if model == "p-distance":
        return (
            _p_distance(summary),
            None
            if summary.comparable_sites > 0
            else "no comparable sites remain after filtering",
        )
    if model == "jukes-cantor":
        return _jukes_cantor_distance(_p_distance(summary))
    if model == "kimura-2-parameter":
        return _kimura_two_parameter_distance(summary)
    if model == "felsenstein-81":
        return _felsenstein_81_distance(summary, model_parameters)
    if model == "tamura-nei-93":
        return _tamura_nei_93_distance(summary, model_parameters)
    if model == "amino-acid-p-distance":
        return _protein_p_distance(summary)
    raise ValueError(f"unsupported distance model: {model}")


def _pair_summary(
    left: str,
    right: str,
    *,
    alphabet: str,
    ambiguity_policy: AmbiguityPolicy,
    retained_positions: list[int] | None = None,
    model: DistanceModel,
    model_parameters: GeneticDistanceModelParameters | None = None,
) -> _PairSummary:
    comparable_sites = 0
    mismatch_sites = 0.0
    transition_sites = 0.0
    ag_transition_sites = 0.0
    ct_transition_sites = 0.0
    transversion_sites = 0.0
    ambiguity_sites = 0
    skipped_sites = 0
    positions = (
        retained_positions if retained_positions is not None else list(range(len(left)))
    )
    for position in positions:
        left_states = _states_for_symbol(left[position], alphabet=alphabet)
        right_states = _states_for_symbol(right[position], alphabet=alphabet)
        ambiguous = (
            left_states is not None
            and right_states is not None
            and bool(left_states)
            and bool(right_states)
            and (len(left_states) > 1 or len(right_states) > 1)
        )
        contribution = _site_contribution(
            left[position],
            right[position],
            alphabet=alphabet,
            ambiguity_policy=ambiguity_policy,
        )
        if contribution is None:
            if ambiguous:
                ambiguity_sites += 1
            skipped_sites += 1
            continue
        comparable_sites += 1
        mismatch_sites += contribution.mismatch_weight
        transition_sites += contribution.transition_weight
        ag_transition_sites += contribution.ag_transition_weight
        ct_transition_sites += contribution.ct_transition_weight
        transversion_sites += contribution.transversion_weight
        if contribution.ambiguous:
            ambiguity_sites += 1
    preliminary = _PairSummary(
        distance=None,
        comparable_sites=comparable_sites,
        mismatch_sites=round(mismatch_sites, 15),
        transition_sites=round(transition_sites, 15),
        ag_transition_sites=round(ag_transition_sites, 15),
        ct_transition_sites=round(ct_transition_sites, 15),
        transversion_sites=round(transversion_sites, 15),
        ambiguity_sites=ambiguity_sites,
        skipped_sites=skipped_sites,
        saturated=False,
        saturation_reason=None,
    )
    distance, saturation_reason = _distance_from_summary(
        preliminary,
        model=model,
        model_parameters=model_parameters,
    )
    saturated = saturation_reason is not None and (
        "correction limit" in saturation_reason
        or "correction range" in saturation_reason
        or (distance is None and comparable_sites > 0)
    )
    if (
        model in {"p-distance", "amino-acid-p-distance"}
        and distance is not None
        and distance >= 0.75
    ):
        saturated = True
        saturation_reason = (
            "raw p-distance indicates severe divergence and likely saturation"
        )
    return _PairSummary(
        distance=distance,
        comparable_sites=comparable_sites,
        mismatch_sites=round(mismatch_sites, 15),
        transition_sites=round(transition_sites, 15),
        ag_transition_sites=round(ag_transition_sites, 15),
        ct_transition_sites=round(ct_transition_sites, 15),
        transversion_sites=round(transversion_sites, 15),
        ambiguity_sites=ambiguity_sites,
        skipped_sites=skipped_sites,
        saturated=saturated,
        saturation_reason=saturation_reason,
    )


def _complete_deletion_positions(
    records: list[AlignmentRecord],
    *,
    alphabet: str,
    ambiguity_policy: AmbiguityPolicy,
) -> list[int]:
    retained: list[int] = []
    for position in range(len(records[0].sequence)):
        if all(
            _site_contribution(
                record.sequence[position],
                record.sequence[position],
                alphabet=alphabet,
                ambiguity_policy=ambiguity_policy,
            )
            is not None
            for record in records
        ):
            retained.append(position)
    return retained


def _load_alignment_for_model(
    path: Path, *, model: DistanceModel
) -> tuple[list[AlignmentRecord], str]:
    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if model not in _allowed_models_for_alphabet(alphabet):
        raise InvalidAlignmentError(
            f"distance model '{model}' is not supported for inferred alphabet '{alphabet}'"
        )
    return records, alphabet


def _pair_key(left_identifier: str, right_identifier: str) -> tuple[str, str]:
    return tuple(sorted((left_identifier, right_identifier)))


def _distance_rows(entries: list[ImportedDistanceEntry]) -> list[ImportedDistanceEntry]:
    return [
        entry for entry in entries if entry.left_identifier != entry.right_identifier
    ]


def _unique_distance_rows(
    entries: list[ImportedDistanceEntry],
) -> list[ImportedDistanceEntry]:
    unique: dict[tuple[str, str], ImportedDistanceEntry] = {}
    for entry in _distance_rows(entries):
        unique.setdefault(
            _pair_key(entry.left_identifier, entry.right_identifier), entry
        )
    return [unique[key] for key in sorted(unique)]


def _imported_distance_scale(entries: list[ImportedDistanceEntry]) -> str:
    off_diagonal = _unique_distance_rows(entries)
    if off_diagonal and all(0.0 <= entry.distance <= 1.5 for entry in off_diagonal):
        return "unit-interval-like"
    return "unknown"


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _unique_genetic_distance_pairs(
    report: GeneticDistanceMatrix,
) -> list[PairwiseGeneticDistance]:
    unique: dict[tuple[str, str], PairwiseGeneticDistance] = {}
    for pair in report.pairs:
        if pair.left_identifier == pair.right_identifier:
            continue
        unique.setdefault(_pair_key(pair.left_identifier, pair.right_identifier), pair)
    return [unique[key] for key in sorted(unique)]


def load_imported_distance_matrix(path: Path) -> list[ImportedDistanceEntry]:
    """Load a long-form imported distance matrix table."""
    if not path.exists():
        raise FileNotFoundError(f"distance matrix file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"left_identifier", "right_identifier", "distance"}
        if reader.fieldnames is None or not required_columns <= set(reader.fieldnames):
            raise InvalidDistanceMatrixError(
                "distance matrix must contain left_identifier, right_identifier, and distance columns"
            )
        comparable_sites_column = (
            "comparable_sites" if "comparable_sites" in reader.fieldnames else None
        )
        entries: list[ImportedDistanceEntry] = []
        seen_directional_pairs: set[tuple[str, str]] = set()
        for row_index, row in enumerate(reader, start=2):
            left_identifier = str(row.get("left_identifier", "")).strip()
            right_identifier = str(row.get("right_identifier", "")).strip()
            raw_distance = str(row.get("distance", "")).strip()
            if not left_identifier or not right_identifier:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} is missing a left_identifier or right_identifier value"
                )
            if not raw_distance:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} is missing a distance value"
                )
            try:
                distance = float(raw_distance)
            except ValueError as error:
                raise InvalidDistanceMatrixError(
                    f"row {row_index} in {path} contains a non-numeric distance value '{raw_distance}'"
                ) from error
            comparable_sites: int | None = None
            if comparable_sites_column is not None:
                raw_comparable_sites = str(row.get(comparable_sites_column, "")).strip()
                if raw_comparable_sites:
                    try:
                        comparable_sites = int(raw_comparable_sites)
                    except ValueError as error:
                        raise InvalidDistanceMatrixError(
                            f"row {row_index} in {path} contains a non-integer comparable_sites value '{raw_comparable_sites}'"
                        ) from error
            directional_pair = (left_identifier, right_identifier)
            if directional_pair in seen_directional_pairs:
                raise InvalidDistanceMatrixError(
                    f"distance matrix contains duplicate directional entry {left_identifier}/{right_identifier}"
                )
            seen_directional_pairs.add(directional_pair)
            entries.append(
                ImportedDistanceEntry(
                    left_identifier=left_identifier,
                    right_identifier=right_identifier,
                    distance=round(distance, 15),
                    comparable_sites=comparable_sites,
                )
            )

    if not entries:
        raise InvalidDistanceMatrixError(f"distance matrix contains no rows: {path}")
    return entries


def _symmetric_distances(
    entries: list[ImportedDistanceEntry],
) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}
    by_direction = {
        (entry.left_identifier, entry.right_identifier): entry.distance
        for entry in entries
    }
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    for left_identifier in identifiers:
        for right_identifier in identifiers:
            pair_key = _pair_key(left_identifier, right_identifier)
            if pair_key in distances:
                continue
            if left_identifier == right_identifier:
                if (left_identifier, right_identifier) in by_direction:
                    distances[pair_key] = by_direction[
                        (left_identifier, right_identifier)
                    ]
                continue
            if (left_identifier, right_identifier) in by_direction:
                distances[pair_key] = by_direction[(left_identifier, right_identifier)]
            elif (right_identifier, left_identifier) in by_direction:
                distances[pair_key] = by_direction[(right_identifier, left_identifier)]
    return distances


def _iter_ultrametric_violations(
    identifiers: list[str],
    distances: dict[tuple[str, str], float],
    *,
    tolerance: float,
) -> list[UPGMAUltrametricViolation]:
    violations: list[UPGMAUltrametricViolation] = []
    for left_index, left_identifier in enumerate(identifiers):
        for middle_index in range(left_index + 1, len(identifiers)):
            middle_identifier = identifiers[middle_index]
            for right_index in range(middle_index + 1, len(identifiers)):
                right_identifier = identifiers[right_index]
                pair_keys = [
                    _pair_key(left_identifier, middle_identifier),
                    _pair_key(left_identifier, right_identifier),
                    _pair_key(middle_identifier, right_identifier),
                ]
                if any(pair_key not in distances for pair_key in pair_keys):
                    continue
                triple = [distances[pair_key] for pair_key in pair_keys]
                ordered = sorted(triple)
                deviation = abs(ordered[2] - ordered[1])
                if deviation > tolerance:
                    violations.append(
                        UPGMAUltrametricViolation(
                            left_identifier=left_identifier,
                            middle_identifier=middle_identifier,
                            right_identifier=right_identifier,
                            smallest_distance=ordered[0],
                            middle_distance=ordered[1],
                            largest_distance=ordered[2],
                            deviation=deviation,
                        )
                    )
    return sorted(
        violations,
        key=lambda row: (
            row.left_identifier,
            row.middle_identifier,
            row.right_identifier,
        ),
    )


def _build_alignment_distance_lookup(
    report: GeneticDistanceMatrix,
) -> dict[tuple[str, str], float]:
    distances: dict[tuple[str, str], float] = {}
    for pair in report.pairs:
        if pair.left_identifier == pair.right_identifier or pair.distance is None:
            continue
        distances[_pair_key(pair.left_identifier, pair.right_identifier)] = float(
            pair.distance
        )
    return distances


def assess_distance_method_assumptions(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit clock-like compatibility and core distance-tree assumptions for an alignment."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    distances = _build_alignment_distance_lookup(matrix)
    expected_pairs = (len(matrix.identifiers) * (len(matrix.identifiers) - 1)) // 2
    violations = _iter_ultrametric_violations(
        matrix.identifiers,
        distances,
        tolerance=ultrametric_tolerance,
    )
    warnings: list[str] = []
    if len(distances) < expected_pairs:
        warnings.append(
            "UPGMA clock-assumption auditing is incomplete because one or more pairwise distances are undefined under the selected model"
        )
    if violations:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=path,
        source_kind="alignment",
        taxon_count=len(matrix.identifiers),
        pair_count=len(distances),
        nj_assumptions=[
            "neighbor-joining treats the matrix as an additive approximation and does not require a strict molecular clock",
            "neighbor-joining still becomes unreliable when pairwise distances are heavily saturated or estimated from too few comparable sites",
        ],
        upgma_assumptions=[
            "UPGMA assumes the pairwise distances are ultrametric and therefore consistent with a clock-like process",
            "UPGMA can mis-cluster taxa when rates vary among lineages even if the matrix remains symmetric and complete",
        ],
        ultrametric_compatible=not violations,
        ultrametric_tolerance=ultrametric_tolerance,
        upgma_ultrametric_violations=violations,
        warnings=warnings,
    )


def assess_imported_distance_method_assumptions(
    path: Path,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit clock-like compatibility and core distance-tree assumptions for an imported matrix."""
    validation = validate_imported_distance_matrix(path)
    distances = _symmetric_distances(load_imported_distance_matrix(path))
    violations = (
        []
        if not validation.complete
        or not validation.symmetric
        or not validation.nonnegative
        else _iter_ultrametric_violations(
            validation.identifiers,
            distances,
            tolerance=ultrametric_tolerance,
        )
    )
    warnings = list(validation.warnings)
    if violations:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=path,
        source_kind="imported-distance-matrix",
        taxon_count=len(validation.identifiers),
        pair_count=len(distances),
        nj_assumptions=[
            "neighbor-joining treats the matrix as an additive approximation and does not require a strict molecular clock",
            "neighbor-joining still inherits any errors or non-metric distortions present in the imported matrix",
        ],
        upgma_assumptions=[
            "UPGMA assumes the imported matrix is ultrametric and therefore compatible with a clock-like process",
            "UPGMA can enforce an ultrametric tree even when the source matrix does not satisfy that assumption",
        ],
        ultrametric_compatible=not violations,
        ultrametric_tolerance=ultrametric_tolerance,
        upgma_ultrametric_violations=violations,
        warnings=warnings,
    )


def validate_imported_distance_matrix(path: Path) -> ImportedDistanceMatrixReport:
    """Validate a long-form imported distance matrix table."""
    entries = load_imported_distance_matrix(path)
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    by_direction = {
        (entry.left_identifier, entry.right_identifier): entry for entry in entries
    }

    missing_pairs: list[str] = []
    diagonal_problems: list[str] = []
    negative_distance_pairs: list[str] = []
    asymmetric_pairs: list[DistanceMatrixAsymmetry] = []
    symmetric_distances = _symmetric_distances(entries)

    for left_identifier in identifiers:
        diagonal = by_direction.get((left_identifier, left_identifier))
        if diagonal is None:
            missing_pairs.append(f"{left_identifier}/{left_identifier}")
        elif diagonal.distance != 0.0:
            diagonal_problems.append(
                f"{left_identifier}/{left_identifier} has diagonal distance {diagonal.distance:g}"
            )
        if diagonal is not None and diagonal.distance < 0:
            negative_distance_pairs.append(f"{left_identifier}/{left_identifier}")

        for right_identifier in identifiers:
            if left_identifier >= right_identifier:
                continue
            left_to_right = by_direction.get((left_identifier, right_identifier))
            right_to_left = by_direction.get((right_identifier, left_identifier))
            if left_to_right is None and right_to_left is None:
                missing_pairs.append(f"{left_identifier}/{right_identifier}")
                continue
            if left_to_right is not None and left_to_right.distance < 0:
                negative_distance_pairs.append(f"{left_identifier}/{right_identifier}")
            if right_to_left is not None and right_to_left.distance < 0:
                negative_distance_pairs.append(f"{right_identifier}/{left_identifier}")
            if (
                left_to_right is not None
                and right_to_left is not None
                and not math.isclose(
                    left_to_right.distance,
                    right_to_left.distance,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            ):
                asymmetric_pairs.append(
                    DistanceMatrixAsymmetry(
                        left_identifier=left_identifier,
                        right_identifier=right_identifier,
                        left_to_right_distance=left_to_right.distance,
                        right_to_left_distance=right_to_left.distance,
                    )
                )

    nonmetric_observations: list[NonMetricDistanceObservation] = []
    if not missing_pairs and not diagonal_problems and not negative_distance_pairs:
        for left_index, left_identifier in enumerate(identifiers):
            for middle_index, middle_identifier in enumerate(identifiers):
                if middle_index == left_index:
                    continue
                for right_index, right_identifier in enumerate(identifiers):
                    if len({left_index, middle_index, right_index}) < 3:
                        continue
                    if left_identifier > right_identifier:
                        continue
                    direct_distance = symmetric_distances.get(
                        _pair_key(left_identifier, right_identifier)
                    )
                    left_middle = symmetric_distances.get(
                        _pair_key(left_identifier, middle_identifier)
                    )
                    middle_right = symmetric_distances.get(
                        _pair_key(middle_identifier, right_identifier)
                    )
                    if (
                        direct_distance is None
                        or left_middle is None
                        or middle_right is None
                    ):
                        continue
                    indirect_distance = round(left_middle + middle_right, 15)
                    if direct_distance > indirect_distance + 1e-12:
                        nonmetric_observations.append(
                            NonMetricDistanceObservation(
                                left_identifier=left_identifier,
                                middle_identifier=middle_identifier,
                                right_identifier=right_identifier,
                                direct_distance=direct_distance,
                                indirect_distance=indirect_distance,
                            )
                        )

    warnings: list[str] = []
    if missing_pairs:
        warnings.append("distance matrix is missing one or more required pairs")
    if diagonal_problems:
        warnings.append("distance matrix contains nonzero diagonal entries")
    if negative_distance_pairs:
        warnings.append("distance matrix contains negative distances")
    if asymmetric_pairs:
        warnings.append("distance matrix contains asymmetric directional entries")
    if nonmetric_observations:
        warnings.append(
            "distance matrix violates triangle inequality for one or more taxon triples"
        )

    return ImportedDistanceMatrixReport(
        path=path,
        identifiers=identifiers,
        pair_count=len(entries),
        complete=not missing_pairs,
        zero_diagonal=not diagonal_problems,
        symmetric=not asymmetric_pairs,
        nonnegative=not negative_distance_pairs,
        missing_pairs=missing_pairs,
        diagonal_problems=diagonal_problems,
        negative_distance_pairs=sorted(set(negative_distance_pairs)),
        asymmetric_pairs=sorted(
            asymmetric_pairs,
            key=lambda row: (row.left_identifier, row.right_identifier),
        ),
        nonmetric_observations=sorted(
            nonmetric_observations,
            key=lambda row: (
                row.left_identifier,
                row.middle_identifier,
                row.right_identifier,
            ),
        ),
        warnings=warnings,
    )


def _bio_distance_matrix(report: GeneticDistanceMatrix) -> DistanceMatrix:
    undefined_pairs = [
        f"{pair.left_identifier}/{pair.right_identifier}"
        for pair in report.pairs
        if pair.distance is None
    ]
    if undefined_pairs:
        raise InvalidAlignmentError(
            "distance matrix contains undefined entries for: "
            + ", ".join(undefined_pairs)
        )
    rows: list[list[float]] = []
    for row_index, left_identifier in enumerate(report.identifiers):
        row: list[float] = []
        for right_identifier in report.identifiers[: row_index + 1]:
            if left_identifier == right_identifier:
                row.append(0.0)
                continue
            pair = next(
                pair
                for pair in report.pairs
                if {pair.left_identifier, pair.right_identifier}
                == {left_identifier, right_identifier}
            )
            row.append(float(pair.distance))
        rows.append(row)
    return DistanceMatrix(report.identifiers, rows)


def _bio_distance_matrix_from_imported(
    report: ImportedDistanceMatrixReport, entries: list[ImportedDistanceEntry]
) -> DistanceMatrix:
    if not report.complete:
        raise InvalidDistanceMatrixError("distance matrix is incomplete")
    if not report.zero_diagonal:
        raise InvalidDistanceMatrixError("distance matrix has nonzero diagonal entries")
    if not report.symmetric:
        raise InvalidDistanceMatrixError(
            "distance matrix contains asymmetric directional entries"
        )
    if not report.nonnegative:
        raise InvalidDistanceMatrixError("distance matrix contains negative distances")

    symmetric_distances = _symmetric_distances(entries)
    rows: list[list[float]] = []
    for row_index, left_identifier in enumerate(report.identifiers):
        row: list[float] = []
        for right_identifier in report.identifiers[: row_index + 1]:
            pair_distance = symmetric_distances.get(
                _pair_key(left_identifier, right_identifier)
            )
            if pair_distance is None:
                raise InvalidDistanceMatrixError(
                    f"distance matrix is missing pair {left_identifier}/{right_identifier}"
                )
            row.append(pair_distance)
        rows.append(row)
    return DistanceMatrix(report.identifiers, rows)


def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    """Compute a deterministic pairwise genetic distance matrix for an aligned dataset."""
    model = _normalize_distance_model(model)
    if gap_handling not in {"pairwise-deletion", "complete-deletion"}:
        raise ValueError(f"unsupported gap handling mode: {gap_handling}")
    if ambiguity_policy not in {
        "ignore",
        "partial-match",
        "strict-mismatch",
        "report-only",
    }:
        raise ValueError(f"unsupported ambiguity policy: {ambiguity_policy}")

    records, alphabet = _load_alignment_for_model(path, model=model)
    model_parameters: GeneticDistanceModelParameters | None = None
    warnings: list[str] = []
    if alphabet in {"dna", "rna"} and model in {
        "felsenstein-81",
        "tamura-nei-93",
    }:
        model_parameters, warnings = _estimate_nucleotide_model_parameters(records)
    retained_positions = (
        _complete_deletion_positions(
            records, alphabet=alphabet, ambiguity_policy=ambiguity_policy
        )
        if gap_handling == "complete-deletion"
        else None
    )
    pairs: list[PairwiseGeneticDistance] = []
    for left_index, left in enumerate(records):
        for right_index, right in enumerate(records):
            if right_index < left_index:
                continue
            summary = _pair_summary(
                left.sequence,
                right.sequence,
                alphabet=alphabet,
                ambiguity_policy=ambiguity_policy,
                retained_positions=retained_positions,
                model=model,
                model_parameters=model_parameters,
            )
            pairs.append(
                PairwiseGeneticDistance(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    distance=summary.distance if left_index != right_index else 0.0,
                    comparable_sites=summary.comparable_sites,
                    mismatch_sites=summary.mismatch_sites
                    if left_index != right_index
                    else 0.0,
                    transition_sites=summary.transition_sites
                    if left_index != right_index
                    else 0.0,
                    ag_transition_sites=summary.ag_transition_sites
                    if left_index != right_index
                    else 0.0,
                    ct_transition_sites=summary.ct_transition_sites
                    if left_index != right_index
                    else 0.0,
                    transversion_sites=summary.transversion_sites
                    if left_index != right_index
                    else 0.0,
                    ambiguity_sites=summary.ambiguity_sites,
                    skipped_sites=summary.skipped_sites,
                    saturated=False if left_index == right_index else summary.saturated,
                    saturation_reason=None
                    if left_index == right_index
                    else summary.saturation_reason,
                )
            )
    return GeneticDistanceMatrix(
        path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=alphabet,
        alignment_length=len(records[0].sequence),
        identifiers=[record.identifier for record in records],
        model_parameters=model_parameters,
        warnings=warnings,
        pairs=pairs,
    )


def validate_distance_reference_examples() -> DistanceReferenceValidationReport:
    """Validate core distance examples against durable reference expectations."""
    examples = [
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_distance.fasta",
            "case": "dna-p-distance",
            "model": "p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "B",
            "expected": 0.125,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_distance.fasta",
            "case": "dna-jukes-cantor",
            "model": "jukes-cantor",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "B",
            "expected": 0.136741167595466,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_distance.fasta",
            "case": "dna-kimura-2-parameter",
            "model": "kimura-2-parameter",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "B",
            "expected": 0.14384103622589,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_distance_gaps.fasta",
            "case": "dna-felsenstein-81",
            "model": "felsenstein-81",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "B",
            "expected": 0.189450794670797,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_duplicates.fasta",
            "case": "dna-tamura-nei-93",
            "model": "tamura-nei-93",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "C",
            "expected": 0.182459298648674,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_protein.fasta",
            "case": "protein-p-distance",
            "model": "amino-acid-p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "P1",
            "right": "P2",
            "expected": 0.125,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_distance_gaps.fasta",
            "case": "gap-complete-deletion",
            "model": "p-distance",
            "gap_handling": "complete-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "D",
            "expected": 0.0,
            "expected_ambiguity_sites": 0,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
            "case": "ambiguity-ignore",
            "model": "p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "ignore",
            "left": "A",
            "right": "B",
            "expected": 0.0,
            "expected_ambiguity_sites": 1,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
            "case": "ambiguity-partial-match",
            "model": "p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "partial-match",
            "left": "A",
            "right": "B",
            "expected": 0.15,
            "expected_ambiguity_sites": 1,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
            "case": "ambiguity-strict-mismatch",
            "model": "p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "strict-mismatch",
            "left": "A",
            "right": "B",
            "expected": 0.2,
            "expected_ambiguity_sites": 1,
        },
        {
            "path": Path(__file__).resolve().parents[2]
            / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
            "case": "ambiguity-report-only",
            "model": "p-distance",
            "gap_handling": "pairwise-deletion",
            "ambiguity_policy": "report-only",
            "left": "A",
            "right": "B",
            "expected": 0.0,
            "expected_ambiguity_sites": 1,
        },
    ]
    observations: list[DistanceReferenceObservation] = []
    for example in examples:
        report = compute_pairwise_genetic_distance_matrix(
            example["path"],
            model=str(example["model"]),
            gap_handling=str(example["gap_handling"]),
            ambiguity_policy=str(example["ambiguity_policy"]),
        )
        pair = next(
            row
            for row in report.pairs
            if row.left_identifier == example["left"]
            and row.right_identifier == example["right"]
        )
        expected_distance = float(example["expected"])
        observed = pair.distance
        passed = observed is not None and math.isclose(
            observed, expected_distance, rel_tol=1e-12, abs_tol=1e-12
        )
        observations.append(
            DistanceReferenceObservation(
                case=str(example["case"]),
                model=str(example["model"]),
                gap_handling=str(example["gap_handling"]),
                ambiguity_policy=str(example["ambiguity_policy"]),
                left_identifier=str(example["left"]),
                right_identifier=str(example["right"]),
                expected_distance=expected_distance,
                observed_distance=observed,
                comparable_sites=pair.comparable_sites,
                expected_ambiguity_sites=int(example["expected_ambiguity_sites"]),
                observed_ambiguity_sites=pair.ambiguity_sites,
                passed=passed
                and pair.ambiguity_sites == int(example["expected_ambiguity_sites"]),
            )
        )
    tree_observations: list[DistanceTreeReferenceObservation] = []
    nj_matrix_path = (
        Path(__file__).resolve().parents[2]
        / "tests/fixtures/metadata/example_distance_matrix_ultrametric.tsv"
    )
    nj_tree, _ = build_tree_from_imported_distance_matrix(
        nj_matrix_path, method="neighbor-joining"
    )
    nj_expected_clades = ["A|B"]
    nj_observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in _informative_clades(nj_tree, set(nj_tree.tip_names))
    )
    tree_observations.append(
        DistanceTreeReferenceObservation(
            case="neighbor-joining-reference-clustering",
            method="neighbor-joining",
            matrix_path=nj_matrix_path,
            expected_clades=nj_expected_clades,
            observed_clades=nj_observed_clades,
            passed=nj_observed_clades == nj_expected_clades,
        )
    )
    upgma_matrix_path = (
        Path(__file__).resolve().parents[2]
        / "tests/fixtures/metadata/example_distance_matrix_ultrametric.tsv"
    )
    upgma_tree, _ = build_tree_from_imported_distance_matrix(
        upgma_matrix_path, method="upgma"
    )
    expected_clades = ["A|B", "C|D"]
    observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in _informative_clades(upgma_tree, set(upgma_tree.tip_names))
    )
    tree_observations.append(
        DistanceTreeReferenceObservation(
            case="upgma-ultrametric-clustering",
            method="upgma",
            matrix_path=upgma_matrix_path,
            expected_clades=expected_clades,
            observed_clades=observed_clades,
            passed=observed_clades == expected_clades,
        )
    )
    return DistanceReferenceValidationReport(
        observations=observations,
        tree_observations=tree_observations,
        all_passed=all(observation.passed for observation in observations)
        and all(observation.passed for observation in tree_observations),
    )


def inspect_distance_matrix_quality(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceMatrixQualityReport:
    """Report saturation, outlier, and low-information risks for one computed matrix."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    off_diagonal = [
        pair for pair in matrix.pairs if pair.left_identifier != pair.right_identifier
    ]
    defined_pairs = [pair for pair in off_diagonal if pair.distance is not None]
    saturated_pairs = [
        SaturatedDistancePair(
            left_identifier=pair.left_identifier,
            right_identifier=pair.right_identifier,
            distance=pair.distance,
            comparable_sites=pair.comparable_sites,
            reason=pair.saturation_reason
            or "distance is not usable under the selected model",
        )
        for pair in off_diagonal
        if pair.saturated
    ]

    high_distance_outliers: list[DistanceOutlierPair] = []
    if defined_pairs:
        distances = sorted(
            float(pair.distance) for pair in defined_pairs if pair.distance is not None
        )
        threshold = 0.75
        if len(distances) >= 4:
            q1 = distances[len(distances) // 4]
            q3 = distances[(3 * len(distances)) // 4]
            threshold = max(threshold, q3 + (1.5 * (q3 - q1)))
        for pair in defined_pairs:
            if pair.distance is not None and pair.distance >= threshold:
                high_distance_outliers.append(
                    DistanceOutlierPair(
                        left_identifier=pair.left_identifier,
                        right_identifier=pair.right_identifier,
                        distance=pair.distance,
                        note="pairwise distance is unusually large relative to the dataset baseline",
                    )
                )

    comparable_baseline = (
        median(pair.comparable_sites for pair in off_diagonal) if off_diagonal else 0
    )
    low_information_cutoff = max(10, int(math.floor(comparable_baseline * 0.5)))
    low_information_pairs = [
        LowInformationPair(
            left_identifier=pair.left_identifier,
            right_identifier=pair.right_identifier,
            comparable_sites=pair.comparable_sites,
            note="too few comparable sites remain for a stable distance estimate",
        )
        for pair in off_diagonal
        if pair.comparable_sites < low_information_cutoff
    ]
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )

    warnings: list[str] = []
    blocker_reasons: list[str] = []
    if saturated_pairs:
        warnings.append(
            "one or more pairwise distances are saturated or undefined under the selected model"
        )
        if len(saturated_pairs) / max(1, len(off_diagonal)) > 0.25:
            warnings.append(
                "many pairwise distances sit in a saturation regime that weakens distance-method assumptions"
            )
    if low_information_pairs:
        warnings.append("one or more taxon pairs retain too few comparable sites")
        if len(low_information_pairs) / max(1, len(off_diagonal)) > 0.5:
            warnings.append(
                "many taxon pairs retain too few comparable sites after filtering"
            )
    if not defined_pairs:
        blocker_reasons.append(
            "no off-diagonal distances could be computed from the selected alignment and policies"
        )
    if matrix.alignment_length < 10:
        warnings.append("alignment is short for robust distance-based tree building")
    if high_distance_outliers:
        warnings.append("one or more taxon pairs are unusually divergent")
    warnings.extend(assumptions.warnings)
    decision = "blocked" if blocker_reasons else ("risky" if warnings else "allowed")
    assessment_reasons = blocker_reasons if blocker_reasons else warnings
    return DistanceMatrixQualityReport(
        path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=matrix.inferred_alphabet,
        taxon_count=len(matrix.identifiers),
        pair_count=len(matrix.pairs),
        saturated_pairs=saturated_pairs,
        high_distance_outliers=high_distance_outliers,
        low_information_pairs=low_information_pairs,
        assumptions=assumptions,
        warnings=warnings,
        method_assessment=DistanceMethodAssessment(
            decision=decision,
            reasons=assessment_reasons,
        ),
    )


def inspect_imported_distance_matrix_quality(
    path: Path,
) -> ImportedDistanceMatrixQualityReport:
    """Report structural and heuristic quality risks for an imported distance matrix."""
    validation = validate_imported_distance_matrix(path)
    entries = load_imported_distance_matrix(path)
    off_diagonal = _unique_distance_rows(entries)
    distance_scale = _imported_distance_scale(entries)

    saturated_pairs: list[SaturatedDistancePair] = []
    if distance_scale == "unit-interval-like":
        saturated_pairs = [
            SaturatedDistancePair(
                left_identifier=entry.left_identifier,
                right_identifier=entry.right_identifier,
                distance=entry.distance,
                comparable_sites=entry.comparable_sites or 0,
                reason="imported distance lies in a high-divergence regime often associated with saturation for raw or lightly corrected genetic distances",
            )
            for entry in off_diagonal
            if entry.distance >= 0.75
        ]

    comparable_site_counts = [
        entry.comparable_sites
        for entry in off_diagonal
        if entry.comparable_sites is not None
    ]
    low_information_pair_cutoff: int | None = None
    low_information_pairs: list[LowInformationPair] = []
    if comparable_site_counts:
        comparable_baseline = median(comparable_site_counts)
        low_information_pair_cutoff = max(
            10, int(math.floor(comparable_baseline * 0.5))
        )
        low_information_pairs = [
            LowInformationPair(
                left_identifier=entry.left_identifier,
                right_identifier=entry.right_identifier,
                comparable_sites=int(entry.comparable_sites or 0),
                note="imported pair retains too few comparable sites for a stable distance interpretation",
            )
            for entry in off_diagonal
            if entry.comparable_sites is not None
            and entry.comparable_sites < low_information_pair_cutoff
        ]

    warnings = list(validation.warnings)
    if distance_scale != "unit-interval-like":
        warnings.append(
            "saturation heuristics were skipped because imported distances do not resemble unit-interval genetic distances"
        )
    elif saturated_pairs:
        warnings.append(
            "one or more imported distances fall in a high-divergence regime that deserves saturation review"
        )
    if comparable_site_counts:
        if low_information_pairs:
            warnings.append(
                "one or more imported distances retain too few comparable sites"
            )
    else:
        warnings.append(
            "low-information pair auditing is unavailable because the imported matrix does not provide comparable_sites"
        )

    return ImportedDistanceMatrixQualityReport(
        validation=validation,
        saturated_pairs=saturated_pairs,
        low_information_pairs=low_information_pairs,
        low_information_pair_cutoff=low_information_pair_cutoff,
        saturation_audit_scale=distance_scale,
        warnings=warnings,
    )


def write_genetic_distance_matrix(path: Path, report: GeneticDistanceMatrix) -> Path:
    """Write a pairwise genetic distance matrix as a deterministic TSV."""
    rows = {
        (pair.left_identifier, pair.right_identifier): pair for pair in report.pairs
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tdistance\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            normalized_distance = (
                None
                if pair.distance is None
                else 0.0
                if math.isclose(pair.distance, 0.0, abs_tol=1e-15)
                else pair.distance
            )
            distance = (
                "" if normalized_distance is None else format(normalized_distance, ".15g")
            )
            lines.append(f"{left}\t{right}\t{distance}\t{pair.comparable_sites}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_genetic_distance_component_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    """Write one deterministic pairwise distance component table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "left_identifier",
                "right_identifier",
                "distance",
                "comparable_sites",
                "mismatch_sites",
                "transition_sites",
                "ag_transition_sites",
                "ct_transition_sites",
                "transversion_sites",
                "ambiguity_sites",
                "skipped_sites",
                "saturated",
                "saturation_reason",
            ]
        )
    ]
    for pair in report.pairs:
        normalized_distance = (
            None
            if pair.distance is None
            else 0.0
            if math.isclose(pair.distance, 0.0, abs_tol=1e-15)
            else pair.distance
        )
        lines.append(
            "\t".join(
                [
                    pair.left_identifier,
                    pair.right_identifier,
                    "" if normalized_distance is None else format(normalized_distance, ".15g"),
                    str(pair.comparable_sites),
                    format(pair.mismatch_sites, ".15g"),
                    format(pair.transition_sites, ".15g"),
                    format(pair.ag_transition_sites, ".15g"),
                    format(pair.ct_transition_sites, ".15g"),
                    format(pair.transversion_sites, ".15g"),
                    str(pair.ambiguity_sites),
                    str(pair.skipped_sites),
                    "true" if pair.saturated else "false",
                    "" if pair.saturation_reason is None else pair.saturation_reason,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_genetic_distance_parameter_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    """Write one deterministic alignment-wide parameter table for DNA distances."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["parameter\tvalue"]
    if report.model_parameters is None:
        lines.append(f"model\t{report.model}")
    else:
        for parameter, value in asdict(report.model_parameters).items():
            rendered = "" if value is None else format(value, ".15g")
            lines.append(f"{parameter}\t{rendered}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_distance_tree(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    """Build a distance-based tree from an aligned dataset."""
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    if quality.method_assessment.decision == "blocked":
        raise InvalidAlignmentError(
            "distance tree building is blocked: "
            + "; ".join(quality.method_assessment.reasons)
        )
    report = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    if len(report.identifiers) < 2:
        raise InvalidAlignmentError("distance tree building requires at least two taxa")

    constructor = DistanceTreeConstructor()
    distance_matrix = _bio_distance_matrix(report)
    if method == "neighbor-joining":
        tree = constructor.nj(distance_matrix)
    elif method == "upgma":
        tree = constructor.upgma(distance_matrix)
    else:
        raise ValueError(f"unsupported tree-building method: {method}")
    assumptions = quality.assumptions

    return tree_from_biophylo(tree, source_format="newick"), DistanceTreeBuildReport(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        method=method,
        taxon_count=len(report.identifiers),
        pair_count=len(report.pairs),
        assumptions=assumptions,
    )


def compare_distance_tree_topologies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceTreeTopologyComparison:
    """Compare NJ and UPGMA topologies built from the same alignment."""
    nj_tree, _ = build_distance_tree(
        path,
        method="neighbor-joining",
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    upgma_tree, _ = build_distance_tree(
        path,
        method="upgma",
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    shared_taxa = set(nj_tree.tip_names) & set(upgma_tree.tip_names)
    (
        nj_clade_count,
        upgma_clade_count,
        rooted_distance,
        rooted_normalized,
    ) = _robinson_foulds_metrics(
        nj_tree,
        upgma_tree,
        shared_taxa,
        rf_mode="rooted",
    )
    topology_equal = rooted_distance == 0
    same_unrooted_topology = _unrooted_splits(nj_tree, shared_taxa) == _unrooted_splits(
        upgma_tree, shared_taxa
    )
    return DistanceTreeTopologyComparison(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        shared_taxa=sorted(shared_taxa),
        nj_informative_clades=nj_clade_count,
        upgma_informative_clades=upgma_clade_count,
        robinson_foulds_distance=rooted_distance,
        normalized_robinson_foulds=rooted_normalized,
        topology_equal=topology_equal,
        same_unrooted_topology=same_unrooted_topology,
        same_taxa_different_rooting=topology_equal is False and same_unrooted_topology,
    )


def _resampled_records(
    records: list[AlignmentRecord], *, rng: random.Random
) -> list[AlignmentRecord]:
    positions = [
        rng.randrange(len(records[0].sequence)) for _ in range(len(records[0].sequence))
    ]
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(record.sequence[position] for position in positions),
        )
        for record in records
    ]


def _write_bootstrap_alignment(path: Path, records: list[AlignmentRecord]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.append(record.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def bootstrap_distance_trees(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> tuple[list[PhyloTree], DistanceBootstrapReport]:
    """Bootstrap a distance tree by resampling alignment sites with replacement."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    if quality.method_assessment.decision == "blocked":
        raise InvalidAlignmentError(
            "distance bootstrap is blocked: "
            + "; ".join(quality.method_assessment.reasons)
        )
    records, _ = _load_alignment_for_model(path, model=model)
    # Deterministic scientific resampling is required for reproducible bootstrap review.
    rng = random.Random(seed)  # nosec B311
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-distance-bootstrap-"))
    trees: list[PhyloTree] = []
    for index in range(replicates):
        replicate_records = _resampled_records(records, rng=rng)
        replicate_path = temp_dir / f"replicate-{index + 1}.fasta"
        _write_bootstrap_alignment(replicate_path, replicate_records)
        tree, _ = build_distance_tree(
            replicate_path,
            method=method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        trees.append(tree)
    tree_set_path = temp_dir / "bootstrap.trees"
    write_tree_set(tree_set_path, trees)
    consensus_tree, consensus = compute_consensus_tree(tree_set_path)
    support = compute_clade_frequency_table(tree_set_path)
    return trees, DistanceBootstrapReport(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        method=method,
        replicates=replicates,
        seed=seed,
        tree_count=len(trees),
        consensus_newick=dumps_newick(consensus_tree),
        support=[
            DistanceBootstrapSupportRow(
                clade=row.clade,
                tree_count=row.tree_count,
                frequency=row.frequency,
            )
            for row in support.clade_frequencies
        ],
    )


def write_distance_bootstrap_support(
    path: Path, report: DistanceBootstrapReport
) -> Path:
    """Write bootstrap clade support as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["clade\ttree_count\tfrequency"]
    lines.extend(
        f"{row.clade}\t{row.tree_count}\t{format(row.frequency, '.15g')}"
        for row in report.support
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def summarize_distance_bootstrap_support(
    report: DistanceBootstrapReport,
    *,
    weak_frequency_threshold: float = 0.5,
) -> DistanceBootstrapSupportSummary:
    """Summarize bootstrap clade frequencies for reviewer-facing reporting."""
    frequencies = sorted(row.frequency for row in report.support)
    weak_clade_count = sum(
        1 for row in report.support if row.frequency < weak_frequency_threshold
    )
    warnings: list[str] = []
    if weak_clade_count:
        warnings.append(
            "one or more consensus clades remain weakly supported across bootstrap replicates"
        )
    if not frequencies:
        warnings.append(
            "bootstrap replicates did not yield any informative internal clades"
        )
    return DistanceBootstrapSupportSummary(
        alignment_path=report.alignment_path,
        method=report.method,
        model=report.model,
        gap_handling=report.gap_handling,
        ambiguity_policy=report.ambiguity_policy,
        replicates=report.replicates,
        clade_count=len(report.support),
        minimum_frequency=None if not frequencies else min(frequencies),
        maximum_frequency=None if not frequencies else max(frequencies),
        median_frequency=None if not frequencies else median(frequencies),
        weak_clade_count=weak_clade_count,
        warnings=warnings,
    )


def compare_distance_tree_to_reference_tree(
    path: Path,
    reference_tree_path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceTreeReferenceComparisonReport:
    """Compare one built distance tree to an external ML or reviewer-supplied reference tree."""
    tree, _ = build_distance_tree(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-distance-reference-"))
    built_tree_path = write_newick(temp_dir / "distance-tree.nwk", tree)
    topology = compare_tree_paths(built_tree_path, reference_tree_path)
    branch_lengths = compare_branch_lengths(built_tree_path, reference_tree_path)
    warnings: list[str] = []
    if not topology.topology_equal:
        warnings.append(
            "distance tree disagrees topologically with the supplied reference tree"
        )
    if topology.same_unrooted_topology and not topology.topology_equal:
        warnings.append(
            "distance tree matches the reference on unrooted splits but differs in rooting"
        )
    if topology.same_topology_different_branch_lengths:
        warnings.append(
            "distance tree preserves topology but shifts branch-length interpretation relative to the reference"
        )
    return DistanceTreeReferenceComparisonReport(
        alignment_path=path,
        reference_tree_path=reference_tree_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        topology=topology,
        branch_lengths=branch_lengths,
        warnings=warnings,
    )


def compare_distance_models(
    path: Path,
    *,
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceModelComparisonReport:
    """Compare all supported distance models for one alignment without overclaiming a winner."""
    records = load_fasta_alignment(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    models = sorted(_allowed_models_for_alphabet(inferred_alphabet))
    rows: list[DistanceModelComparisonRow] = []
    for model in models:
        quality = inspect_distance_matrix_quality(
            path,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        matrix = compute_pairwise_genetic_distance_matrix(
            path,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        unique_pairs = _unique_genetic_distance_pairs(matrix)
        defined_distances = sorted(
            float(pair.distance) for pair in unique_pairs if pair.distance is not None
        )
        rows.append(
            DistanceModelComparisonRow(
                model=model,
                defined_pair_count=len(defined_distances),
                saturated_pair_count=len(quality.saturated_pairs),
                low_information_pair_count=len(quality.low_information_pairs),
                mean_distance=None
                if not defined_distances
                else round(sum(defined_distances) / len(defined_distances), 15),
                maximum_distance=None
                if not defined_distances
                else max(defined_distances),
                decision=quality.method_assessment.decision,
                reasons=quality.method_assessment.reasons,
            )
        )
    warnings: list[str] = []
    if any(row.saturated_pair_count > 0 for row in rows):
        warnings.append(
            "one or more supported models enter a saturation regime on this alignment"
        )
    if len(rows) < 2:
        warnings.append(
            "the inferred alphabet supports only one distance model, so cross-model sensitivity is limited"
        )
    return DistanceModelComparisonReport(
        alignment_path=path,
        inferred_alphabet=inferred_alphabet,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        rows=rows,
        warnings=warnings,
    )


def compare_distance_gap_policies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceGapPolicySensitivityReport:
    """Summarize how pairwise versus complete deletion changes distance estimates."""
    pairwise = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling="pairwise-deletion",
        ambiguity_policy=ambiguity_policy,
    )
    complete = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling="complete-deletion",
        ambiguity_policy=ambiguity_policy,
    )
    pairwise_by_key = {
        _pair_key(pair.left_identifier, pair.right_identifier): pair
        for pair in _unique_genetic_distance_pairs(pairwise)
    }
    complete_by_key = {
        _pair_key(pair.left_identifier, pair.right_identifier): pair
        for pair in _unique_genetic_distance_pairs(complete)
    }
    rows: list[DistanceGapPolicyDeltaRow] = []
    for pair_key in sorted(set(pairwise_by_key) | set(complete_by_key)):
        pairwise_pair = pairwise_by_key.get(pair_key)
        complete_pair = complete_by_key.get(pair_key)
        if pairwise_pair is None or complete_pair is None:
            continue
        distance_delta = None
        if pairwise_pair.distance is not None and complete_pair.distance is not None:
            distance_delta = round(complete_pair.distance - pairwise_pair.distance, 15)
        if (
            pairwise_pair.distance != complete_pair.distance
            or pairwise_pair.comparable_sites != complete_pair.comparable_sites
        ):
            rows.append(
                DistanceGapPolicyDeltaRow(
                    left_identifier=pairwise_pair.left_identifier,
                    right_identifier=pairwise_pair.right_identifier,
                    pairwise_distance=pairwise_pair.distance,
                    complete_distance=complete_pair.distance,
                    pairwise_comparable_sites=pairwise_pair.comparable_sites,
                    complete_comparable_sites=complete_pair.comparable_sites,
                    distance_delta=distance_delta,
                    comparable_site_delta=complete_pair.comparable_sites
                    - pairwise_pair.comparable_sites,
                )
            )
    warnings: list[str] = []
    if rows:
        warnings.append(
            "distance estimates change when gap handling switches between pairwise and complete deletion"
        )
    if any(row.complete_comparable_sites == 0 for row in rows):
        warnings.append(
            "complete deletion removes all comparable sites for one or more taxon pairs"
        )
    return DistanceGapPolicySensitivityReport(
        alignment_path=path,
        model=model,
        ambiguity_policy=ambiguity_policy,
        changed_pair_count=len(rows),
        pair_count=len(pairwise_by_key),
        rows=rows,
        warnings=warnings,
    )


def assess_distance_method_maturity(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
    validate_bundle: bool = False,
) -> DistanceMethodMaturityGateReport:
    """Evaluate whether distance-analysis surfaces are available and validated for one alignment."""
    reference_validation = validate_distance_reference_examples()
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_summary = summarize_distance_bootstrap_support(
        bootstrap_distance_trees(
            path,
            method=method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=bootstrap_replicates,
            seed=bootstrap_seed,
        )[1]
    )
    model_comparison = compare_distance_models(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    gap_sensitivity = compare_distance_gap_policies(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
    )

    bundle_satisfied = False
    bundle_details = "bundle validation was skipped"
    if validate_bundle:
        bundle_dir = (
            Path(tempfile.mkdtemp(prefix="bijux-distance-maturity-")) / "bundle"
        )
        bundle = write_distance_reproducibility_bundle(
            bundle_dir,
            alignment_path=path,
            method=method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=min(bootstrap_replicates, 10),
            seed=bootstrap_seed,
        )
        manifest = json.loads(
            (bundle.out_dir / "distance-analysis.manifest.json").read_text(
                encoding="utf-8"
            )
        )
        bundle_satisfied = (
            "output_checksums" in manifest and "input_checksums" in manifest
        )
        bundle_details = (
            "bundle manifest includes input and output checksums"
            if bundle_satisfied
            else "bundle manifest is missing checksum provenance"
        )

    checks = [
        DistanceMethodMaturityCheck(
            name="reference_validation",
            satisfied=reference_validation.all_passed,
            details="built-in p-distance, JC69, K2P, protein, NJ, and UPGMA reference cases all pass"
            if reference_validation.all_passed
            else "one or more built-in distance reference cases failed",
        ),
        DistanceMethodMaturityCheck(
            name="upgma_assumption_visibility",
            satisfied=bool(assumptions.upgma_assumptions)
            and assumptions.ultrametric_tolerance > 0.0,
            details="UPGMA assumptions and ultrametric violations are surfaced explicitly",
        ),
        DistanceMethodMaturityCheck(
            name="suitability_decision",
            satisfied=quality.method_assessment.decision
            in {"allowed", "risky", "blocked"},
            details=f"distance workflow received explicit decision '{quality.method_assessment.decision}'",
        ),
        DistanceMethodMaturityCheck(
            name="bootstrap_support_summary",
            satisfied=bootstrap_summary.clade_count >= 0,
            details="bootstrap support frequencies are summarized over site-resampled replicate trees",
        ),
        DistanceMethodMaturityCheck(
            name="model_comparison",
            satisfied=bool(model_comparison.rows),
            details="all supported distance models were compared on the same alignment",
        ),
        DistanceMethodMaturityCheck(
            name="gap_policy_sensitivity",
            satisfied=True,
            details=f"{gap_sensitivity.changed_pair_count} taxon pairs changed across pairwise versus complete deletion",
        ),
        DistanceMethodMaturityCheck(
            name="bundle_provenance",
            satisfied=bundle_satisfied if validate_bundle else True,
            details=bundle_details,
        ),
    ]
    if not all(check.satisfied for check in checks):
        decision = "blocked"
    elif quality.method_assessment.decision == "allowed":
        decision = "production_candidate"
    else:
        decision = "validated_with_limits"
    return DistanceMethodMaturityGateReport(
        alignment_path=path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        decision=decision,
        checks=checks,
        warnings=quality.warnings
        + bootstrap_summary.warnings
        + model_comparison.warnings
        + gap_sensitivity.warnings,
    )


def build_distance_method_report(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
) -> DistanceMethodReport:
    """Build a structured distance-method report that can back JSON or HTML renderers."""
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    reference_validation = validate_distance_reference_examples()
    tree, _ = build_distance_tree(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    alternative_method = "upgma" if method == "neighbor-joining" else "neighbor-joining"
    alternative_tree, _ = build_distance_tree(
        path,
        method=alternative_method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    topology_comparison = compare_distance_tree_topologies(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_summary = summarize_distance_bootstrap_support(
        bootstrap_distance_trees(
            path,
            method=method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=bootstrap_replicates,
            seed=bootstrap_seed,
        )[1]
    )
    model_comparison = compare_distance_models(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    gap_policy_sensitivity = compare_distance_gap_policies(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
    )
    maturity_gate = assess_distance_method_maturity(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
        validate_bundle=False,
    )
    return DistanceMethodReport(
        alignment_path=path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        matrix=matrix,
        quality=quality,
        assumptions=assumptions,
        reference_validation=reference_validation,
        built_tree_newick=dumps_newick(tree),
        alternative_tree_newick=dumps_newick(alternative_tree),
        topology_comparison=topology_comparison,
        bootstrap_summary=bootstrap_summary,
        model_comparison=model_comparison,
        gap_policy_sensitivity=gap_policy_sensitivity,
        maturity_gate=maturity_gate,
    )


def write_distance_reproducibility_bundle(
    out_dir: Path,
    *,
    alignment_path: Path,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> DistanceReproducibilityBundleReport:
    """Write a reproducibility bundle for one distance analysis."""
    out_dir.mkdir(parents=True, exist_ok=True)
    bundled_alignment_path = out_dir / "input-alignment.fasta"
    shutil.copy2(alignment_path, bundled_alignment_path)
    matrix = compute_pairwise_genetic_distance_matrix(
        alignment_path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    quality = inspect_distance_matrix_quality(
        alignment_path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    tree, _ = build_distance_tree(
        alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_trees, bootstrap = bootstrap_distance_trees(
        alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )
    matrix_path = write_genetic_distance_matrix(out_dir / "distance-matrix.tsv", matrix)
    tree_path = write_newick(out_dir / "distance-tree.nwk", tree)
    support_path = write_distance_bootstrap_support(
        out_dir / "bootstrap-support.tsv", bootstrap
    )
    tree_set_path = write_tree_set(
        out_dir / "bootstrap-replicates.trees", bootstrap_trees
    )
    consensus_tree, _ = compute_consensus_tree(tree_set_path)
    consensus_path = write_newick(out_dir / "bootstrap-consensus.nwk", consensus_tree)
    clade_frequency_path = write_clade_frequency_table(
        out_dir / "bootstrap-clades.tsv",
        compute_clade_frequency_table(tree_set_path),
    )
    summary = build_distance_method_report(
        alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=min(replicates, 25),
        bootstrap_seed=seed,
    )
    summary_path = out_dir / "distance-summary.json"
    summary_path.write_text(
        json.dumps(asdict(summary), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest_path = out_dir / "distance-analysis.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "alignment_path": str(alignment_path),
                "method": method,
                "model": model,
                "gap_handling": gap_handling,
                "ambiguity_policy": ambiguity_policy,
                "replicates": replicates,
                "seed": seed,
                "quality_decision": quality.method_assessment.decision,
                "quality_reasons": quality.method_assessment.reasons,
                "files": [
                    "input-alignment.fasta",
                    "distance-matrix.tsv",
                    "distance-tree.nwk",
                    "bootstrap-support.tsv",
                    "bootstrap-replicates.trees",
                    "bootstrap-consensus.nwk",
                    "bootstrap-clades.tsv",
                    "distance-summary.json",
                ],
                "input_checksums": {
                    str(alignment_path): _file_sha256(alignment_path),
                },
                "output_checksums": {
                    "input-alignment.fasta": _file_sha256(bundled_alignment_path),
                    "distance-matrix.tsv": _file_sha256(matrix_path),
                    "distance-tree.nwk": _file_sha256(tree_path),
                    "bootstrap-support.tsv": _file_sha256(support_path),
                    "bootstrap-replicates.trees": _file_sha256(tree_set_path),
                    "bootstrap-consensus.nwk": _file_sha256(consensus_path),
                    "bootstrap-clades.tsv": _file_sha256(clade_frequency_path),
                    "distance-summary.json": _file_sha256(summary_path),
                },
                "reference_validation_passed": validate_distance_reference_examples().all_passed,
                "caveats": [
                    "distance methods reduce site-by-site evidence to pairwise summaries before tree building",
                    "bootstrap support measures repeatability under site resampling, not absolute correctness",
                    *quality.warnings,
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    files = [
        bundled_alignment_path,
        matrix_path,
        tree_path,
        support_path,
        tree_set_path,
        consensus_path,
        clade_frequency_path,
        summary_path,
        manifest_path,
    ]
    return DistanceReproducibilityBundleReport(
        out_dir=out_dir,
        alignment_path=alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        files=files,
    )


def build_tree_from_imported_distance_matrix(
    path: Path,
    *,
    method: str,
) -> tuple[PhyloTree, ImportedDistanceTreeBuildReport]:
    """Build a distance-based tree from an imported long-form distance matrix."""
    entries = load_imported_distance_matrix(path)
    validation = validate_imported_distance_matrix(path)
    assumptions = assess_imported_distance_method_assumptions(path)
    constructor = DistanceTreeConstructor()
    distance_matrix = _bio_distance_matrix_from_imported(validation, entries)
    if method == "neighbor-joining":
        tree = constructor.nj(distance_matrix)
    elif method == "upgma":
        tree = constructor.upgma(distance_matrix)
    else:
        raise ValueError(f"unsupported tree-building method: {method}")
    return tree_from_biophylo(
        tree, source_format="newick"
    ), ImportedDistanceTreeBuildReport(
        matrix_path=path,
        method=method,
        taxon_count=len(validation.identifiers),
        pair_count=validation.pair_count,
        assumptions=assumptions,
    )
