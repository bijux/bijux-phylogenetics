from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import math
from pathlib import Path
from statistics import median

from Bio.Phylo.TreeConstruction import DistanceMatrix

from bijux_phylogenetics.compare.topology import (
    BranchLengthComparisonReport,
    TreeComparisonReport,
)
from bijux_phylogenetics.core.alignment import AlignmentRecord, DnaBinAlignment
from bijux_phylogenetics.core.clade_sets import (
    informative_rooted_clades,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    UnsupportedDistanceTreeMethodError,
)
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_dna_bin_alignment,
    load_fasta_alignment,
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
_PACKAGE_ROOT = Path(__file__).resolve().parents[3]


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
    method_policy: DistanceTreeMethodPolicy
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
    method_policy: DistanceTreeMethodPolicy
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


@dataclass(frozen=True, slots=True)
class DistanceTreeMethodPolicy:
    """Stable support policy for one distance-tree method surface."""

    method: str
    supported: bool
    reference_surface: str | None
    support_scope: str
    summary: str
    limitations: list[str]


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
    method_policy: DistanceTreeMethodPolicy
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
    method_policy: DistanceTreeMethodPolicy
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


def list_distance_tree_method_policies() -> list[DistanceTreeMethodPolicy]:
    """Return the governed distance-tree method support policy surface."""
    return [
        DistanceTreeMethodPolicy(
            method="neighbor-joining",
            supported=True,
            reference_surface="ape::nj",
            support_scope="owned-runtime-and-live-ape-parity",
            summary="Neighbor-Joining is fully supported and covered by the governed live ape parity lane.",
            limitations=[
                "neighbor-joining remains a distance-summary method rather than a full likelihood inference",
            ],
        ),
        DistanceTreeMethodPolicy(
            method="upgma",
            supported=True,
            reference_surface=None,
            support_scope="owned-runtime",
            summary="UPGMA is supported as an owned runtime method for ultrametric reviewer workflows.",
            limitations=[
                "upgma assumes an ultrametric clock-like process and can misplace taxa when rates vary among lineages",
            ],
        ),
        DistanceTreeMethodPolicy(
            method="bionj",
            supported=False,
            reference_surface="ape::bionj",
            support_scope="explicitly-excluded",
            summary="BIONJ is explicitly out of scope for this round, and ape::bionj does not have an owned Bijux runtime or governed ape parity lane.",
            limitations=[
                "request neighbor-joining when a governed ape-parity distance-tree workflow is required in this round",
            ],
        ),
    ]


def resolve_distance_tree_method_policy(method: str) -> DistanceTreeMethodPolicy:
    """Resolve one distance-tree method name to its owned support policy."""
    normalized = method.strip().lower()
    for policy in list_distance_tree_method_policies():
        if policy.method == normalized:
            return policy
    supported = [
        policy.method
        for policy in list_distance_tree_method_policies()
        if policy.supported
    ]
    raise UnsupportedDistanceTreeMethodError(
        f"unsupported tree-building method '{method}'; supported methods are {', '.join(supported)}",
        details={
            "requested_method": method,
            "supported_methods": supported,
        },
    )


def _require_supported_distance_tree_method(method: str) -> DistanceTreeMethodPolicy:
    """Resolve one distance-tree method and reject explicitly excluded methods."""
    policy = resolve_distance_tree_method_policy(method)
    if policy.supported:
        return policy
    supported = [
        row.method for row in list_distance_tree_method_policies() if row.supported
    ]
    raise UnsupportedDistanceTreeMethodError(
        policy.summary,
        details={
            "requested_method": policy.method,
            "supported_methods": supported,
            "excluded_method": policy.method,
            "reference_surface": policy.reference_surface,
        },
    )


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
    counts = dict.fromkeys("ACGT", 0)
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
        (pi_r * pi_y) - ((pi_a * pi_g * pi_y) / pi_r) - ((pi_c * pi_t * pi_r) / pi_y)
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
    normalized_model = _normalize_distance_model(model)
    if normalized_model != "amino-acid-p-distance":
        matrix = load_dna_bin_alignment(path, normalize_uracil=True)
        if normalized_model not in _allowed_models_for_alphabet(matrix.source_alphabet):
            raise InvalidAlignmentError(
                f"distance model '{normalized_model}' is not supported for inferred alphabet '{matrix.source_alphabet}'"
            )
        records = [
            AlignmentRecord(
                identifier=record.identifier,
                sequence=record.sequence.upper(),
            )
            for record in matrix.records
        ]
        return records, matrix.source_alphabet

    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if normalized_model not in _allowed_models_for_alphabet(alphabet):
        raise InvalidAlignmentError(
            f"distance model '{normalized_model}' is not supported for inferred alphabet '{alphabet}'"
        )
    return records, alphabet


def compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    """Compute a deterministic nucleotide distance matrix from one DNAbin-compatible alignment."""
    model = _normalize_distance_model(model)
    if model not in _allowed_models_for_alphabet(alignment.source_alphabet):
        raise InvalidAlignmentError(
            f"distance model '{model}' is not supported for inferred alphabet '{alignment.source_alphabet}'"
        )
    if model == "amino-acid-p-distance":
        raise InvalidAlignmentError(
            "dnabin-compatible nucleotide distance loading does not support amino-acid distances"
        )
    if gap_handling not in {"pairwise-deletion", "complete-deletion"}:
        raise ValueError(f"unsupported gap handling mode: {gap_handling}")
    if ambiguity_policy not in {
        "ignore",
        "partial-match",
        "strict-mismatch",
        "report-only",
    }:
        raise ValueError(f"unsupported ambiguity policy: {ambiguity_policy}")

    records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=record.sequence.upper(),
        )
        for record in alignment.records
    ]
    alphabet = alignment.source_alphabet
    model_parameters: GeneticDistanceModelParameters | None = None
    warnings: list[str] = []
    if model in {"felsenstein-81", "tamura-nei-93"}:
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
        path=alignment.path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        inferred_alphabet=alphabet,
        alignment_length=alignment.alignment_length,
        identifiers=[record.identifier for record in records],
        model_parameters=model_parameters,
        warnings=warnings,
        pairs=pairs,
    )


def _pair_key(left_identifier: str, right_identifier: str) -> tuple[str, str]:
    return tuple(sorted((left_identifier, right_identifier)))

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


def assess_distance_method_assumptions_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    """Audit core distance-tree assumptions from one in-memory distance matrix."""
    distances = _build_alignment_distance_lookup(report)
    expected_pairs = (len(report.identifiers) * (len(report.identifiers) - 1)) // 2
    violations = _iter_ultrametric_violations(
        report.identifiers,
        distances,
        tolerance=ultrametric_tolerance,
    )
    warnings = list(report.warnings)
    if len(distances) < expected_pairs:
        warnings.append(
            "UPGMA clock-assumption auditing is incomplete because one or more pairwise distances are undefined under the selected model"
        )
    if violations:
        warnings.append(
            "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
        )
    return DistanceMethodAssumptionReport(
        source_path=report.path,
        source_kind="alignment",
        taxon_count=len(report.identifiers),
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


def _distance_lookup(report: GeneticDistanceMatrix) -> dict[tuple[str, str], float]:
    lookup: dict[tuple[str, str], float] = {}
    for identifier in report.identifiers:
        lookup[(identifier, identifier)] = 0.0
    for pair in report.pairs:
        if pair.distance is None:
            raise InvalidAlignmentError(
                "distance matrix contains undefined entries for: "
                f"{pair.left_identifier}/{pair.right_identifier}"
            )
        lookup[(pair.left_identifier, pair.right_identifier)] = float(pair.distance)
        lookup[(pair.right_identifier, pair.left_identifier)] = float(pair.distance)
    return lookup


def load_imported_distance_matrix(path: Path) -> list[ImportedDistanceEntry]:
    from .imported import load_imported_distance_matrix as load_impl

    return load_impl(path)


def assess_imported_distance_method_assumptions(
    path: Path,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    from .imported import (
        assess_imported_distance_method_assumptions as assess_imported_impl,
    )

    return assess_imported_impl(
        path, ultrametric_tolerance=ultrametric_tolerance
    )


def validate_imported_distance_matrix(path: Path) -> ImportedDistanceMatrixReport:
    from .imported import validate_imported_distance_matrix as validate_impl

    return validate_impl(path)


def inspect_imported_distance_matrix_quality(
    path: Path,
) -> ImportedDistanceMatrixQualityReport:
    from .imported import (
        inspect_imported_distance_matrix_quality as inspect_quality_impl,
    )

    return inspect_quality_impl(path)


def build_tree_from_imported_distance_matrix(
    path: Path,
    *,
    method: str,
) -> tuple[PhyloTree, ImportedDistanceTreeBuildReport]:
    from .imported import (
        build_tree_from_imported_distance_matrix as build_imported_tree_impl,
    )

    return build_imported_tree_impl(path, method=method)


def build_distance_tree(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    from .tree_inference import build_distance_tree as build_tree_impl

    return build_tree_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def build_distance_tree_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    method: str,
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    from .tree_inference import (
        build_distance_tree_from_genetic_distance_matrix as build_from_matrix_impl,
    )

    return build_from_matrix_impl(report, method=method)


def compare_distance_tree_topologies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceTreeTopologyComparison:
    from .tree_inference import compare_distance_tree_topologies as compare_trees_impl

    return compare_trees_impl(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


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
    from .tree_inference import bootstrap_distance_trees as bootstrap_impl

    return bootstrap_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )


def write_distance_bootstrap_support(
    path: Path, report: DistanceBootstrapReport
) -> Path:
    from .tree_inference import (
        write_distance_bootstrap_support as write_support_impl,
    )

    return write_support_impl(path, report)


def summarize_distance_bootstrap_support(
    report: DistanceBootstrapReport,
    *,
    weak_frequency_threshold: float = 0.5,
) -> DistanceBootstrapSupportSummary:
    from .tree_inference import (
        summarize_distance_bootstrap_support as summarize_support_impl,
    )

    return summarize_support_impl(
        report, weak_frequency_threshold=weak_frequency_threshold
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
    from .tree_inference import (
        compare_distance_tree_to_reference_tree as compare_reference_impl,
    )

    return compare_reference_impl(
        path,
        reference_tree_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )



def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    """Compute a deterministic pairwise genetic distance matrix for an aligned dataset."""
    model = _normalize_distance_model(model)
    if model != "amino-acid-p-distance":
        matrix = load_dna_bin_alignment(path, normalize_uracil=True)
        return compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
            matrix,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )

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
                model_parameters=None,
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
        model_parameters=None,
        warnings=[],
        pairs=pairs,
    )


def validate_distance_reference_examples() -> DistanceReferenceValidationReport:
    """Validate core distance examples against durable reference expectations."""
    examples = [
        {
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
            "path": _PACKAGE_ROOT
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
        _PACKAGE_ROOT / "tests/fixtures/metadata/example_distance_matrix_ultrametric.tsv"
    )
    nj_tree, _ = build_tree_from_imported_distance_matrix(
        nj_matrix_path, method="neighbor-joining"
    )
    nj_expected_clades = ["A|B"]
    nj_observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in informative_rooted_clades(nj_tree, set(nj_tree.tip_names))
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
        _PACKAGE_ROOT / "tests/fixtures/metadata/example_distance_matrix_ultrametric.tsv"
    )
    upgma_tree, _ = build_tree_from_imported_distance_matrix(
        upgma_matrix_path, method="upgma"
    )
    expected_clades = ["A|B", "C|D"]
    observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in informative_rooted_clades(upgma_tree, set(upgma_tree.tip_names))
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
                ""
                if normalized_distance is None
                else format(normalized_distance, ".15g")
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
                    ""
                    if normalized_distance is None
                    else format(normalized_distance, ".15g"),
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


def compare_distance_models(
    path: Path,
    *,
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceModelComparisonReport:
    from .reporting import compare_distance_models as compare_models_impl

    return compare_models_impl(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def compare_distance_gap_policies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceGapPolicySensitivityReport:
    from .reporting import (
        compare_distance_gap_policies as compare_gap_policies_impl,
    )

    return compare_gap_policies_impl(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
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
    from .reporting import (
        assess_distance_method_maturity as assess_maturity_impl,
    )

    return assess_maturity_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
        validate_bundle=validate_bundle,
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
    from .reporting import build_distance_method_report as build_report_impl

    return build_report_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
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
    from .reporting import (
        write_distance_reproducibility_bundle as write_bundle_impl,
    )

    return write_bundle_impl(
        out_dir,
        alignment_path=alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )
