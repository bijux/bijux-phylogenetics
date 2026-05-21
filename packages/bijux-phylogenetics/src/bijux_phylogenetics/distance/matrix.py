from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceMatrix

from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.matrix import load_dna_bin_alignment
from bijux_phylogenetics.phylo.alignment import AlignmentRecord, DnaBinAlignment
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .models import (
    AmbiguityPolicy,
    DistanceModel,
    GapHandlingMode,
    GeneticDistanceMatrix,
    GeneticDistanceModelParameters,
    PairwiseGeneticDistance,
)
from .shared import (
    _allowed_models_for_alphabet,
    _normalize_distance_model,
    _pair_key,
)
from .genetic_distance_matrix.site_policy import (
    _complete_deletion_positions,
    _site_contribution,
    _states_for_symbol,
)


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
