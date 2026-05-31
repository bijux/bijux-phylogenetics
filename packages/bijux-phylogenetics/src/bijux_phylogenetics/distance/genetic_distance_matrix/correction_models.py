from __future__ import annotations

import math

from bijux_phylogenetics.phylo.alignment import AlignmentRecord

from ..models import DistanceModel, GeneticDistanceModelParameters
from .pairwise_summary import _PairSummary
from .site_policy import _states_for_symbol


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
