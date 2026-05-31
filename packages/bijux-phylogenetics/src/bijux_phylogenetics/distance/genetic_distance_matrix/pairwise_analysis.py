from __future__ import annotations

from ..models import AmbiguityPolicy, DistanceModel, GeneticDistanceModelParameters
from .correction_models import _distance_from_summary
from .pairwise_summary import _PairSummary
from .site_policy import _site_contribution, _states_for_symbol


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
