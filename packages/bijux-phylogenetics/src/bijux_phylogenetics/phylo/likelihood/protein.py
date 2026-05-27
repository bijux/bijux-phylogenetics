from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.runtime.errors import AlignmentTaxonMismatchError
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

PROTEIN_STATE_ORDER = (
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
)
PROTEIN_STATE_INDEX = {
    state: index for index, state in enumerate(PROTEIN_STATE_ORDER)
}
UNIFORM_PROTEIN_ROOT_PRIOR = numpy.full(
    len(PROTEIN_STATE_ORDER),
    1.0 / len(PROTEIN_STATE_ORDER),
    dtype=float,
)
PROTEIN_GAP_CHARACTER = "-"
PROTEIN_MISSING_CHARACTER = "?"
_PROTEIN_ALLOWED_POLICY_VALUES = frozenset({"treat-as-missing", "reject"})


def normalize_unambiguous_protein_records(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> list[AlignmentRecord]:
    """Uppercase aligned protein records and apply explicit gap and missing policy."""
    validated_gap_policy = validate_protein_gap_policy(
        gap_policy,
        model_name=model_name,
    )
    validated_missing_policy = validate_protein_missing_policy(
        missing_policy,
        model_name=model_name,
    )
    normalized_records: list[AlignmentRecord] = []
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_states: set[str] = set()
        for state in normalized_sequence:
            if state in PROTEIN_STATE_INDEX:
                continue
            if state == PROTEIN_GAP_CHARACTER:
                if validated_gap_policy == "reject":
                    raise InvalidAlignmentError(
                        f"{model_name} likelihood gap policy rejects '-' in record '{record.identifier}'"
                    )
                continue
            if state == PROTEIN_MISSING_CHARACTER:
                if validated_missing_policy == "reject":
                    raise InvalidAlignmentError(
                        f"{model_name} likelihood missing-state policy rejects '?' in record '{record.identifier}'"
                    )
                continue
            invalid_states.add(state)
        if invalid_states:
            joined_states = ", ".join(sorted(invalid_states))
            raise InvalidAlignmentError(
                f"{model_name} likelihood currently requires unambiguous amino-acid states plus explicit gap/missing policy; "
                f"record '{record.identifier}' contains {joined_states}"
            )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def validate_protein_gap_policy(
    gap_policy: str,
    *,
    model_name: str,
) -> str:
    if gap_policy not in _PROTEIN_ALLOWED_POLICY_VALUES:
        raise ValueError(
            f"{model_name} protein gap policy must be one of {sorted(_PROTEIN_ALLOWED_POLICY_VALUES)}"
        )
    return gap_policy


def validate_protein_missing_policy(
    missing_policy: str,
    *,
    model_name: str,
) -> str:
    if missing_policy not in _PROTEIN_ALLOWED_POLICY_VALUES:
        raise ValueError(
            f"{model_name} protein missing-state policy must be one of {sorted(_PROTEIN_ALLOWED_POLICY_VALUES)}"
        )
    return missing_policy


def protein_leaf_likelihood_vector(
    states_by_taxon: dict[str, str],
    *,
    model_name: str,
    node_name: str | None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> numpy.ndarray:
    """Return a one-hot or all-states-allowed vector for one protein tip state."""
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires named tree tips for alignment lookup"
        )
    state = states_by_taxon[node_name]
    if state in PROTEIN_STATE_INDEX:
        vector = numpy.zeros(len(PROTEIN_STATE_ORDER), dtype=float)
        vector[PROTEIN_STATE_INDEX[state]] = 1.0
        return vector
    if state == PROTEIN_GAP_CHARACTER:
        if validate_protein_gap_policy(gap_policy, model_name=model_name) == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood gap policy rejects '-' at taxon '{node_name}'"
            )
        return numpy.ones(len(PROTEIN_STATE_ORDER), dtype=float)
    if state == PROTEIN_MISSING_CHARACTER:
        if (
            validate_protein_missing_policy(
                missing_policy,
                model_name=model_name,
            )
            == "reject"
        ):
            raise InvalidAlignmentError(
                f"{model_name} likelihood missing-state policy rejects '?' at taxon '{node_name}'"
            )
        return numpy.ones(len(PROTEIN_STATE_ORDER), dtype=float)
    raise InvalidAlignmentError(
        f"{model_name} likelihood encountered unsupported amino-acid state '{state}' at taxon '{node_name}'"
    )
