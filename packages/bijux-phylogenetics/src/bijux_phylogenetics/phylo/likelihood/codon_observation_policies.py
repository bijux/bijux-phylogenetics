from __future__ import annotations

from itertools import product

import numpy

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

_CODON_DNA_BASES = frozenset({"A", "C", "G", "T"})
_CODON_EXPLICIT_MISSING_SYMBOLS = frozenset({"?"})
_CODON_GAP_SYMBOL = "-"
_SUPPORTED_CODON_OBSERVATION_POLICIES = frozenset(
    {"reject", "treat-as-missing", "ambiguity-vector"}
)
_CODON_AMBIGUITY_BASES_BY_SYMBOL = {
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "N": ("A", "C", "G", "T"),
    "R": ("A", "G"),
    "S": ("C", "G"),
    "V": ("A", "C", "G"),
    "W": ("A", "T"),
    "Y": ("C", "T"),
}


def validate_codon_observation_policy(
    observation_policy: str,
    *,
    owner_name: str,
) -> str:
    normalized_policy = observation_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_CODON_OBSERVATION_POLICIES:
        raise ValueError(
            f"{owner_name} observation_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_CODON_OBSERVATION_POLICIES))
        )
    return normalized_policy


def validate_codon_observation_state(
    codon: str,
    *,
    codon_index: int,
    genetic_code_name: str,
    owner_name: str,
    observation_policy: str,
    record_identifier: str,
    stop_codons: frozenset[str],
) -> None:
    validated_policy = validate_codon_observation_policy(
        observation_policy,
        owner_name=owner_name,
    )
    for symbol in codon:
        if (
            symbol not in _CODON_DNA_BASES
            and symbol not in _CODON_AMBIGUITY_BASES_BY_SYMBOL
            and symbol not in _CODON_EXPLICIT_MISSING_SYMBOLS
            and symbol != _CODON_GAP_SYMBOL
        ):
            raise InvalidAlignmentError(
                f"{owner_name} does not recognize codon observation state '{codon}' in record '{record_identifier}' at codon site {codon_index}"
            )
    if set(codon) <= _CODON_DNA_BASES and codon in stop_codons:
        raise InvalidAlignmentError(
            f"{owner_name} excludes stop codon states under the {genetic_code_name} code; "
            f"record '{record_identifier}' contains stop codon '{codon}' at codon site {codon_index}"
        )
    if validated_policy == "reject" and any(
        symbol in _CODON_EXPLICIT_MISSING_SYMBOLS or symbol == _CODON_GAP_SYMBOL
        for symbol in codon
    ):
        raise InvalidAlignmentError(
            f"{owner_name} observation policy 'reject' does not allow missing or gap codon state '{codon}' in record '{record_identifier}' at codon site {codon_index}"
        )
    if (
        validated_policy == "reject"
        and codon not in stop_codons
        and codon
        not in (
            _compatible_sense_codons(
                codon,
                state_index={codon: 0},
            )
        )
    ):
        raise InvalidAlignmentError(
            f"{owner_name} observation policy 'reject' requires resolved sense codons only; "
            f"record '{record_identifier}' contains '{codon}' at codon site {codon_index}"
        )


def resolve_codon_observation_leaf_vector(
    codon: str,
    *,
    model_name: str,
    node_name: str,
    observation_policy: str,
    state_order: tuple[str, ...],
    state_index: dict[str, int],
) -> numpy.ndarray:
    validated_policy = validate_codon_observation_policy(
        observation_policy,
        owner_name=f"{model_name} likelihood",
    )
    vector = numpy.zeros(len(state_order), dtype=float)
    if codon in state_index:
        vector[state_index[codon]] = 1.0
        return vector
    if any(
        symbol in _CODON_EXPLICIT_MISSING_SYMBOLS or symbol == _CODON_GAP_SYMBOL
        for symbol in codon
    ):
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow missing or gap codon state '{codon}' at taxon '{node_name}'"
            )
        vector[:] = 1.0
        return vector
    if validated_policy == "reject":
        raise InvalidAlignmentError(
            f"{model_name} likelihood observation policy 'reject' does not allow ambiguous codon state '{codon}' at taxon '{node_name}'"
        )
    if validated_policy == "treat-as-missing":
        vector[:] = 1.0
        return vector
    compatible_codons = _compatible_sense_codons(
        codon,
        state_index=state_index,
    )
    if not compatible_codons:
        raise InvalidAlignmentError(
            f"{model_name} likelihood ambiguity policy could not resolve '{codon}' to any supported sense codon at taxon '{node_name}'"
        )
    for compatible_codon in compatible_codons:
        vector[state_index[compatible_codon]] = 1.0
    return vector


def _compatible_sense_codons(
    codon: str,
    *,
    state_index: dict[str, int],
) -> tuple[str, ...]:
    allowed_bases_by_position: list[tuple[str, ...]] = []
    for symbol in codon:
        if symbol in _CODON_DNA_BASES:
            allowed_bases_by_position.append((symbol,))
            continue
        ambiguity_bases = _CODON_AMBIGUITY_BASES_BY_SYMBOL.get(symbol)
        if ambiguity_bases is None:
            return ()
        allowed_bases_by_position.append(ambiguity_bases)
    return tuple(
        candidate
        for candidate in (
            "".join(bases) for bases in product(*allowed_bases_by_position)
        )
        if candidate in state_index
    )


__all__ = [
    "resolve_codon_observation_leaf_vector",
    "validate_codon_observation_policy",
    "validate_codon_observation_state",
]
