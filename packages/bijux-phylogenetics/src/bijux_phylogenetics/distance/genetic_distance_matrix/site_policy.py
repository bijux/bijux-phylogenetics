from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.alignment import AlignmentRecord

from ..models import AmbiguityPolicy

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


@dataclass(frozen=True, slots=True)
class _SiteContribution:
    comparable: bool
    mismatch_weight: float
    transition_weight: float
    ag_transition_weight: float
    ct_transition_weight: float
    transversion_weight: float
    ambiguous: bool


def _normalize_residue(residue: str) -> str:
    upper = residue.upper()
    if upper == "U":
        return "T"
    return upper


def _states_for_symbol(symbol: str, *, alphabet: str) -> set[str] | None:
    normalized = _normalize_residue(symbol)
    if normalized in _MISSING_OR_GAP:
        return set()
    if alphabet == "protein":
        return _PROTEIN_AMBIGUITIES.get(normalized)
    return _NUCLEOTIDE_AMBIGUITIES.get(normalized)


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
