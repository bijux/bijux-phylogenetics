from __future__ import annotations

from collections import defaultdict

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import TaxonIdentityAudit, TaxonLabelPair


def _identity_key(label: str) -> str:
    return "".join(character.lower() for character in label if character.isalnum())


def _space_underscore_key(label: str) -> str:
    collapsed = label.strip().replace("_", " ")
    return " ".join(collapsed.split()).lower()


def _suspicious_near_duplicate(label: str, other: str) -> bool:
    left_key = _identity_key(label)
    right_key = _identity_key(other)
    minimum_length = min(len(left_key), len(right_key))
    if minimum_length < 5:
        return False
    distance = _levenshtein_distance(left_key, right_key)
    if distance == 0 or distance > 2:
        return False
    return distance / max(len(left_key), len(right_key)) <= 0.2


def _levenshtein_distance(left: str, right: str) -> int:
    if left == right:
        return 0
    if not left:
        return len(right)
    if not right:
        return len(left)
    previous = list(range(len(right) + 1))
    for index, left_char in enumerate(left, start=1):
        current = [index]
        for right_index, right_char in enumerate(right, start=1):
            insert_cost = current[right_index - 1] + 1
            delete_cost = previous[right_index] + 1
            substitute_cost = previous[right_index - 1] + (
                0 if left_char == right_char else 1
            )
            current.append(min(insert_cost, delete_cost, substitute_cost))
        previous = current
    return previous[-1]


def _sorted_pairs(pairs: set[tuple[str, str]]) -> list[TaxonLabelPair]:
    return [
        TaxonLabelPair(left_label=left, right_label=right)
        for left, right in sorted(pairs)
    ]


def inspect_tree_taxon_identity(tree: PhyloTree) -> TaxonIdentityAudit:
    """Audit tree tip labels for likely biological identity conflicts."""
    labels = sorted(tree.tip_names)
    spelling_variants: set[tuple[str, str]] = set()
    whitespace_variants: set[tuple[str, str]] = set()
    underscore_space_collisions: set[tuple[str, str]] = set()
    case_collisions: set[tuple[str, str]] = set()
    suspicious_near_duplicates: set[tuple[str, str]] = set()

    by_identity: dict[str, list[str]] = defaultdict(list)
    by_space_key: dict[str, list[str]] = defaultdict(list)
    by_casefold: dict[str, list[str]] = defaultdict(list)
    for label in labels:
        by_identity[_identity_key(label)].append(label)
        by_space_key[_space_underscore_key(label)].append(label)
        by_casefold[label.casefold()].append(label)

    for grouped in by_identity.values():
        if len(set(grouped)) < 2:
            continue
        ordered = sorted(set(grouped))
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                if left.replace(" ", "") != right.replace(" ", ""):
                    spelling_variants.add((left, right))

    for grouped in by_space_key.values():
        if len(set(grouped)) < 2:
            continue
        ordered = sorted(set(grouped))
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                if left != right:
                    if (
                        " ".join(left.split()).lower()
                        == " ".join(right.split()).lower()
                    ):
                        whitespace_variants.add((left, right))
                    if (
                        left.replace("_", " ").casefold()
                        == right.replace("_", " ").casefold()
                    ):
                        underscore_space_collisions.add((left, right))

    for grouped in by_casefold.values():
        if len(set(grouped)) < 2:
            continue
        ordered = sorted(set(grouped))
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                case_collisions.add((left, right))

    for index, left in enumerate(labels):
        for right in labels[index + 1 :]:
            if left.casefold() == right.casefold():
                continue
            if _identity_key(left) == _identity_key(right):
                continue
            if _suspicious_near_duplicate(left, right):
                suspicious_near_duplicates.add((left, right))

    return TaxonIdentityAudit(
        spelling_variants=_sorted_pairs(spelling_variants),
        whitespace_variants=_sorted_pairs(whitespace_variants),
        underscore_space_collisions=_sorted_pairs(underscore_space_collisions),
        case_collisions=_sorted_pairs(case_collisions),
        suspicious_near_duplicates=_sorted_pairs(suspicious_near_duplicates),
    )
