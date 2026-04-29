from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode


@dataclass(slots=True)
class TaxonomyReference:
    """Placeholder taxonomy reference used by future adapters."""

    authority: str
    identifier: str


@dataclass(frozen=True, slots=True)
class TaxonRename:
    """Explicit record of one taxon label normalization."""

    raw_label: str
    normalized_label: str


@dataclass(slots=True)
class TaxonNormalizationReport:
    """Result of explicit taxon renaming."""

    policy: str
    renamed_taxa: list[TaxonRename]


@dataclass(frozen=True, slots=True)
class UnsafeTaxonName:
    """One taxon label and the downstream-safety issues it triggers."""

    raw_label: str
    normalized_label: str
    reasons: list[str]


@dataclass(frozen=True, slots=True)
class TaxonNormalizationCollision:
    """Two or more raw labels that normalize to the same downstream label."""

    normalized_label: str
    raw_labels: list[str]


@dataclass(slots=True)
class TaxonSafetyReport:
    """Audit report for raw tree tip labels under an explicit normalization policy."""

    policy: str
    unsafe_taxa: list[UnsafeTaxonName]
    collisions: list[TaxonNormalizationCollision]


@dataclass(frozen=True, slots=True)
class TaxonLabelPair:
    left_label: str
    right_label: str


@dataclass(slots=True)
class TaxonIdentityAudit:
    """Audit likely identity conflicts among tip labels."""

    spelling_variants: list[TaxonLabelPair]
    whitespace_variants: list[TaxonLabelPair]
    underscore_space_collisions: list[TaxonLabelPair]
    case_collisions: list[TaxonLabelPair]
    suspicious_near_duplicates: list[TaxonLabelPair]


def _normalize_label(label: str, *, policy: str) -> str:
    if policy == "spaces-to-underscores":
        return "_".join(label.split())
    raise ValueError(f"unsupported taxon normalization policy: {policy}")


def _unsafe_reasons(label: str, *, normalized_label: str) -> list[str]:
    reasons: list[str] = []
    if any(character.isspace() for character in label):
        reasons.append("contains whitespace")
    if "'" in label or '"' in label:
        reasons.append("contains quote characters")
    if "/" in label or "\\" in label:
        reasons.append("contains slash characters")
    if not normalized_label:
        reasons.append("normalizes to an empty label")
    return reasons


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
            substitute_cost = previous[right_index - 1] + (0 if left_char == right_char else 1)
            current.append(min(insert_cost, delete_cost, substitute_cost))
        previous = current
    return previous[-1]


def _sorted_pairs(pairs: set[tuple[str, str]]) -> list[TaxonLabelPair]:
    return [TaxonLabelPair(left_label=left, right_label=right) for left, right in sorted(pairs)]


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
                    if " ".join(left.split()).lower() == " ".join(right.split()).lower():
                        whitespace_variants.add((left, right))
                    if left.replace("_", " ").casefold() == right.replace("_", " ").casefold():
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


def inspect_tree_taxa_safety(tree: PhyloTree, *, policy: str) -> TaxonSafetyReport:
    """Inspect tree tip labels for downstream-unsafe names and normalization collisions."""
    grouped_labels: dict[str, list[str]] = defaultdict(list)
    unsafe_taxa: list[UnsafeTaxonName] = []

    for label in tree.tip_names:
        normalized_label = _normalize_label(label, policy=policy)
        grouped_labels[normalized_label].append(label)
        reasons = _unsafe_reasons(label, normalized_label=normalized_label)
        if reasons:
            unsafe_taxa.append(
                UnsafeTaxonName(
                    raw_label=label,
                    normalized_label=normalized_label,
                    reasons=reasons,
                )
            )

    collisions = [
        TaxonNormalizationCollision(
            normalized_label=normalized_label,
            raw_labels=sorted(raw_labels),
        )
        for normalized_label, raw_labels in sorted(grouped_labels.items())
        if len(set(raw_labels)) > 1
    ]

    collision_labels = {label for collision in collisions for label in collision.raw_labels}
    if collision_labels:
        augmented_unsafe_taxa: list[UnsafeTaxonName] = []
        for entry in unsafe_taxa:
            reasons = list(entry.reasons)
            if entry.raw_label in collision_labels:
                reasons.append("collides with another label after normalization")
            augmented_unsafe_taxa.append(
                UnsafeTaxonName(
                    raw_label=entry.raw_label,
                    normalized_label=entry.normalized_label,
                    reasons=reasons,
                )
            )
        unsafe_taxa = augmented_unsafe_taxa

        existing = {entry.raw_label for entry in unsafe_taxa}
        for label in sorted(collision_labels - existing):
            normalized_label = _normalize_label(label, policy=policy)
            unsafe_taxa.append(
                UnsafeTaxonName(
                    raw_label=label,
                    normalized_label=normalized_label,
                    reasons=["collides with another label after normalization"],
                )
            )

    return TaxonSafetyReport(
        policy=policy,
        unsafe_taxa=sorted(unsafe_taxa, key=lambda item: item.raw_label),
        collisions=collisions,
    )


def _normalize_node(node: TreeNode, *, policy: str, renames: list[TaxonRename]) -> TreeNode:
    normalized_name = node.name
    if node.name is not None:
        normalized_name = _normalize_label(node.name, policy=policy)
        if normalized_name != node.name:
            renames.append(TaxonRename(raw_label=node.name, normalized_label=normalized_name))
    return TreeNode(
        name=normalized_name,
        branch_length=node.branch_length,
        children=[_normalize_node(child, policy=policy, renames=renames) for child in node.children],
    )


def normalize_tree_taxa(tree: PhyloTree, *, policy: str) -> tuple[PhyloTree, TaxonNormalizationReport]:
    """Apply an explicit taxon normalization policy to a tree."""
    renames: list[TaxonRename] = []
    normalized_tree = PhyloTree(
        root=_normalize_node(tree.root, policy=policy, renames=renames),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    return normalized_tree, TaxonNormalizationReport(policy=policy, renamed_taxa=renames)


def write_taxon_mapping(path: Path, renames: list[TaxonRename]) -> Path:
    """Write a TSV mapping file for explicit taxon renames."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["raw_label\tnormalized_label"]
    lines.extend(f"{rename.raw_label}\t{rename.normalized_label}" for rename in renames)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
