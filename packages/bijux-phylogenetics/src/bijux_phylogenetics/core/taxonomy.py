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
    )
    return normalized_tree, TaxonNormalizationReport(policy=policy, renamed_taxa=renames)


def write_taxon_mapping(path: Path, renames: list[TaxonRename]) -> Path:
    """Write a TSV mapping file for explicit taxon renames."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["raw_label\tnormalized_label"]
    lines.extend(f"{rename.raw_label}\t{rename.normalized_label}" for rename in renames)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
