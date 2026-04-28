from __future__ import annotations

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


def _normalize_label(label: str, *, policy: str) -> str:
    if policy == "spaces-to-underscores":
        return "_".join(label.split())
    raise ValueError(f"unsupported taxon normalization policy: {policy}")


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
