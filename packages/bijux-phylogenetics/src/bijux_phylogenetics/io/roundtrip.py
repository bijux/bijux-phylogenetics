from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from Bio import Phylo

from bijux_phylogenetics.compare.topology import (
    _compare_branch_lengths_for_trees,
    _compare_tree_objects,
)
from bijux_phylogenetics.io.biopython import tree_to_biophylo
from bijux_phylogenetics.io.trees import TreeFormat, detect_tree_format, load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


@dataclass(slots=True)
class SupportNormalizationAuditRow:
    node: str
    raw_label: str
    raw_value: float
    inferred_scale: str
    normalized_probability: float | None
    confidence_of_inference: str


@dataclass(slots=True)
class SemanticLossObservation:
    feature: str
    severity: str
    message: str


@dataclass(slots=True)
class TreeRoundtripValidationReport:
    source_path: Path
    source_format: str
    target_format: str
    preserved_taxa: bool
    preserved_topology: bool
    preserved_branch_lengths: bool
    preserved_support_labels: bool
    rooted_topology_preserved: bool
    semantic_loss: list[SemanticLossObservation]
    support_audit: list[SupportNormalizationAuditRow]


def _parse_internal_label_numeric(label: str) -> float | None:
    try:
        return float(label)
    except ValueError:
        return None


def _support_rows(tree: PhyloTree) -> list[SupportNormalizationAuditRow]:
    rows: list[SupportNormalizationAuditRow] = []
    mixed_fraction = False
    mixed_percent = False
    staged: list[tuple[str, str, float]] = []

    def descendant_taxa(node: TreeNode) -> list[str]:
        if node.is_leaf():
            return [node.name] if node.name is not None else []
        taxa: list[str] = []
        for child in node.children:
            taxa.extend(descendant_taxa(child))
        return sorted(taxa)

    def node_id(node: TreeNode) -> str:
        return node.node_id or (
            "|".join(descendant_taxa(node))
            if descendant_taxa(node)
            else (node.name or "<unnamed>")
        )

    for node in tree.iter_nodes():
        if node.is_leaf() or node.name is None or not node.name.strip():
            continue
        numeric = _parse_internal_label_numeric(node.name.strip())
        if numeric is None:
            continue
        staged.append((node_id(node), node.name, numeric))
        if 0.0 <= numeric <= 1.0:
            mixed_fraction = True
        elif 1.0 < numeric <= 100.0:
            mixed_percent = True

    for signature, label, numeric in staged:
        inferred_scale = "out-of-range"
        normalized = None
        confidence = "low"
        if 0.0 <= numeric <= 1.0:
            inferred_scale = "fraction"
            normalized = round(numeric, 15)
            confidence = "medium" if mixed_percent else "high"
            if numeric in {0.0, 1.0}:
                confidence = "medium" if not mixed_percent else "low"
        elif 1.0 < numeric <= 100.0:
            inferred_scale = "percentage"
            normalized = round(numeric / 100.0, 15)
            confidence = "medium" if mixed_fraction else "high"
        rows.append(
            SupportNormalizationAuditRow(
                node=signature,
                raw_label=label,
                raw_value=numeric,
                inferred_scale=inferred_scale,
                normalized_probability=normalized,
                confidence_of_inference=confidence,
            )
        )
    return sorted(rows, key=lambda row: (row.node, row.raw_label))


def _write_tree(path: Path, tree: PhyloTree, *, target_format: TreeFormat) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if target_format == "newick":
        from bijux_phylogenetics.io.newick import write_newick

        return write_newick(path, tree)
    biophylo_tree = tree_to_biophylo(tree)
    Phylo.write(biophylo_tree, path, target_format)
    return path


def _semantic_loss(
    path: Path, *, source_format: str, target_format: str
) -> list[SemanticLossObservation]:
    text = path.read_text(encoding="utf-8")
    lowered = text.lower()
    observations: list[SemanticLossObservation] = []

    if source_format == "nexus":
        if (
            "begin taxa" in lowered
            or "begin data" in lowered
            or "begin characters" in lowered
        ):
            observations.append(
                SemanticLossObservation(
                    feature="nexus-metadata-blocks",
                    severity="warning",
                    message="nexus taxa or character metadata blocks are not represented in the local tree model and cannot be proven preserved",
                )
            )
        if "translate" in lowered:
            observations.append(
                SemanticLossObservation(
                    feature="nexus-translate",
                    severity="info",
                    message="nexus translate token mappings are normalized to resolved labels and the original token table is not preserved",
                )
            )
    if source_format == "phyloxml":
        feature_map = {
            "<property": "phyloxml properties are not preserved in the local tree model",
            "<taxonomy": "phyloxml taxonomy annotations are not preserved in the local tree model",
            "<sequence": "phyloxml sequence annotations are not preserved in the local tree model",
            "<distribution": "phyloxml distribution annotations are not preserved in the local tree model",
            "<events": "phyloxml event annotations are not preserved in the local tree model",
            "<date": "phyloxml date annotations are not preserved in the local tree model",
        }
        for marker, message in feature_map.items():
            if marker in lowered:
                observations.append(
                    SemanticLossObservation(
                        feature=marker.removeprefix("<"),
                        severity="warning",
                        message=message,
                    )
                )
    if target_format == "newick" and source_format in {"nexus", "phyloxml"}:
        observations.append(
            SemanticLossObservation(
                feature="target-newick-structure",
                severity="warning",
                message="newick output cannot carry structured annotations beyond labels and branch lengths",
            )
        )
    return observations


def validate_tree_roundtrip(
    source_path: Path,
    *,
    target_format: TreeFormat,
    out_path: Path | None = None,
) -> TreeRoundtripValidationReport:
    """Convert a tree into another format, reload it, and validate preserved structure."""
    source_format = detect_tree_format(source_path)
    original = load_tree(source_path, source_format=source_format)
    target_path = out_path
    if target_path is None:
        suffix = {"newick": ".nwk", "nexus": ".nex", "phyloxml": ".phyloxml"}[
            target_format
        ]
        temp_dir = Path(tempfile.mkdtemp(prefix="bijux-tree-roundtrip-"))
        target_path = temp_dir / f"roundtrip{suffix}"
    _write_tree(target_path, original, target_format=target_format)
    reloaded = load_tree(target_path, source_format=target_format)
    topology = _compare_tree_objects(original, reloaded)
    branch_lengths = _compare_branch_lengths_for_trees(
        source_path,
        target_path,
        original,
        reloaded,
        taxon_overlap_policy="require-identical",
    )
    original_support = _support_rows(original)
    reloaded_support = _support_rows(reloaded)
    return TreeRoundtripValidationReport(
        source_path=source_path,
        source_format=source_format,
        target_format=target_format,
        preserved_taxa=original.tip_names == reloaded.tip_names,
        preserved_topology=topology.same_unrooted_topology,
        preserved_branch_lengths=all(
            row.left_length == row.right_length for row in branch_lengths.shared_splits
        ),
        preserved_support_labels=[
            (row.node, row.raw_value, row.inferred_scale, row.normalized_probability)
            for row in original_support
        ]
        == [
            (row.node, row.raw_value, row.inferred_scale, row.normalized_probability)
            for row in reloaded_support
        ],
        rooted_topology_preserved=topology.topology_equal,
        semantic_loss=_semantic_loss(
            source_path, source_format=source_format, target_format=target_format
        ),
        support_audit=original_support,
    )
