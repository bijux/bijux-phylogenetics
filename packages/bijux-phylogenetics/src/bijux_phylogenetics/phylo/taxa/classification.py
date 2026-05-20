from __future__ import annotations

from collections import Counter
import re

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    TaxonNamespaceAssignment,
    TaxonNamespaceReport,
    TaxonRankAssignment,
    TaxonRankConsistencyReport,
)

_SAMPLE_TOKENS = ("sample", "specimen", "individual", "indiv", "voucher")
_ISOLATE_TOKENS = ("isolate", "strain", "clone", "cultivar", "lineage")
_ACCESSION_PATTERNS = (
    re.compile(r"^[A-Z]{1,4}_?\d{4,}(?:\.\d+)?$"),
    re.compile(r"^[A-Z]{2}\d{6,}(?:\.\d+)?$"),
)
_SPECIES_NAME_PATTERN = re.compile(r"^[A-Z][a-z]+(?:[ _][a-z][a-z-]+){1,2}$")
_GENUS_NAME_PATTERN = re.compile(r"^[A-Z][a-z-]+$")
_POPULATION_TOKENS = ("population", "pop", "subpopulation", "lineage", "clade")


def infer_taxon_namespace(label: str) -> TaxonNamespaceAssignment:
    """Classify one taxon label into a pragmatic namespace heuristic."""
    normalized = label.strip()
    lowered = normalized.lower()
    for pattern in _ACCESSION_PATTERNS:
        if pattern.fullmatch(normalized):
            return TaxonNamespaceAssignment(
                label=label,
                namespace="accession_id",
                evidence="label matches a compact accession-like identifier pattern",
            )
    if _SPECIES_NAME_PATTERN.fullmatch(normalized):
        return TaxonNamespaceAssignment(
            label=label,
            namespace="species_name",
            evidence="label matches a Latin-like binomial or trinomial pattern",
        )
    if any(token in lowered for token in _ISOLATE_TOKENS):
        return TaxonNamespaceAssignment(
            label=label,
            namespace="isolate_id",
            evidence="label contains isolate or strain-like tokens",
        )
    if any(token in lowered for token in _SAMPLE_TOKENS) or re.fullmatch(
        r"[A-Za-z]+[-_]\d{2,}", normalized
    ):
        return TaxonNamespaceAssignment(
            label=label,
            namespace="sample_id",
            evidence="label contains sample-style tokens or a sample-like alphanumeric pattern",
        )
    return TaxonNamespaceAssignment(
        label=label,
        namespace="user_defined_label",
        evidence="label does not match accession, species, sample, or isolate heuristics",
    )


def inspect_tree_taxon_namespaces(tree: PhyloTree) -> TaxonNamespaceReport:
    """Audit tree tips for namespace consistency and mixed identifier styles."""
    assignments = sorted(
        (infer_taxon_namespace(label) for label in tree.tip_names),
        key=lambda row: row.label,
    )
    namespace_counts = dict(
        sorted(Counter(row.namespace for row in assignments).items())
    )
    explicit_namespaces = {
        namespace for namespace in namespace_counts if namespace != "user_defined_label"
    }
    mixed_namespaces = len(explicit_namespaces) > 1
    warnings: list[str] = []
    if mixed_namespaces:
        warnings.append(
            "tree mixes multiple explicit taxon namespaces without an explicit mapping layer"
        )
    if (
        explicit_namespaces
        and "user_defined_label" in namespace_counts
        and len(namespace_counts) > 1
    ):
        warnings.append(
            "tree mixes structured namespaces with arbitrary user-defined labels"
        )
    dominant_namespace = (
        None
        if not namespace_counts
        else max(namespace_counts.items(), key=lambda item: (item[1], item[0]))[0]
    )
    return TaxonNamespaceReport(
        assignments=assignments,
        namespace_counts=namespace_counts,
        dominant_namespace=dominant_namespace,
        mixed_namespaces=mixed_namespaces,
        warnings=warnings,
    )


def infer_taxon_rank(label: str) -> TaxonRankAssignment:
    """Infer a pragmatic biological naming level for one taxon label."""
    namespace = infer_taxon_namespace(label)
    normalized = label.strip()
    lowered = normalized.lower()
    if namespace.namespace == "accession_id":
        return TaxonRankAssignment(
            label=label, inferred_rank="accession", evidence=namespace.evidence
        )
    if namespace.namespace in {"sample_id", "isolate_id"}:
        return TaxonRankAssignment(
            label=label, inferred_rank="sample", evidence=namespace.evidence
        )
    if _SPECIES_NAME_PATTERN.fullmatch(normalized):
        return TaxonRankAssignment(
            label=label,
            inferred_rank="species",
            evidence="label matches a Latin-like species or subspecies name",
        )
    if _GENUS_NAME_PATTERN.fullmatch(normalized):
        return TaxonRankAssignment(
            label=label,
            inferred_rank="genus",
            evidence="label matches a single capitalized genus-like token",
        )
    if any(token in lowered for token in _POPULATION_TOKENS):
        return TaxonRankAssignment(
            label=label,
            inferred_rank="population",
            evidence="label contains population-like grouping tokens",
        )
    return TaxonRankAssignment(
        label=label,
        inferred_rank="unknown",
        evidence="label does not match accession, sample, genus, species, or population heuristics",
    )


def inspect_tree_taxon_rank_consistency(tree: PhyloTree) -> TaxonRankConsistencyReport:
    """Flag datasets that mix species, genus, accession, sample, or population-level labels."""
    assignments = sorted(
        (infer_taxon_rank(label) for label in tree.tip_names), key=lambda row: row.label
    )
    rank_counts = dict(
        sorted(Counter(row.inferred_rank for row in assignments).items())
    )
    explicit_ranks = {rank for rank in rank_counts if rank != "unknown"}
    mixed_ranks = len(explicit_ranks) > 1
    warnings: list[str] = []
    if mixed_ranks:
        warnings.append(
            "tree mixes multiple biological naming levels without an explicit crosswalk"
        )
    if explicit_ranks and "unknown" in rank_counts and len(rank_counts) > 1:
        warnings.append(
            "tree mixes interpretable biological ranks with labels of unknown rank semantics"
        )
    dominant_rank = (
        None
        if not rank_counts
        else max(rank_counts.items(), key=lambda item: (item[1], item[0]))[0]
    )
    return TaxonRankConsistencyReport(
        assignments=assignments,
        rank_counts=rank_counts,
        dominant_rank=dominant_rank,
        mixed_ranks=mixed_ranks,
        warnings=warnings,
    )
