from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .identity import inspect_tree_taxon_identity
from .models import (
    DuplicateBiologicalIdentityCandidate,
    DuplicateBiologicalIdentityReport,
    TaxonMappingConflictReport,
    TaxonMappingConflictRow,
)
from .normalization import inspect_tree_taxa_safety
from .synonyms import (
    _canonical_taxon_key,
    audit_tree_taxon_synonyms,
    export_tree_accepted_names,
)


def detect_duplicate_biological_identities(
    tree: PhyloTree,
    *,
    synonym_table_path: Path | None = None,
) -> DuplicateBiologicalIdentityReport:
    """Identify raw labels likely referring to the same biological taxon or sample."""
    candidates: set[tuple[str, str, str, str | None]] = set()
    identity = inspect_tree_taxon_identity(tree)
    for pair in identity.whitespace_variants:
        candidates.add((pair.left_label, pair.right_label, "whitespace_variant", None))
    for pair in identity.underscore_space_collisions:
        candidates.add(
            (pair.left_label, pair.right_label, "underscore_space_collision", None)
        )
    for pair in identity.case_collisions:
        candidates.add((pair.left_label, pair.right_label, "case_collision", None))
    for pair in identity.suspicious_near_duplicates:
        candidates.add((pair.left_label, pair.right_label, "near_duplicate", None))

    warnings: list[str] = []
    if synonym_table_path is not None:
        accepted = export_tree_accepted_names(tree, synonym_table_path)
        by_accepted: dict[str, list[str]] = defaultdict(list)
        for row in accepted.rows:
            if row.status == "ambiguous":
                continue
            by_accepted[_canonical_taxon_key(row.accepted_label)].append(row.raw_label)
        for accepted_key, raw_labels in sorted(by_accepted.items()):
            unique_raw = sorted(set(raw_labels))
            if len(unique_raw) < 2:
                continue
            accepted_label = next(
                row.accepted_label
                for row in accepted.rows
                if _canonical_taxon_key(row.accepted_label) == accepted_key
            )
            for index, left in enumerate(unique_raw):
                for right in unique_raw[index + 1 :]:
                    candidates.add(
                        (left, right, "shared_accepted_name", accepted_label)
                    )
        warnings.extend(accepted.warnings)

    ordered = [
        DuplicateBiologicalIdentityCandidate(
            left_label=left,
            right_label=right,
            evidence=evidence,
            accepted_label=accepted_label,
        )
        for left, right, evidence, accepted_label in sorted(candidates)
    ]
    if ordered:
        warnings.append(
            "tree contains labels that may represent duplicate biological identities"
        )
    return DuplicateBiologicalIdentityReport(
        candidates=ordered, warnings=sorted(dict.fromkeys(warnings))
    )


def build_taxon_mapping_conflict_report(
    tree: PhyloTree,
    *,
    synonym_table_path: Path | None = None,
) -> TaxonMappingConflictReport:
    """Collect explicit mapping conflicts across normalization, synonyms, and identity heuristics."""
    rows: list[TaxonMappingConflictRow] = []
    safety = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    for collision in safety.collisions:
        rows.append(
            TaxonMappingConflictRow(
                conflict_type="normalization_collision",
                raw_labels=collision.raw_labels,
                candidate_labels=[collision.normalized_label],
                detail="multiple raw labels collapse to the same normalized downstream-safe label",
            )
        )

    identity = inspect_tree_taxon_identity(tree)
    for pair in identity.suspicious_near_duplicates:
        rows.append(
            TaxonMappingConflictRow(
                conflict_type="near_duplicate_identity",
                raw_labels=[pair.left_label, pair.right_label],
                candidate_labels=[],
                detail="labels are close enough to suggest one biological identity may be duplicated",
            )
        )

    warnings: list[str] = []
    if synonym_table_path is not None:
        synonym_audit = audit_tree_taxon_synonyms(tree, synonym_table_path)
        for mapping in synonym_audit.ambiguous_mappings:
            rows.append(
                TaxonMappingConflictRow(
                    conflict_type="ambiguous_synonym",
                    raw_labels=[mapping.raw_label],
                    candidate_labels=mapping.accepted_labels,
                    detail="one raw label maps to multiple accepted taxa in the configured synonym table",
                )
            )
        accepted = export_tree_accepted_names(tree, synonym_table_path)
        by_accepted: dict[str, list[str]] = defaultdict(list)
        for row in accepted.rows:
            if row.status == "ambiguous":
                continue
            by_accepted[_canonical_taxon_key(row.accepted_label)].append(row.raw_label)
        for accepted_key, raw_labels in sorted(by_accepted.items()):
            unique_raw = sorted(set(raw_labels))
            if len(unique_raw) < 2:
                continue
            accepted_label = next(
                row.accepted_label
                for row in accepted.rows
                if _canonical_taxon_key(row.accepted_label) == accepted_key
            )
            rows.append(
                TaxonMappingConflictRow(
                    conflict_type="accepted_name_collision",
                    raw_labels=unique_raw,
                    candidate_labels=[accepted_label],
                    detail="multiple tree labels resolve to the same accepted name and may require deduplication",
                )
            )
        warnings.extend(synonym_audit.warnings)

    if rows:
        warnings.append(
            "one or more taxon mapping conflicts require reviewer attention"
        )
    return TaxonMappingConflictReport(
        rows=rows, warnings=sorted(dict.fromkeys(warnings))
    )
