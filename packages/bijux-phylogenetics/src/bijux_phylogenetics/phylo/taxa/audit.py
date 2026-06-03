from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .classification import (
    inspect_tree_taxon_namespaces,
    inspect_tree_taxon_rank_consistency,
)
from .conflicts import (
    build_taxon_mapping_conflict_report,
    detect_duplicate_biological_identities,
)
from .identity import inspect_tree_taxon_identity
from .models import TaxonAuditReport
from .normalization import inspect_tree_taxa_safety
from .synonyms import (
    audit_tree_taxon_synonyms,
    export_tree_accepted_names,
)


def build_taxon_audit_report(
    tree: PhyloTree,
    *,
    synonym_table_path: Path | None = None,
) -> TaxonAuditReport:
    """Build a reviewer-readable audit of taxon naming, safety, and mapping trust."""
    identity = inspect_tree_taxon_identity(tree)
    safety = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    namespace_report = inspect_tree_taxon_namespaces(tree)
    rank_consistency = inspect_tree_taxon_rank_consistency(tree)
    synonym_audit = (
        None
        if synonym_table_path is None
        else audit_tree_taxon_synonyms(tree, synonym_table_path)
    )
    duplicate_identities = detect_duplicate_biological_identities(
        tree, synonym_table_path=synonym_table_path
    )
    mapping_conflicts = build_taxon_mapping_conflict_report(
        tree, synonym_table_path=synonym_table_path
    )
    accepted_name_export = (
        None
        if synonym_table_path is None
        else export_tree_accepted_names(tree, synonym_table_path)
    )

    summary: list[str] = []
    warnings = (
        list(namespace_report.warnings)
        + list(rank_consistency.warnings)
        + list(duplicate_identities.warnings)
        + list(mapping_conflicts.warnings)
    )
    identity_variant_pairs = (
        len(identity.spelling_variants)
        + len(identity.whitespace_variants)
        + len(identity.underscore_space_collisions)
        + len(identity.case_collisions)
    )
    if identity_variant_pairs:
        warnings.append(
            "one or more taxon labels differ only by spacing, punctuation, or case and should be reviewed before normalization"
        )
        summary.append(
            f"{identity_variant_pairs} spelling or formatting collision pairs were detected across taxon labels"
        )
    if identity.suspicious_near_duplicates:
        warnings.append(
            "one or more near-duplicate taxon labels may represent the same biological identity"
        )
        summary.append(
            f"{len(identity.suspicious_near_duplicates)} near-duplicate taxon label pairs need manual review"
        )
    if safety.unsafe_taxa:
        warnings.append(
            "one or more taxon labels are unsafe for downstream external engines"
        )
        summary.append(
            f"{len(safety.unsafe_taxa)} labels need downstream-safe normalization or quoting"
        )
    if namespace_report.mixed_namespaces:
        summary.append(
            "tree mixes explicit taxon namespaces and needs a crosswalk before automation"
        )
    if rank_consistency.mixed_ranks:
        summary.append(
            "tree mixes biological naming levels and should not be treated as one uniform taxon namespace"
        )
    if mapping_conflicts.rows:
        summary.append(
            f"{len(mapping_conflicts.rows)} mapping conflicts need manual review"
        )
    if duplicate_identities.candidates:
        summary.append(
            f"{len(duplicate_identities.candidates)} duplicate-identity candidates were detected"
        )
    if synonym_audit is not None and synonym_audit.candidates:
        summary.append(
            f"{len(synonym_audit.candidates)} labels have accepted-name synonym candidates"
        )

    status = "ready"
    if (
        namespace_report.mixed_namespaces
        or rank_consistency.mixed_ranks
        or mapping_conflicts.rows
    ):
        status = "needs_review"
    if safety.collisions:
        status = "blocked"

    return TaxonAuditReport(
        tree_tip_count=tree.tip_count,
        identity_audit=identity,
        safety_report=safety,
        namespace_report=namespace_report,
        rank_consistency=rank_consistency,
        synonym_audit=synonym_audit,
        duplicate_identities=duplicate_identities,
        mapping_conflicts=mapping_conflicts,
        accepted_name_export=accepted_name_export,
        status=status,
        summary=summary,
        warnings=sorted(dict.fromkeys(warnings)),
    )
