from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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
    original_tip_count: int
    normalized_tip_count: int
    unchanged_taxa: list[str]
    topology_preserved: bool
    branch_lengths_preserved: bool


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


@dataclass(frozen=True, slots=True)
class TaxonSynonymRow:
    """One configured synonym mapping row."""

    raw_label: str
    accepted_label: str
    provenance: str
    authority: str | None


@dataclass(frozen=True, slots=True)
class SynonymCandidate:
    """One tree label likely referring to an accepted taxon through a synonym table."""

    raw_label: str
    accepted_label: str
    provenance: str
    authority: str | None


@dataclass(frozen=True, slots=True)
class AmbiguousSynonymMapping:
    """One raw label that maps to multiple accepted taxa."""

    raw_label: str
    accepted_labels: list[str]
    provenances: list[str]


@dataclass(slots=True)
class TaxonSynonymAudit:
    """Audit likely synonym candidates for one tree."""

    synonym_table_path: Path
    candidates: list[SynonymCandidate]
    ambiguous_mappings: list[AmbiguousSynonymMapping]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TaxonSynonymResolutionRow:
    """One reversible synonym normalization mapping."""

    raw_label: str
    resolved_label: str
    accepted_label: str
    provenance: str
    authority: str | None


@dataclass(slots=True)
class TaxonSynonymResolutionReport:
    """Controlled synonym-resolution report for one tree."""

    synonym_table_path: Path
    resolution_policy: str
    renamed_taxa: list[TaxonSynonymResolutionRow]
    unchanged_taxa: list[str]
    ambiguous_mappings: list[AmbiguousSynonymMapping]
    duplicate_resolved_labels: list[TaxonNormalizationCollision]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TaxonNamespaceAssignment:
    """Namespace classification for one tree tip label."""

    label: str
    namespace: str
    evidence: str


@dataclass(slots=True)
class TaxonNamespaceReport:
    """Explicit namespace audit across the tip labels of one tree."""

    assignments: list[TaxonNamespaceAssignment]
    namespace_counts: dict[str, int]
    dominant_namespace: str | None
    mixed_namespaces: bool
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TaxonRankAssignment:
    """Rank-like interpretation for one taxon label."""

    label: str
    inferred_rank: str
    evidence: str


@dataclass(slots=True)
class TaxonRankConsistencyReport:
    """Audit whether a tree mixes biologically incompatible naming levels."""

    assignments: list[TaxonRankAssignment]
    rank_counts: dict[str, int]
    dominant_rank: str | None
    mixed_ranks: bool
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class DuplicateBiologicalIdentityCandidate:
    """Two labels likely representing the same biological identity."""

    left_label: str
    right_label: str
    evidence: str
    accepted_label: str | None


@dataclass(slots=True)
class DuplicateBiologicalIdentityReport:
    """Likely duplicate biological identities across raw tree labels."""

    candidates: list[DuplicateBiologicalIdentityCandidate]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AcceptedNameRow:
    """One accepted-name export row for a raw tree label."""

    raw_label: str
    accepted_label: str
    status: str
    provenance: str
    authority: str | None


@dataclass(slots=True)
class AcceptedNameExport:
    """Accepted-name export for all tree labels under a synonym policy."""

    synonym_table_path: Path
    rows: list[AcceptedNameRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TaxonMappingConflictRow:
    """One taxon mapping conflict a reviewer should inspect."""

    conflict_type: str
    raw_labels: list[str]
    candidate_labels: list[str]
    detail: str


@dataclass(slots=True)
class TaxonMappingConflictReport:
    """Explicit conflict table across synonym, identity, and normalization mapping."""

    rows: list[TaxonMappingConflictRow]
    warnings: list[str]


@dataclass(slots=True)
class TaxonAuditReport:
    """Reviewer-readable taxon audit across identity, safety, synonyms, and mapping."""

    tree_tip_count: int
    identity_audit: TaxonIdentityAudit
    safety_report: TaxonSafetyReport
    namespace_report: TaxonNamespaceReport
    rank_consistency: TaxonRankConsistencyReport
    synonym_audit: TaxonSynonymAudit | None
    duplicate_identities: DuplicateBiologicalIdentityReport
    mapping_conflicts: TaxonMappingConflictReport
    accepted_name_export: AcceptedNameExport | None
    status: str
    summary: list[str]
    warnings: list[str]
