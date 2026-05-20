from __future__ import annotations

from collections import Counter, defaultdict
import csv
from pathlib import Path
import re

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode, normalize_taxon_key
from bijux_phylogenetics.runtime.errors import MetadataJoinError
from .identity import (
    _identity_key,
    _levenshtein_distance,
    _sorted_pairs,
    _space_underscore_key,
    _suspicious_near_duplicate,
    inspect_tree_taxon_identity,
)
from .models import (
    AcceptedNameExport,
    AcceptedNameRow,
    AmbiguousSynonymMapping,
    DuplicateBiologicalIdentityCandidate,
    DuplicateBiologicalIdentityReport,
    SynonymCandidate,
    TaxonAuditReport,
    TaxonIdentityAudit,
    TaxonLabelPair,
    TaxonMappingConflictReport,
    TaxonMappingConflictRow,
    TaxonNamespaceAssignment,
    TaxonNamespaceReport,
    TaxonNormalizationCollision,
    TaxonNormalizationReport,
    TaxonRankAssignment,
    TaxonRankConsistencyReport,
    TaxonRename,
    TaxonSafetyReport,
    TaxonSynonymAudit,
    TaxonSynonymResolutionReport,
    TaxonSynonymResolutionRow,
    TaxonSynonymRow,
    TaxonomyReference,
    UnsafeTaxonName,
)
from .normalization import (
    _normalize_label,
    _normalize_node,
    _unsafe_reasons,
    inspect_tree_taxa_safety,
    normalize_tree_taxa,
    write_taxon_mapping,
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


def _detect_delimiter(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ",", "csv"
    if suffix == ".tsv":
        return "\t", "tsv"
    header_line = (
        path.read_text(encoding="utf-8").splitlines()[0] if path.exists() else ""
    )
    if "\t" in header_line:
        return "\t", "tsv"
    return ",", "csv"


def _canonical_taxon_key(label: str) -> str:
    return normalize_taxon_key(label).casefold()


def _resolve_column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    normalized = {column.strip().lower(): column for column in columns}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    return None


def load_taxon_synonym_rows(path: Path) -> list[TaxonSynonymRow]:
    """Load a synonym table without requiring unique raw labels."""
    if not path.exists():
        raise FileNotFoundError(f"synonym table not found: {path}")
    delimiter, _ = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise MetadataJoinError(f"synonym table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        raw_column = _resolve_column(
            columns, ("raw_label", "synonym", "label", "taxon")
        )
        accepted_column = _resolve_column(
            columns, ("accepted_label", "accepted_name", "accepted", "canonical_label")
        )
        provenance_column = _resolve_column(
            columns, ("provenance", "source", "authority_source", "note")
        )
        authority_column = _resolve_column(
            columns, ("authority", "namespace", "catalog")
        )
        if raw_column is None or accepted_column is None:
            raise MetadataJoinError(
                "synonym table must contain raw_label/synonym and accepted_label/accepted_name columns"
            )
        rows: list[TaxonSynonymRow] = []
        for row_index, row in enumerate(reader, start=2):
            raw_label = str(row.get(raw_column, "")).strip()
            accepted_label = str(row.get(accepted_column, "")).strip()
            if not raw_label or not accepted_label:
                raise MetadataJoinError(
                    f"row {row_index} in {path} must contain non-empty synonym and accepted-name values"
                )
            rows.append(
                TaxonSynonymRow(
                    raw_label=raw_label,
                    accepted_label=accepted_label,
                    provenance=""
                    if provenance_column is None
                    else str(row.get(provenance_column, "")).strip(),
                    authority=None
                    if authority_column is None
                    else (str(row.get(authority_column, "")).strip() or None),
                )
            )
    return rows


def _group_synonym_rows(
    rows: list[TaxonSynonymRow],
) -> tuple[dict[str, list[TaxonSynonymRow]], list[AmbiguousSynonymMapping]]:
    grouped: dict[str, list[TaxonSynonymRow]] = defaultdict(list)
    ambiguous: list[AmbiguousSynonymMapping] = []
    for row in rows:
        grouped[_canonical_taxon_key(row.raw_label)].append(row)
    for _key, grouped_rows in sorted(grouped.items()):
        accepted_labels = sorted({row.accepted_label for row in grouped_rows})
        if len(accepted_labels) > 1:
            ambiguous.append(
                AmbiguousSynonymMapping(
                    raw_label=min((row.raw_label for row in grouped_rows), key=len),
                    accepted_labels=accepted_labels,
                    provenances=sorted(
                        {row.provenance for row in grouped_rows if row.provenance}
                    ),
                )
            )
    return grouped, ambiguous


def audit_tree_taxon_synonyms(
    tree: PhyloTree, synonym_table_path: Path
) -> TaxonSynonymAudit:
    """Audit tree tips against a configurable synonym table."""
    synonym_rows = load_taxon_synonym_rows(synonym_table_path)
    grouped_rows, ambiguous_mappings = _group_synonym_rows(synonym_rows)
    ambiguous_keys = {
        _canonical_taxon_key(entry.raw_label) for entry in ambiguous_mappings
    }
    candidates: list[SynonymCandidate] = []
    for label in sorted(tree.tip_names):
        rows = grouped_rows.get(_canonical_taxon_key(label), [])
        if not rows or _canonical_taxon_key(label) in ambiguous_keys:
            continue
        row = rows[0]
        if _canonical_taxon_key(row.raw_label) == _canonical_taxon_key(
            row.accepted_label
        ):
            continue
        candidates.append(
            SynonymCandidate(
                raw_label=label,
                accepted_label=row.accepted_label,
                provenance=row.provenance,
                authority=row.authority,
            )
        )
    warnings: list[str] = []
    if ambiguous_mappings:
        warnings.append(
            "synonym table contains one or more raw labels that map to multiple accepted taxa"
        )
    if candidates:
        warnings.append(
            "tree contains one or more labels with configured synonym candidates"
        )
    return TaxonSynonymAudit(
        synonym_table_path=synonym_table_path,
        candidates=candidates,
        ambiguous_mappings=ambiguous_mappings,
        warnings=warnings,
    )


def _resolve_node_synonyms(
    node: TreeNode,
    *,
    grouped_rows: dict[str, list[TaxonSynonymRow]],
    ambiguous_keys: set[str],
    renames: list[TaxonSynonymResolutionRow],
) -> TreeNode:
    resolved_name = node.name
    if node.name is not None:
        key = _canonical_taxon_key(node.name)
        if key not in ambiguous_keys and key in grouped_rows:
            row = grouped_rows[key][0]
            resolved_name = row.accepted_label
            if _canonical_taxon_key(node.name) != _canonical_taxon_key(resolved_name):
                renames.append(
                    TaxonSynonymResolutionRow(
                        raw_label=node.name,
                        resolved_label=resolved_name,
                        accepted_label=row.accepted_label,
                        provenance=row.provenance,
                        authority=row.authority,
                    )
                )
    return TreeNode(
        name=resolved_name,
        branch_length=node.branch_length,
        children=[
            _resolve_node_synonyms(
                child,
                grouped_rows=grouped_rows,
                ambiguous_keys=ambiguous_keys,
                renames=renames,
            )
            for child in node.children
        ],
    )


def resolve_tree_taxon_synonyms(
    tree: PhyloTree,
    *,
    synonym_table_path: Path,
    resolution_policy: str = "reject-ambiguous",
) -> tuple[PhyloTree, TaxonSynonymResolutionReport]:
    """Resolve tree tip synonyms to accepted labels with explicit provenance."""
    if resolution_policy != "reject-ambiguous":
        raise ValueError(f"unsupported synonym-resolution policy: {resolution_policy}")
    synonym_rows = load_taxon_synonym_rows(synonym_table_path)
    grouped_rows, ambiguous_mappings = _group_synonym_rows(synonym_rows)
    ambiguous_keys = {
        _canonical_taxon_key(entry.raw_label) for entry in ambiguous_mappings
    }
    renames: list[TaxonSynonymResolutionRow] = []
    resolved_tree = PhyloTree(
        root=_resolve_node_synonyms(
            tree.root,
            grouped_rows=grouped_rows,
            ambiguous_keys=ambiguous_keys,
            renames=renames,
        ),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    grouped_resolved: dict[str, list[str]] = defaultdict(list)
    for label in resolved_tree.tip_names:
        grouped_resolved[_canonical_taxon_key(label)].append(label)
    duplicate_resolved_labels = [
        TaxonNormalizationCollision(
            normalized_label=labels[0],
            raw_labels=sorted(labels),
        )
        for labels in grouped_resolved.values()
        if len(labels) > 1
    ]
    warnings: list[str] = []
    if ambiguous_mappings:
        warnings.append(
            "one or more synonym rows were left unresolved because they map to multiple accepted taxa"
        )
    if duplicate_resolved_labels:
        warnings.append(
            "synonym resolution collapses two or more tree tips to the same accepted label"
        )
    renamed_labels = {row.raw_label for row in renames}
    return resolved_tree, TaxonSynonymResolutionReport(
        synonym_table_path=synonym_table_path,
        resolution_policy=resolution_policy,
        renamed_taxa=sorted(renames, key=lambda row: row.raw_label),
        unchanged_taxa=sorted(
            label for label in tree.tip_names if label not in renamed_labels
        ),
        ambiguous_mappings=ambiguous_mappings,
        duplicate_resolved_labels=duplicate_resolved_labels,
        warnings=warnings,
    )


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


def export_tree_accepted_names(
    tree: PhyloTree, synonym_table_path: Path
) -> AcceptedNameExport:
    """Export raw tree labels to accepted names with explicit resolution status."""
    synonym_rows = load_taxon_synonym_rows(synonym_table_path)
    grouped_rows, ambiguous_mappings = _group_synonym_rows(synonym_rows)
    ambiguous_by_key = {
        _canonical_taxon_key(entry.raw_label): entry for entry in ambiguous_mappings
    }
    rows: list[AcceptedNameRow] = []
    for label in sorted(tree.tip_names):
        key = _canonical_taxon_key(label)
        ambiguous = ambiguous_by_key.get(key)
        if ambiguous is not None:
            rows.append(
                AcceptedNameRow(
                    raw_label=label,
                    accepted_label=" | ".join(ambiguous.accepted_labels),
                    status="ambiguous",
                    provenance=" | ".join(ambiguous.provenances),
                    authority=None,
                )
            )
            continue
        candidates = grouped_rows.get(key, [])
        if not candidates:
            rows.append(
                AcceptedNameRow(
                    raw_label=label,
                    accepted_label=label,
                    status="unchanged",
                    provenance="",
                    authority=None,
                )
            )
            continue
        candidate = candidates[0]
        rows.append(
            AcceptedNameRow(
                raw_label=label,
                accepted_label=candidate.accepted_label,
                status=(
                    "resolved"
                    if _canonical_taxon_key(label)
                    != _canonical_taxon_key(candidate.accepted_label)
                    else "already_accepted"
                ),
                provenance=candidate.provenance,
                authority=candidate.authority,
            )
        )
    warnings: list[str] = []
    if ambiguous_mappings:
        warnings.append(
            "accepted-name export includes ambiguous synonym rows that require manual review"
        )
    return AcceptedNameExport(
        synonym_table_path=synonym_table_path, rows=rows, warnings=warnings
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


def write_synonym_resolution_mapping(
    path: Path, report: TaxonSynonymResolutionReport
) -> Path:
    """Write a TSV mapping file for synonym resolution with provenance."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["raw_label\tresolved_label\taccepted_label\tprovenance\tauthority"]
    lines.extend(
        "\t".join(
            [
                row.raw_label,
                row.resolved_label,
                row.accepted_label,
                row.provenance,
                "" if row.authority is None else row.authority,
            ]
        )
        for row in report.renamed_taxa
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_accepted_name_mapping(path: Path, report: AcceptedNameExport) -> Path:
    """Write the accepted-name export as a TSV mapping file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["raw_label\taccepted_label\tstatus\tprovenance\tauthority"]
    lines.extend(
        "\t".join(
            [
                row.raw_label,
                row.accepted_label,
                row.status,
                row.provenance,
                "" if row.authority is None else row.authority,
            ]
        )
        for row in report.rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
