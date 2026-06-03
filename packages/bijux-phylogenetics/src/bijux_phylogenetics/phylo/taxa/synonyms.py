from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import (
    PhyloTree,
    TreeNode,
    normalize_taxon_key,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError

from .models import (
    AcceptedNameExport,
    AcceptedNameRow,
    AmbiguousSynonymMapping,
    SynonymCandidate,
    TaxonNormalizationCollision,
    TaxonSynonymAudit,
    TaxonSynonymResolutionReport,
    TaxonSynonymResolutionRow,
    TaxonSynonymRow,
)


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
    for grouped_rows in grouped.values():
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
    return grouped, sorted(ambiguous, key=lambda row: row.raw_label)


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
