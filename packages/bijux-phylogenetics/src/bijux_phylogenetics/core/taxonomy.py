from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import re

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode, normalize_taxon_key
from bijux_phylogenetics.errors import MetadataJoinError


_SAMPLE_TOKENS = ("sample", "specimen", "individual", "indiv", "voucher")
_ISOLATE_TOKENS = ("isolate", "strain", "clone", "cultivar", "lineage")
_ACCESSION_PATTERNS = (
    re.compile(r"^[A-Z]{1,4}_?\d{4,}(?:\.\d+)?$"),
    re.compile(r"^[A-Z]{2}\d{6,}(?:\.\d+)?$"),
)
_SPECIES_NAME_PATTERN = re.compile(r"^[A-Z][a-z]+(?:[ _][a-z][a-z-]+){1,2}$")


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


def _detect_delimiter(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ",", "csv"
    if suffix == ".tsv":
        return "\t", "tsv"
    header_line = path.read_text(encoding="utf-8").splitlines()[0] if path.exists() else ""
    if "\t" in header_line:
        return "\t", "tsv"
    return ",", "csv"


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
        raw_column = _resolve_column(columns, ("raw_label", "synonym", "label", "taxon"))
        accepted_column = _resolve_column(columns, ("accepted_label", "accepted_name", "accepted", "canonical_label"))
        provenance_column = _resolve_column(columns, ("provenance", "source", "authority_source", "note"))
        authority_column = _resolve_column(columns, ("authority", "namespace", "catalog"))
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
                    provenance="" if provenance_column is None else str(row.get(provenance_column, "")).strip(),
                    authority=None if authority_column is None else (str(row.get(authority_column, "")).strip() or None),
                )
            )
    return rows


def _group_synonym_rows(rows: list[TaxonSynonymRow]) -> tuple[dict[str, list[TaxonSynonymRow]], list[AmbiguousSynonymMapping]]:
    grouped: dict[str, list[TaxonSynonymRow]] = defaultdict(list)
    ambiguous: list[AmbiguousSynonymMapping] = []
    for row in rows:
        grouped[_canonical_taxon_key(row.raw_label)].append(row)
    for key, grouped_rows in sorted(grouped.items()):
        accepted_labels = sorted({row.accepted_label for row in grouped_rows})
        if len(accepted_labels) > 1:
            ambiguous.append(
                AmbiguousSynonymMapping(
                    raw_label=min((row.raw_label for row in grouped_rows), key=len),
                    accepted_labels=accepted_labels,
                    provenances=sorted({row.provenance for row in grouped_rows if row.provenance}),
                )
            )
    return grouped, ambiguous


def audit_tree_taxon_synonyms(tree: PhyloTree, synonym_table_path: Path) -> TaxonSynonymAudit:
    """Audit tree tips against a configurable synonym table."""
    synonym_rows = load_taxon_synonym_rows(synonym_table_path)
    grouped_rows, ambiguous_mappings = _group_synonym_rows(synonym_rows)
    ambiguous_keys = {_canonical_taxon_key(entry.raw_label) for entry in ambiguous_mappings}
    candidates: list[SynonymCandidate] = []
    for label in sorted(tree.tip_names):
        rows = grouped_rows.get(_canonical_taxon_key(label), [])
        if not rows or _canonical_taxon_key(label) in ambiguous_keys:
            continue
        row = rows[0]
        if _canonical_taxon_key(row.raw_label) == _canonical_taxon_key(row.accepted_label):
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
        warnings.append("synonym table contains one or more raw labels that map to multiple accepted taxa")
    if candidates:
        warnings.append("tree contains one or more labels with configured synonym candidates")
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
    ambiguous_keys = {_canonical_taxon_key(entry.raw_label) for entry in ambiguous_mappings}
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
        warnings.append("one or more synonym rows were left unresolved because they map to multiple accepted taxa")
    if duplicate_resolved_labels:
        warnings.append("synonym resolution collapses two or more tree tips to the same accepted label")
    renamed_labels = {row.raw_label for row in renames}
    return resolved_tree, TaxonSynonymResolutionReport(
        synonym_table_path=synonym_table_path,
        resolution_policy=resolution_policy,
        renamed_taxa=sorted(renames, key=lambda row: row.raw_label),
        unchanged_taxa=sorted(label for label in tree.tip_names if label not in renamed_labels),
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
    if any(token in lowered for token in _SAMPLE_TOKENS) or re.fullmatch(r"[A-Za-z]+[-_]\d{2,}", normalized):
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
    assignments = sorted((infer_taxon_namespace(label) for label in tree.tip_names), key=lambda row: row.label)
    namespace_counts = dict(sorted(Counter(row.namespace for row in assignments).items()))
    explicit_namespaces = {namespace for namespace in namespace_counts if namespace != "user_defined_label"}
    mixed_namespaces = len(explicit_namespaces) > 1
    warnings: list[str] = []
    if mixed_namespaces:
        warnings.append("tree mixes multiple explicit taxon namespaces without an explicit mapping layer")
    if explicit_namespaces and "user_defined_label" in namespace_counts and len(namespace_counts) > 1:
        warnings.append("tree mixes structured namespaces with arbitrary user-defined labels")
    dominant_namespace = None if not namespace_counts else max(namespace_counts.items(), key=lambda item: (item[1], item[0]))[0]
    return TaxonNamespaceReport(
        assignments=assignments,
        namespace_counts=namespace_counts,
        dominant_namespace=dominant_namespace,
        mixed_namespaces=mixed_namespaces,
        warnings=warnings,
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


def write_synonym_resolution_mapping(path: Path, report: TaxonSynonymResolutionReport) -> Path:
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
