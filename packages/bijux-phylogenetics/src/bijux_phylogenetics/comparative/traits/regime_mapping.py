from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.render.tree_svg import TreeRenderResult, render_tree_svg
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_DEFAULT_REGIME_COLORS = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


@dataclass(slots=True)
class BranchIdentityMetadata:
    """Deterministic identity and exposure context for one non-root branch."""

    branch_id: str
    child_node_name: str | None
    is_tip_branch: bool
    branch_length: float
    descendant_taxa: list[str]
    analyzed_descendant_taxa: list[str]
    contributes_to_analysis: bool


@dataclass(slots=True)
class TraitRegimeExclusion:
    """One taxon excluded before regime mapping from tip states."""

    taxon: str
    reason: str


@dataclass(slots=True)
class TraitRegimeNodeRow:
    """One node-level regime assignment from tip-state reconstruction."""

    node_id: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    regime: str
    candidate_regimes: list[str]
    assignment_confidence: float
    ambiguous_assignment: bool
    state_probabilities: dict[str, float]


@dataclass(slots=True)
class TraitRegimeBranchRow:
    """One normalized branch-to-regime assignment."""

    branch_id: str
    child_node_name: str | None
    is_tip_branch: bool
    branch_length: float
    regime: str
    candidate_regimes: list[str]
    assignment_confidence: float
    ambiguous_assignment: bool
    assignment_origin: str
    descendant_taxa: list[str]
    analyzed_descendant_taxa: list[str]
    contributes_to_analysis: bool


@dataclass(slots=True)
class TraitRegimeMappingReport:
    """Reviewer-facing regime map for comparative evolutionary workflows."""

    tree_path: Path
    source_path: Path
    source_kind: str
    trait: str | None
    taxon_column: str | None
    reconstruction_model: str | None
    state_ordering: str | None
    ordered_states: list[str]
    branch_id_column: str
    regime_column: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[TraitRegimeExclusion]
    observed_regimes: list[str]
    branch_rows: list[TraitRegimeBranchRow]
    node_rows: list[TraitRegimeNodeRow]
    ambiguous_branch_count: int
    warnings: list[str]
    analysis_tree_newick: str | None


def summarize_trait_regime_mapping(
    tree_path: Path,
    *,
    tip_states_path: Path | None = None,
    regime_map_path: Path | None = None,
    trait: str | None = None,
    taxon_column: str | None = None,
    reconstruction_model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    branch_id_column: str | None = None,
    regime_column: str = "regime",
) -> TraitRegimeMappingReport:
    """Reconstruct or normalize a branch regime map for comparative workflows."""
    if (tip_states_path is None) == (regime_map_path is None):
        raise ComparativeMethodError(
            "trait regime mapping requires exactly one source: tip states or user-provided regime map"
        )
    if tip_states_path is not None:
        if trait is None:
            raise ComparativeMethodError(
                "tip-state regime mapping requires a discrete trait column name"
            )
        if reconstruction_model == "fitch" and state_ordering != "unordered":
            raise ComparativeMethodError(
                "ordered regime reconstruction requires a likelihood model"
            )
        return _summarize_trait_regime_mapping_from_tip_states(
            tree_path,
            tip_states_path,
            trait=trait,
            taxon_column=taxon_column,
            reconstruction_model=reconstruction_model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
    return _summarize_trait_regime_mapping_from_map(
        tree_path,
        regime_map_path,
        branch_id_column=branch_id_column,
        regime_column=regime_column,
    )


def build_branch_identity_lookup(
    tree: PhyloTree,
    *,
    analyzed_taxa: list[str] | None = None,
) -> dict[str, BranchIdentityMetadata]:
    """Build stable branch identities for every non-root branch in a tree."""
    analyzed_set = None if analyzed_taxa is None else set(analyzed_taxa)
    rows: dict[str, BranchIdentityMetadata] = {}
    for node in tree.iter_nodes():
        if node is tree.root:
            continue
        branch_id = node_signature(node)
        descendants = node_descendant_taxa(node)
        analyzed_descendants = (
            list(descendants)
            if analyzed_set is None
            else [taxon for taxon in descendants if taxon in analyzed_set]
        )
        rows[branch_id] = BranchIdentityMetadata(
            branch_id=branch_id,
            child_node_name=node.name,
            is_tip_branch=node.is_leaf(),
            branch_length=float(node.branch_length or 0.0),
            descendant_taxa=descendants,
            analyzed_descendant_taxa=analyzed_descendants,
            contributes_to_analysis=bool(analyzed_descendants),
        )
    return rows


def resolve_branch_regime_id_column(path: Path, *, requested: str | None) -> str:
    """Resolve the branch-identity column for a regime map table."""
    if requested is not None:
        return requested
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    delimiter = "\t" if ("\t" in first_line or path.suffix.lower() == ".tsv") else ","
    columns = [
        column.strip() for column in next(csv.reader([first_line], delimiter=delimiter))
    ]
    for candidate in ("branch_id", "branch", "node"):
        if candidate in columns:
            return candidate
    return columns[0]


def render_trait_regime_map(
    report: TraitRegimeMappingReport,
    *,
    out_path: Path,
    layout: str = "cladogram",
) -> TreeRenderResult:
    """Render one SVG tree that visualizes the branch regime assignment."""
    tip_regimes = {
        row.child_node_name: row.regime
        for row in report.branch_rows
        if row.is_tip_branch and row.child_node_name is not None
    }
    internal_annotations = {
        row.branch_id: row.regime
        if not row.ambiguous_assignment
        else "/".join(row.candidate_regimes)
        for row in report.branch_rows
        if not row.is_tip_branch
    }
    palette = {
        regime: _DEFAULT_REGIME_COLORS[index % len(_DEFAULT_REGIME_COLORS)]
        for index, regime in enumerate(sorted(report.observed_regimes))
    }
    internal_annotation_colors = {
        row.branch_id: palette.get(row.regime, "#6d28d9")
        for row in report.branch_rows
        if not row.is_tip_branch
    }
    render_tree_path = report.tree_path
    temporary_render_path: Path | None = None
    if report.analysis_tree_newick is not None:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            suffix=".nwk",
            delete=False,
        ) as handle:
            handle.write(report.analysis_tree_newick)
            temporary_render_path = Path(handle.name)
        render_tree_path = temporary_render_path
    try:
        return render_tree_svg(
            render_tree_path,
            out_path=out_path,
            layout=layout,
            categorical_traits=tip_regimes,
            internal_annotations=internal_annotations,
            internal_annotation_colors=internal_annotation_colors,
        )
    finally:
        if temporary_render_path is not None:
            temporary_render_path.unlink(missing_ok=True)


def _summarize_trait_regime_mapping_from_tip_states(
    tree_path: Path,
    tip_states_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    reconstruction_model: str,
    state_ordering: str,
    ordered_states: list[str] | None,
) -> TraitRegimeMappingReport:
    source_tree = load_tree(tree_path)
    table = load_taxon_table(tip_states_path, taxon_column=taxon_column)
    dataset = load_discrete_dataset(
        tree_path,
        tip_states_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    reconstruction = reconstruct_discrete_ancestral_states(
        tree_path,
        tip_states_path,
        trait=trait,
        taxon_column=taxon_column,
        model=reconstruction_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    branch_lookup = build_branch_identity_lookup(
        dataset.tree,
        analyzed_taxa=dataset.taxa,
    )
    node_rows = [
        TraitRegimeNodeRow(
            node_id=estimate.node,
            node_name=estimate.node_name,
            is_tip=estimate.is_tip,
            descendant_taxa=list(estimate.descendant_taxa),
            regime=estimate.most_likely_state,
            candidate_regimes=list(estimate.state_set),
            assignment_confidence=estimate.confidence,
            ambiguous_assignment=estimate.ambiguous,
            state_probabilities=dict(estimate.state_probabilities),
        )
        for estimate in reconstruction.estimates
    ]
    branch_rows = [
        TraitRegimeBranchRow(
            branch_id=estimate.node,
            child_node_name=estimate.node_name,
            is_tip_branch=estimate.is_tip,
            branch_length=branch_lookup[estimate.node].branch_length,
            regime=estimate.most_likely_state,
            candidate_regimes=list(estimate.state_set),
            assignment_confidence=estimate.confidence,
            ambiguous_assignment=estimate.ambiguous,
            assignment_origin="tip-state-reconstruction",
            descendant_taxa=list(estimate.descendant_taxa),
            analyzed_descendant_taxa=list(
                branch_lookup[estimate.node].analyzed_descendant_taxa
            ),
            contributes_to_analysis=branch_lookup[
                estimate.node
            ].contributes_to_analysis,
        )
        for estimate in reconstruction.estimates
        if estimate.node in branch_lookup
    ]
    warnings = list(reconstruction.warnings)
    state_collisions = sorted(
        {
            state
            for state in reconstruction.observed_states
            if any(token in state for token in ("|", "/", ";"))
        }
    )
    if state_collisions:
        warnings.append(
            "one or more observed regime states contain reserved display separators and may be harder to interpret in branch tables"
        )
    return TraitRegimeMappingReport(
        tree_path=tree_path,
        source_path=tip_states_path,
        source_kind="tip-state-reconstruction",
        trait=trait,
        taxon_column=dataset.taxon_column,
        reconstruction_model=reconstruction_model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        branch_id_column="branch_id",
        regime_column="regime",
        tree_taxon_count=len(source_tree.tip_names),
        analyzed_taxa=list(dataset.taxa),
        analyzed_taxon_count=len(dataset.taxa),
        excluded_taxa=_build_tip_state_exclusions(
            source_tree=source_tree,
            table=table,
            dataset=dataset,
            trait=trait,
        ),
        observed_regimes=list(reconstruction.observed_states),
        branch_rows=sorted(branch_rows, key=lambda row: row.branch_id),
        node_rows=sorted(node_rows, key=lambda row: row.node_id),
        ambiguous_branch_count=sum(
            1 for row in branch_rows if row.ambiguous_assignment
        ),
        warnings=warnings,
        analysis_tree_newick=reconstruction.analysis_tree_newick,
    )


def _summarize_trait_regime_mapping_from_map(
    tree_path: Path,
    regime_map_path: Path,
    *,
    branch_id_column: str | None,
    regime_column: str,
) -> TraitRegimeMappingReport:
    tree = load_tree(tree_path)
    branch_lookup = build_branch_identity_lookup(tree, analyzed_taxa=tree.tip_names)
    resolved_branch_id_column = resolve_branch_regime_id_column(
        regime_map_path,
        requested=branch_id_column,
    )
    table = load_taxon_table(regime_map_path, taxon_column=resolved_branch_id_column)
    if regime_column not in table.columns:
        raise ComparativeMethodError(
            f"regime map does not contain column '{regime_column}'"
        )
    mapped_branch_ids = {row[table.taxon_column] for row in table.rows}
    expected_branch_ids = set(branch_lookup)
    missing = sorted(expected_branch_ids - mapped_branch_ids)
    if missing:
        raise ComparativeMethodError(
            "regime map is missing one or more non-root branches: " + ", ".join(missing)
        )
    extra = sorted(mapped_branch_ids - expected_branch_ids)
    if extra:
        raise ComparativeMethodError(
            "regime map contains branches absent from the tree: " + ", ".join(extra)
        )
    branch_rows: list[TraitRegimeBranchRow] = []
    for row in table.rows:
        branch = branch_lookup[row[table.taxon_column]]
        regime = row[regime_column]
        if not regime:
            raise ComparativeMethodError(
                f"regime map branch '{branch.branch_id}' has an empty '{regime_column}' value"
            )
        branch_rows.append(
            TraitRegimeBranchRow(
                branch_id=branch.branch_id,
                child_node_name=branch.child_node_name,
                is_tip_branch=branch.is_tip_branch,
                branch_length=branch.branch_length,
                regime=regime,
                candidate_regimes=[regime],
                assignment_confidence=1.0,
                ambiguous_assignment=False,
                assignment_origin="user-provided-map",
                descendant_taxa=list(branch.descendant_taxa),
                analyzed_descendant_taxa=list(branch.analyzed_descendant_taxa),
                contributes_to_analysis=branch.contributes_to_analysis,
            )
        )
    return TraitRegimeMappingReport(
        tree_path=tree_path,
        source_path=regime_map_path,
        source_kind="user-provided-map",
        trait=None,
        taxon_column=None,
        reconstruction_model=None,
        state_ordering=None,
        ordered_states=[],
        branch_id_column=resolved_branch_id_column,
        regime_column=regime_column,
        tree_taxon_count=len(tree.tip_names),
        analyzed_taxa=list(tree.tip_names),
        analyzed_taxon_count=len(tree.tip_names),
        excluded_taxa=[],
        observed_regimes=sorted({row.regime for row in branch_rows}),
        branch_rows=sorted(branch_rows, key=lambda row: row.branch_id),
        node_rows=[],
        ambiguous_branch_count=0,
        warnings=[],
        analysis_tree_newick=None,
    )


def _build_tip_state_exclusions(
    *,
    source_tree: PhyloTree,
    table,
    dataset,
    trait: str,
) -> list[TraitRegimeExclusion]:
    rows: list[TraitRegimeExclusion] = []
    rows.extend(
        TraitRegimeExclusion(taxon=taxon, reason="missing_from_state_table")
        for taxon in sorted(
            set(source_tree.tip_names) - {row[table.taxon_column] for row in table.rows}
        )
    )
    rows.extend(
        TraitRegimeExclusion(taxon=taxon, reason="missing_state_value")
        for taxon in dataset.dropped_missing_taxa
    )
    rows.extend(
        TraitRegimeExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in sorted(set(table.taxa) - set(source_tree.tip_names))
    )
    return rows


def write_trait_regime_summary_table(
    path: Path,
    report: TraitRegimeMappingReport,
) -> Path:
    """Write one summary ledger for a trait regime map."""
    return write_taxon_rows(
        path,
        columns=[
            "source_kind",
            "trait",
            "taxon_column",
            "reconstruction_model",
            "state_ordering",
            "ordered_states",
            "branch_id_column",
            "regime_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "branch_count",
            "regime_count",
            "ambiguous_branch_count",
            "node_count",
        ],
        rows=[
            {
                "source_kind": report.source_kind,
                "trait": report.trait or "",
                "taxon_column": report.taxon_column or "",
                "reconstruction_model": report.reconstruction_model or "",
                "state_ordering": report.state_ordering or "",
                "ordered_states": ",".join(report.ordered_states),
                "branch_id_column": report.branch_id_column,
                "regime_column": report.regime_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "branch_count": len(report.branch_rows),
                "regime_count": len(report.observed_regimes),
                "ambiguous_branch_count": report.ambiguous_branch_count,
                "node_count": len(report.node_rows),
            }
        ],
    )


def write_trait_regime_branch_table(
    path: Path,
    report: TraitRegimeMappingReport,
) -> Path:
    """Write one normalized branch regime table."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "child_node_name",
            "is_tip_branch",
            "branch_length",
            "regime",
            "candidate_regimes",
            "assignment_confidence",
            "ambiguous_assignment",
            "assignment_origin",
            "descendant_taxa",
            "analyzed_descendant_taxa",
            "contributes_to_analysis",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "child_node_name": row.child_node_name or "",
                "is_tip_branch": str(row.is_tip_branch).lower(),
                "branch_length": format(row.branch_length, ".15g"),
                "regime": row.regime,
                "candidate_regimes": ",".join(row.candidate_regimes),
                "assignment_confidence": format(row.assignment_confidence, ".15g"),
                "ambiguous_assignment": str(row.ambiguous_assignment).lower(),
                "assignment_origin": row.assignment_origin,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "analyzed_descendant_taxa": ",".join(row.analyzed_descendant_taxa),
                "contributes_to_analysis": str(row.contributes_to_analysis).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_trait_regime_node_table(
    path: Path,
    report: TraitRegimeMappingReport,
) -> Path:
    """Write one node-level state reconstruction ledger for a regime map."""
    return write_taxon_rows(
        path,
        columns=[
            "node_id",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "regime",
            "candidate_regimes",
            "assignment_confidence",
            "ambiguous_assignment",
            "state_probabilities",
        ],
        rows=[
            {
                "node_id": row.node_id,
                "node_name": row.node_name or "",
                "is_tip": str(row.is_tip).lower(),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "regime": row.regime,
                "candidate_regimes": ",".join(row.candidate_regimes),
                "assignment_confidence": format(row.assignment_confidence, ".15g"),
                "ambiguous_assignment": str(row.ambiguous_assignment).lower(),
                "state_probabilities": json.dumps(
                    row.state_probabilities,
                    sort_keys=True,
                ),
            }
            for row in report.node_rows
        ],
    )


def write_trait_regime_exclusion_table(
    path: Path,
    report: TraitRegimeMappingReport,
) -> Path:
    """Write one excluded-taxon ledger for a tip-state regime map."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )
