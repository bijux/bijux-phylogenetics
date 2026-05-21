from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.pruning import drop_tree_taxa
from bijux_phylogenetics.datasets.study_inputs import (
    link_tree_to_traits,
    validate_traits_table,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError
from bijux_phylogenetics.io.fasta import load_fasta_alignment, write_fasta_alignment
from bijux_phylogenetics.io.newick import write_newick

from ..models import (
    CatarrhineDataQualityStressPanelDataset,
    CatarrhineDataQualityStressPanelWorkflowReport,
    DataQualityRepairAction,
    TraitDuplicateResolution,
)
from ..panel import TREE_BRANCH_FLOOR
from ..traits import (
    detect_missing_traits,
    load_permissive_trait_rows,
    resolve_duplicate_traits,
    selected_trait_rows,
)
from .sequence_processing import prepare_sequence_surfaces
from .shared import apply_branch_length_floor


def build_catarrhine_data_quality_stress_panel_workflow_report(
    dataset: CatarrhineDataQualityStressPanelDataset,
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    out_dir.mkdir(parents=True, exist_ok=True)
    assembled_root = out_dir / "assembled"
    assembled_root.mkdir(parents=True, exist_ok=True)

    sequence_surfaces = prepare_sequence_surfaces(
        dataset=dataset,
        assembled_root=assembled_root,
    )
    raw_trait_mismatch_linkage = link_tree_to_traits(
        dataset.raw_tree_path,
        dataset.raw_trait_mismatch_path,
        strict=False,
    )
    raw_trait_mismatch_error: str | None = None
    try:
        link_tree_to_traits(
            dataset.raw_tree_path,
            dataset.raw_trait_mismatch_path,
            strict=True,
        )
    except MetadataJoinError as error:
        raw_trait_mismatch_error = str(error)

    raw_trait_rows = load_permissive_trait_rows(dataset.raw_traits_path)
    trait_duplicates = resolve_duplicate_traits(raw_trait_rows)
    duplicate_lookup = {row.taxon: row for row in trait_duplicates}
    cleaned_trait_candidates = selected_trait_rows(raw_trait_rows)
    missing_traits = detect_missing_traits(
        raw_trait_rows,
        required_traits=set(dataset.required_traits),
        duplicate_lookup=duplicate_lookup,
    )

    missing_required_trait_taxa = sorted(
        {
            row.taxon
            for row in missing_traits
            if row.required_for_analysis
            and row.action == "drop_taxon_from_cleaned_traits"
        }
    )
    excluded_taxa = sorted(
        set(sequence_surfaces.sequence_outlier_taxa)
        | set(sequence_surfaces.tree_outlier_taxa)
        | set(missing_required_trait_taxa)
    )

    cleaned_trait_rows = [
        dict(row)
        for row in cleaned_trait_candidates
        if row["taxon"] not in excluded_taxa
        and all(row[trait] for trait in dataset.required_traits)
    ]
    cleaned_traits_path = write_taxon_rows(
        assembled_root / "cleaned-traits.csv",
        columns=list(cleaned_trait_rows[0].keys()),
        rows=cleaned_trait_rows,
    )

    cleaned_tree, _ = drop_tree_taxa(dataset.raw_tree_path, excluded_taxa)
    repaired_branch_nodes = apply_branch_length_floor(
        cleaned_tree.root,
        floor=TREE_BRANCH_FLOOR,
    )
    cleaned_tree_path = write_newick(assembled_root / "cleaned-tree.nwk", cleaned_tree)

    cleaned_alignment_records = [
        record
        for record in load_fasta_alignment(dataset.raw_alignment_path)
        if record.identifier not in excluded_taxa
    ]
    cleaned_alignment_path = write_fasta_alignment(
        assembled_root / "cleaned-alignment.fasta",
        cleaned_alignment_records,
    )

    cleaned_trait_validation = validate_traits_table(cleaned_traits_path)
    cleaned_tree_validation = sequence_surfaces.validate_cleaned_tree(cleaned_tree_path)
    cleaned_alignment_validation = sequence_surfaces.validate_cleaned_alignment(
        cleaned_alignment_path
    )
    cleaned_linkage = link_tree_to_traits(
        cleaned_tree_path,
        cleaned_traits_path,
        strict=True,
    )
    cleaned_taxa = sorted(cleaned_linkage.usable_taxa)

    repair_actions = build_repair_actions(
        raw_sequence_input_repair=sequence_surfaces.raw_sequence_input_repair,
        raw_sequence_length_outliers=sequence_surfaces.raw_sequence_length_outliers,
        coding_sequence_preparation=sequence_surfaces.coding_sequence_preparation,
        raw_trait_mismatch_error=raw_trait_mismatch_error,
        trait_duplicates=trait_duplicates,
        missing_required_trait_taxa=missing_required_trait_taxa,
        sequence_outlier_taxa=sequence_surfaces.sequence_outlier_taxa,
        tree_outlier_taxa=sequence_surfaces.tree_outlier_taxa,
        repaired_branch_nodes=repaired_branch_nodes,
    )
    return CatarrhineDataQualityStressPanelWorkflowReport(
        dataset=dataset,
        raw_sequence_input_validation=sequence_surfaces.raw_sequence_input_validation,
        raw_sequence_input_repair=sequence_surfaces.raw_sequence_input_repair,
        raw_sequence_length_outliers=sequence_surfaces.raw_sequence_length_outliers,
        repaired_sequence_input_validation=(
            sequence_surfaces.repaired_sequence_input_validation
        ),
        coding_sequence_preparation=sequence_surfaces.coding_sequence_preparation,
        raw_trait_mismatch_linkage=raw_trait_mismatch_linkage,
        raw_trait_mismatch_error=raw_trait_mismatch_error,
        raw_alignment_validation=sequence_surfaces.raw_alignment_validation,
        sequence_outliers=sequence_surfaces.sequence_outliers,
        raw_tree_inspection=sequence_surfaces.raw_tree_inspection,
        raw_tree_validation=sequence_surfaces.raw_tree_validation,
        trait_duplicates=trait_duplicates,
        missing_traits=missing_traits,
        cleaned_trait_validation=cleaned_trait_validation,
        cleaned_tree_validation=cleaned_tree_validation,
        cleaned_linkage=cleaned_linkage,
        cleaned_alignment_validation=cleaned_alignment_validation,
        cleaned_alignment_records=cleaned_alignment_records,
        repaired_sequence_input_path=sequence_surfaces.repaired_sequence_input_path,
        prepared_coding_sequences_path=sequence_surfaces.prepared_coding_sequences_path,
        cleaned_tree_path=cleaned_tree_path,
        cleaned_traits_path=cleaned_traits_path,
        cleaned_alignment_path=cleaned_alignment_path,
        cleaned_taxa=cleaned_taxa,
        dropped_taxa=excluded_taxa,
        repair_actions=repair_actions,
        repaired_branch_nodes=repaired_branch_nodes,
    )
def build_repair_actions(
    *,
    raw_sequence_input_repair,
    raw_sequence_length_outliers,
    coding_sequence_preparation,
    raw_trait_mismatch_error: str | None,
    trait_duplicates: list[TraitDuplicateResolution],
    missing_required_trait_taxa: list[str],
    sequence_outlier_taxa: list[str],
    tree_outlier_taxa: list[str],
    repaired_branch_nodes: list[str],
) -> list[DataQualityRepairAction]:
    actions: list[DataQualityRepairAction] = []
    if raw_sequence_input_repair.normalized_identifiers:
        actions.append(
            DataQualityRepairAction(
                action_kind="normalize_duplicate_sequence_identifiers",
                affected_taxa=[
                    row.repaired_identifier
                    for row in raw_sequence_input_repair.normalized_identifiers
                ],
                affected_nodes=[],
                reason="raw FASTA identifiers were duplicated and required deterministic renaming",
                result="rewrote retained duplicate identifiers into unique names",
            )
        )
    if raw_sequence_input_repair.removed_records:
        actions.append(
            DataQualityRepairAction(
                action_kind="remove_invalid_sequence_records",
                affected_taxa=[
                    row.identifier for row in raw_sequence_input_repair.removed_records
                ],
                affected_nodes=[],
                reason="raw FASTA records contained empty sequences or illegal characters",
                result="removed invalid records before downstream reuse",
            )
        )
    if raw_sequence_length_outliers:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_sequence_length_outliers",
                affected_taxa=[row.identifier for row in raw_sequence_length_outliers],
                affected_nodes=[],
                reason="raw FASTA sequence lengths deviated too far from the dataset baseline",
                result="excluded length outliers from the repaired FASTA subset",
            )
        )
    if coding_sequence_preparation.excluded_sequences:
        actions.append(
            DataQualityRepairAction(
                action_kind="exclude_invalid_coding_sequences",
                affected_taxa=[
                    row.identifier
                    for row in coding_sequence_preparation.excluded_sequences
                ],
                affected_nodes=[],
                reason="coding-sequence preparation found frame errors or premature stop codons",
                result="retained only translation-ready coding sequences",
            )
        )
    if raw_trait_mismatch_error:
        actions.append(
            DataQualityRepairAction(
                action_kind="reject_raw_trait_linkage_mismatch",
                affected_taxa=[],
                affected_nodes=[],
                reason="raw trait linkage omitted one tree taxon and introduced one extra trait taxon",
                result="strict raw linkage failed and the mismatched table was kept only as a review surface",
            )
        )
    for row in trait_duplicates:
        actions.append(
            DataQualityRepairAction(
                action_kind="resolve_duplicate_trait_taxon",
                affected_taxa=[row.taxon],
                affected_nodes=[],
                reason="multiple raw trait rows were present for one taxon",
                result=(
                    f"kept row {row.selected_row_number} and discarded "
                    + ",".join(str(value) for value in row.discarded_row_numbers)
                ),
            )
        )
    if missing_required_trait_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_taxa_for_missing_required_traits",
                affected_taxa=missing_required_trait_taxa,
                affected_nodes=[],
                reason="required comparative traits remained missing after duplicate resolution",
                result="excluded taxa from the cleaned comparative subset",
            )
        )
    if sequence_outlier_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_sequence_outlier_taxa",
                affected_taxa=sequence_outlier_taxa,
                affected_nodes=[],
                reason="alignment composition outliers exceeded the governed stress threshold",
                result="excluded taxa from the cleaned alignment and downstream subset",
            )
        )
    if tree_outlier_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_tree_branch_outlier_taxa",
                affected_taxa=tree_outlier_taxa,
                affected_nodes=[],
                reason="terminal branches were extreme relative to the raw tree baseline",
                result="excluded taxa from the cleaned tree and linked surfaces",
            )
        )
    if repaired_branch_nodes:
        actions.append(
            DataQualityRepairAction(
                action_kind="floor_nonpositive_branch_lengths",
                affected_taxa=[],
                affected_nodes=repaired_branch_nodes,
                reason="nonpositive branch lengths block comparative and tree-distance assumptions",
                result=f"replaced nonpositive branch lengths with {TREE_BRANCH_FLOOR}",
            )
        )
    return actions
