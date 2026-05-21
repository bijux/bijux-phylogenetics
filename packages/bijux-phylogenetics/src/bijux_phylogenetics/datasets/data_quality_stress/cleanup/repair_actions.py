from __future__ import annotations

from bijux_phylogenetics.phylo.alignment import (
    CodingSequencePreparationReport,
    FastaRepairReport,
    SequenceLengthOutlier,
)

from ..models import DataQualityRepairAction, TraitDuplicateResolution
from ..panel import TREE_BRANCH_FLOOR


def build_repair_actions(
    *,
    raw_sequence_input_repair: FastaRepairReport,
    raw_sequence_length_outliers: list[SequenceLengthOutlier],
    coding_sequence_preparation: CodingSequencePreparationReport,
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
