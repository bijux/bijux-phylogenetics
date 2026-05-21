from __future__ import annotations

from bijux_phylogenetics.compare.topology import compare_tree_paths

from ..models import (
    RabiesMethodSensitivityPreprocessingComparisonRow,
    RabiesMethodSensitivityVariant,
    RabiesMethodSensitivityVariantRun,
)


def _build_preprocessing_comparison_rows(
    variant_runs: list[RabiesMethodSensitivityVariantRun],
) -> list[RabiesMethodSensitivityPreprocessingComparisonRow]:
    rows: list[RabiesMethodSensitivityPreprocessingComparisonRow] = []
    for index, left in enumerate(variant_runs):
        for right in variant_runs[index + 1 :]:
            comparison = compare_tree_paths(
                left.rooted_iqtree_path, right.rooted_iqtree_path
            )
            rows.append(
                RabiesMethodSensitivityPreprocessingComparisonRow(
                    left_variant_id=left.config.variant_id,
                    right_variant_id=right.config.variant_id,
                    comparison_axis=_comparison_axis(left.config, right.config),
                    robinson_foulds_distance=comparison.robinson_foulds_distance,
                    normalized_robinson_foulds=comparison.normalized_robinson_foulds,
                    same_taxa_different_rooting=comparison.same_taxa_different_rooting,
                )
            )
    return rows


def _comparison_axis(
    left: RabiesMethodSensitivityVariant, right: RabiesMethodSensitivityVariant
) -> str:
    alignment_changed = left.alignment_mode != right.alignment_mode
    trimming_changed = left.trimming_mode != right.trimming_mode
    if alignment_changed and not trimming_changed:
        return "alignment_mode"
    if trimming_changed and not alignment_changed:
        return "trimming_mode"
    return "combined_preprocessing"
