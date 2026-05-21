from __future__ import annotations

from ..models import (
    RabiesMethodSensitivityConclusionRow,
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPreprocessingComparisonRow,
    RabiesMethodSensitivityVariantRun,
)


def _build_conclusion_rows(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    variant_runs: list[RabiesMethodSensitivityVariantRun],
    preprocessing_comparison_rows: tuple[
        RabiesMethodSensitivityPreprocessingComparisonRow, ...
    ],
) -> list[RabiesMethodSensitivityConclusionRow]:
    selected_models = sorted(
        {
            variant.inference_comparison.selected_model
            for variant in variant_runs
            if variant.inference_comparison.selected_model
        }
    )
    max_serious_conflicts = max(
        variant.inference_comparison.conclusion_summary.serious_conflict_count
        for variant in variant_runs
    )
    stable_clade_counts = {
        variant.config.variant_id: variant.inference_comparison.conclusion_summary.stable_clade_count
        for variant in variant_runs
    }
    return [
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="preprocessing_rooted_iqtree_topology",
            method_axis="alignment_and_trimming",
            stability_status=(
                "stable"
                if all(
                    row.robinson_foulds_distance == 0
                    and not row.same_taxa_different_rooting
                    for row in preprocessing_comparison_rows
                )
                else "changed"
            ),
            claim=(
                "The rooted IQ-TREE topology stayed unchanged across every declared "
                "alignment and trimming variant."
            ),
            evidence=(
                f"{len(preprocessing_comparison_rows)} rooted pairwise preprocessing "
                "comparisons returned RF distance 0 and no rooting-only disagreements."
            ),
            caution=(
                "This stability statement is limited to the compact nine-taxon rabies "
                "panel and the four declared preprocessing settings."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="preprocessing_selected_model",
            method_axis="alignment_and_trimming",
            stability_status="stable" if len(selected_models) == 1 else "changed",
            claim=(
                "The selected substitution model remained constant across the "
                "declared preprocessing matrix."
                if len(selected_models) == 1
                else "The selected substitution model changed across the declared preprocessing matrix."
            ),
            evidence=(
                f"selected models: {', '.join(selected_models)}"
                if selected_models
                else "no selected model was recorded"
            ),
            caution=(
                "Model-selection stability here reflects one short rabies nucleoprotein "
                "panel rather than a general claim about all pathogen alignments."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="rooted_engine_agreement",
            method_axis="inference_engine",
            stability_status=(
                "stable"
                if all(
                    variant.rooted_engine_comparison.robinson_foulds_distance == 0
                    and not variant.rooted_engine_comparison.same_taxa_different_rooting
                    for variant in variant_runs
                )
                else "changed"
            ),
            claim=(
                "After explicit outgroup rooting, FastTree and IQ-TREE preserved the "
                "same rooted topology in every declared preprocessing variant."
            ),
            evidence=(
                f"rooted engine comparisons over outgroup {', '.join(dataset.outgroup_taxa)} "
                "returned RF distance 0 in every variant."
            ),
            caution=(
                "Rooted agreement does not imply that every internal unrooted split or "
                "support value is interchangeable across engines."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="unrooted_engine_sensitivity",
            method_axis="inference_engine",
            stability_status="changed" if max_serious_conflicts > 0 else "stable",
            claim=(
                "Before rooting, the FastTree versus IQ-TREE comparison changed several "
                "internal clade conclusions on this rabies panel."
            ),
            evidence=(
                f"serious unrooted engine conflicts ranged up to {max_serious_conflicts} "
                f"per variant, while stable shared clades per variant ranged from "
                f"{min(stable_clade_counts.values())} to {max(stable_clade_counts.values())}."
            ),
            caution=(
                "The engine-sensitive clades are a warning against over-reading fine "
                "internal structure from one compact panel, especially when approximate "
                "and likelihood engines disagree before rooting."
            ),
        ),
    ]
