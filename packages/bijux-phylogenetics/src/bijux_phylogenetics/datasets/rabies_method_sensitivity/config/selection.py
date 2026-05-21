from __future__ import annotations

from dataclasses import replace

from ..models import RabiesMethodSensitivityPanelDataset, RabiesMethodSensitivityVariant


def _resolve_selected_variant_dataset(
    dataset: RabiesMethodSensitivityPanelDataset,
    *,
    variant_ids: tuple[str, ...] | None,
) -> RabiesMethodSensitivityPanelDataset:
    """Return either the full dataset or an explicit variant-scoped subset."""
    if variant_ids is None:
        return dataset
    if not variant_ids:
        raise ValueError("variant_ids must not be empty when provided")
    variants_by_id = {variant.variant_id: variant for variant in dataset.variants}
    selected_variants: list[RabiesMethodSensitivityVariant] = []
    seen_variant_ids: set[str] = set()
    for variant_id in variant_ids:
        if variant_id in seen_variant_ids:
            raise ValueError(f"duplicate variant_id requested: {variant_id}")
        seen_variant_ids.add(variant_id)
        variant = variants_by_id.get(variant_id)
        if variant is None:
            known = ", ".join(sorted(variants_by_id))
            raise ValueError(
                f"unknown variant_id '{variant_id}'; known variants: {known}"
            )
        selected_variants.append(variant)
    return replace(dataset, variants=tuple(selected_variants))
