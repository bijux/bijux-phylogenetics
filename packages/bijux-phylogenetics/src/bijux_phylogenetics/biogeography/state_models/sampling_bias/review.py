from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)

from .comparison_assembly import (
    build_exclusion_rows,
    build_node_rows,
    build_summary,
    build_transition_rows,
)
from .contracts import GeographicSamplingBiasReport
from .model_execution import (
    public_model_alias,
    resolve_internal_model,
    run_sampling_bias_model,
)
from .weighting_policy import (
    build_count_rows,
    build_warnings,
    included_region_counts,
    resolve_region_weights,
)


def summarize_geographic_sampling_bias(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
    weights_path: Path | None = None,
    region_column: str = "region",
    weight_column: str = "weight",
) -> GeographicSamplingBiasReport:
    """Review how explicit region weights change biogeographic state inference."""
    internal_model = resolve_internal_model(model)
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_regions,
    )
    exclusion_rows = build_exclusion_rows(audit)
    included_counts = included_region_counts(audit)
    weights, weighting_mode = resolve_region_weights(
        included_counts,
        weights_path=weights_path,
        region_column=region_column,
        weight_column=weight_column,
    )
    baseline = run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=None,
    )
    weighted = run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=weights,
    )
    count_rows = build_count_rows(included_counts, weights)
    node_rows = build_node_rows(baseline, weighted)
    transition_rows = build_transition_rows(baseline, weighted)
    warnings = build_warnings(count_rows, weighting_mode, node_rows)
    summary = build_summary(
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=public_model_alias(internal_model),
        internal_model=internal_model,
        analyzed_taxon_count=sum(included_counts.values()),
        excluded_taxon_count=len(exclusion_rows),
        weighting_mode=weighting_mode,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        warning_count=len(warnings),
    )
    return GeographicSamplingBiasReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=summary.model,
        internal_model=internal_model,
        weighting_mode=weighting_mode,
        summary=summary,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )
