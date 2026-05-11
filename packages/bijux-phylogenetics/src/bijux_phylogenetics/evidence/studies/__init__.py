"""Study-specific evidence-book generators."""

from .primate_longevity_signal import (
    build_primate_claim_registry,
    build_primate_family_index,
    build_primate_parity_policy,
    build_primate_scalar_parity_table,
    build_primate_source_fragment_map,
    build_primate_summary_bundle_claims,
    render_primate_scalar_parity_table_markdown,
)
from .primate_pcm1_component_bundles import (
    build_primate_data_preparation_bundle_index,
    build_primate_pcm1_component_bundles,
    build_primate_structural_parity_table,
    render_primate_data_preparation_bundle_index_markdown,
    render_primate_structural_parity_table_markdown,
)
from .primate_pgls_and_signal import (
    build_primate_pgls_signal_bundles,
    build_primate_pgls_signal_claim_registry,
    build_primate_pgls_signal_external_sources,
    build_primate_pgls_signal_family_index,
    build_primate_pgls_signal_parity_policy,
    build_primate_pgls_signal_scalar_parity_table,
    build_primate_pgls_signal_source_fragment_map,
    render_primate_pgls_signal_scalar_parity_table_markdown,
)

__all__ = [
    "build_primate_data_preparation_bundle_index",
    "build_primate_pcm1_component_bundles",
    "build_primate_claim_registry",
    "build_primate_family_index",
    "build_primate_parity_policy",
    "build_primate_pgls_signal_bundles",
    "build_primate_pgls_signal_claim_registry",
    "build_primate_pgls_signal_external_sources",
    "build_primate_pgls_signal_family_index",
    "build_primate_pgls_signal_parity_policy",
    "build_primate_pgls_signal_scalar_parity_table",
    "build_primate_pgls_signal_source_fragment_map",
    "build_primate_scalar_parity_table",
    "build_primate_summary_bundle_claims",
    "build_primate_source_fragment_map",
    "build_primate_structural_parity_table",
    "render_primate_data_preparation_bundle_index_markdown",
    "render_primate_pgls_signal_scalar_parity_table_markdown",
    "render_primate_structural_parity_table_markdown",
    "render_primate_scalar_parity_table_markdown",
]
