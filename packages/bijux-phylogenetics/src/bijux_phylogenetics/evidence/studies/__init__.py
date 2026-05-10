"""Study-specific evidence-book generators."""

from .primate_longevity_signal import (
    build_primate_claim_registry,
    build_primate_family_index,
    build_primate_source_fragment_map,
)

__all__ = [
    "build_primate_claim_registry",
    "build_primate_family_index",
    "build_primate_source_fragment_map",
]
