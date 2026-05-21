"""Owned publication package surface for annotated trait trees."""

from __future__ import annotations

from .contracts import (
    AnnotatedTraitTreeCoverageRow,
    AnnotatedTraitTreePackageResult,
    AnnotatedTraitTreePublicationAudit,
    AnnotatedTraitTreeSummaryRow,
)
from .inputs import (
    build_annotation_strips,
    build_full_label_map,
    build_numeric_map,
    build_string_map,
    require_table,
)
from .summaries import (
    build_coverage_row,
    build_heatmap_summary_row,
    build_label_summary_row,
    build_numeric_summary_row,
    build_string_summary_row,
)

__all__ = [
    "AnnotatedTraitTreeCoverageRow",
    "AnnotatedTraitTreePackageResult",
    "AnnotatedTraitTreePublicationAudit",
    "AnnotatedTraitTreeSummaryRow",
    "build_annotation_strips",
    "build_coverage_row",
    "build_full_label_map",
    "build_heatmap_summary_row",
    "build_label_summary_row",
    "build_numeric_map",
    "build_numeric_summary_row",
    "build_string_map",
    "build_string_summary_row",
    "require_table",
]
