from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.render.tree_figure_package import TreeFigurePackageResult


@dataclass(frozen=True, slots=True)
class AnnotatedTraitTreeCoverageRow:
    """Coverage audit for one annotation surface on a publication trait tree."""

    surface: str
    source_kind: str
    required: bool
    visible_taxon_count: int
    observed_taxon_count: int
    covered_taxon_count: int
    complete: bool
    missing_taxa: list[str]
    extra_taxa: list[str]


@dataclass(frozen=True, slots=True)
class AnnotatedTraitTreeSummaryRow:
    """Stable summary for one rendered label, trait, or metadata surface."""

    surface: str
    source_kind: str
    value_kind: str
    observed_taxon_count: int
    missing_taxon_count: int
    distinct_value_count: int
    minimum_numeric_value: float | None
    maximum_numeric_value: float | None
    example_values: list[str]


@dataclass(frozen=True, slots=True)
class AnnotatedTraitTreePublicationAudit:
    """Reviewer-facing publication audit for an annotated trait tree package."""

    publication_ready: bool
    required_surface_count: int
    complete_surface_count: int
    missing_surface_count: int
    legend_entry_count: int
    caption_ready: bool
    legible: bool
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class AnnotatedTraitTreePackageResult:
    """Materialized annotated trait tree package and its review artifacts."""

    output_dir: Path
    figure_package: TreeFigurePackageResult
    coverage_path: Path
    summary_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    coverage_rows: list[AnnotatedTraitTreeCoverageRow]
    summary_rows: list[AnnotatedTraitTreeSummaryRow]
    audit: AnnotatedTraitTreePublicationAudit
