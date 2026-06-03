from __future__ import annotations

from bijux_phylogenetics.render.tree_svg import AnnotationStrip

from .contracts import (
    AnnotatedTraitTreeCoverageRow,
    AnnotatedTraitTreeSummaryRow,
)


def build_coverage_row(
    *,
    surface: str,
    source_kind: str,
    required: bool,
    taxa: list[str],
    observed_taxa: set[str],
) -> AnnotatedTraitTreeCoverageRow:
    """Build one coverage audit row for a rendered annotation surface."""
    tree_taxa = set(taxa)
    missing_taxa = sorted(taxon for taxon in taxa if taxon not in observed_taxa)
    extra_taxa = sorted(observed_taxa - tree_taxa)
    covered_taxa = len(tree_taxa & observed_taxa)
    return AnnotatedTraitTreeCoverageRow(
        surface=surface,
        source_kind=source_kind,
        required=required,
        visible_taxon_count=len(taxa),
        observed_taxon_count=len(observed_taxa),
        covered_taxon_count=covered_taxa,
        complete=not missing_taxa and not extra_taxa,
        missing_taxa=missing_taxa,
        extra_taxa=extra_taxa,
    )


def build_label_summary_row(
    *,
    labels: dict[str, str],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    """Build the summary row for rendered publication labels."""
    unique_labels = sorted(set(labels.values()))
    return AnnotatedTraitTreeSummaryRow(
        surface="labels",
        source_kind="labels",
        value_kind="text",
        observed_taxon_count=len(taxa),
        missing_taxon_count=0,
        distinct_value_count=len(unique_labels),
        minimum_numeric_value=None,
        maximum_numeric_value=None,
        example_values=unique_labels[:3],
    )


def build_string_summary_row(
    *,
    surface: str,
    source_kind: str,
    values: dict[str, str],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    """Build the summary row for one categorical annotation surface."""
    unique_values = sorted(set(values.values()))
    return AnnotatedTraitTreeSummaryRow(
        surface=surface,
        source_kind=source_kind,
        value_kind="categorical",
        observed_taxon_count=len(values),
        missing_taxon_count=len(taxa)
        - len({taxon for taxon in taxa if taxon in values}),
        distinct_value_count=len(unique_values),
        minimum_numeric_value=None,
        maximum_numeric_value=None,
        example_values=unique_values[:4],
    )


def build_numeric_summary_row(
    *,
    surface: str,
    source_kind: str,
    values: dict[str, float],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    """Build the summary row for one numeric annotation surface."""
    observed_values = list(values.values())
    return AnnotatedTraitTreeSummaryRow(
        surface=surface,
        source_kind=source_kind,
        value_kind="numeric",
        observed_taxon_count=len(observed_values),
        missing_taxon_count=len(taxa)
        - len({taxon for taxon in taxa if taxon in values}),
        distinct_value_count=len(set(observed_values)),
        minimum_numeric_value=min(observed_values) if observed_values else None,
        maximum_numeric_value=max(observed_values) if observed_values else None,
        example_values=[
            format(value, ".6g") for value in sorted(set(observed_values))[:4]
        ],
    )


def build_heatmap_summary_row(
    *,
    column: AnnotationStrip,
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    """Build the summary row for one rendered heatmap column."""
    observed_values = [value for value in column.values.values() if value]
    if not observed_values:
        return AnnotatedTraitTreeSummaryRow(
            surface=column.name,
            source_kind="heatmap",
            value_kind="empty",
            observed_taxon_count=0,
            missing_taxon_count=len(taxa),
            distinct_value_count=0,
            minimum_numeric_value=None,
            maximum_numeric_value=None,
            example_values=[],
        )
    try:
        numeric_values = [float(value) for value in observed_values]
    except ValueError:
        unique_values = sorted(set(observed_values))
        return AnnotatedTraitTreeSummaryRow(
            surface=column.name,
            source_kind="heatmap",
            value_kind="categorical",
            observed_taxon_count=len(observed_values),
            missing_taxon_count=len(taxa)
            - len({taxon for taxon in taxa if taxon in column.values}),
            distinct_value_count=len(unique_values),
            minimum_numeric_value=None,
            maximum_numeric_value=None,
            example_values=unique_values[:4],
        )
    unique_values = sorted(set(numeric_values))
    return AnnotatedTraitTreeSummaryRow(
        surface=column.name,
        source_kind="heatmap",
        value_kind="numeric",
        observed_taxon_count=len(observed_values),
        missing_taxon_count=len(taxa)
        - len({taxon for taxon in taxa if taxon in column.values}),
        distinct_value_count=len(unique_values),
        minimum_numeric_value=min(numeric_values),
        maximum_numeric_value=max(numeric_values),
        example_values=[format(value, ".6g") for value in unique_values[:4]],
    )
