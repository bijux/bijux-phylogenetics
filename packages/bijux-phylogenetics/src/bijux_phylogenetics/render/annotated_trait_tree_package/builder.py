from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.tree_figure_package import build_tree_figure_package

from .artifacts import (
    build_package_manifest,
    sha256,
    write_annotation_coverage_table,
    write_annotation_summary_table,
    write_package_manifest,
    write_package_reproducibility_manifest,
)
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
from .presentation import write_review_report
from .summaries import (
    build_coverage_row,
    build_heatmap_summary_row,
    build_label_summary_row,
    build_numeric_summary_row,
    build_string_summary_row,
)


def build_annotated_trait_tree_package(
    tree_path: Path,
    *,
    out_dir: Path,
    metadata_path: Path | None = None,
    traits_path: Path | None = None,
    taxon_column: str | None = None,
    label_column: str | None = None,
    categorical_column: str | None = None,
    continuous_column: str | None = None,
    metadata_strip_columns: list[str] | None = None,
    heatmap_columns: list[str] | None = None,
    layout: str = "phylogram",
    show_support_values: bool = True,
    title: str = "Bijux Annotated Trait Tree",
) -> AnnotatedTraitTreePackageResult:
    """Build a publication package for one annotated trait tree."""
    out_dir.mkdir(parents=True, exist_ok=True)
    coverage_path = out_dir / "annotation-coverage.tsv"
    summary_path = out_dir / "annotation-surface-summary.tsv"
    review_path = out_dir / "annotated-trait-tree-review.html"
    manifest_path = out_dir / "annotated-trait-tree-package-manifest.json"
    reproducibility_manifest_path = (
        out_dir / "annotated-trait-tree-reproducibility.manifest.json"
    )

    metadata_strip_columns = metadata_strip_columns or []
    heatmap_columns = heatmap_columns or []
    if (
        categorical_column is None
        and continuous_column is None
        and not metadata_strip_columns
        and not heatmap_columns
    ):
        raise ValueError(
            "annotated trait tree package requires at least one trait or metadata annotation surface"
        )

    tree = load_tree(tree_path)
    taxa = tree.tip_names
    metadata_table = (
        load_taxon_table(metadata_path, taxon_column=taxon_column)
        if metadata_path is not None
        else None
    )
    traits_table = (
        load_taxon_table(traits_path, taxon_column=taxon_column)
        if traits_path is not None
        else None
    )

    labels = build_full_label_map(
        taxa=taxa,
        metadata_table=metadata_table,
        label_column=label_column,
    )
    categorical_traits: dict[str, str] = {}
    if categorical_column is not None:
        categorical_traits = build_string_map(
            require_table(
                traits_table,
                path=traits_path,
                surface="a trait table for categorical trait rendering",
            ),
            categorical_column,
        )
    continuous_traits: dict[str, float] = {}
    if continuous_column is not None:
        continuous_traits = build_numeric_map(
            require_table(
                traits_table,
                path=traits_path,
                surface="a trait table for continuous trait rendering",
            ),
            continuous_column,
        )
    metadata_strips = (
        build_annotation_strips(
            require_table(
                metadata_table,
                path=metadata_path,
                surface="a metadata table for metadata strip rendering",
            ),
            metadata_strip_columns,
        )
        if metadata_strip_columns
        else []
    )
    heatmap_strips = (
        build_annotation_strips(
            require_table(
                traits_table,
                path=traits_path,
                surface="a trait table for heatmap rendering",
            ),
            heatmap_columns,
        )
        if heatmap_columns
        else []
    )

    figure_package = build_tree_figure_package(
        tree_path,
        out_dir=out_dir,
        title=title,
        labels=labels,
        layout=layout,
        show_support_values=show_support_values,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_strips,
    )

    coverage_rows: list[AnnotatedTraitTreeCoverageRow] = [
        build_coverage_row(
            surface="labels",
            source_kind="labels",
            required=True,
            taxa=taxa,
            observed_taxa=set(labels),
        )
    ]
    summary_rows: list[AnnotatedTraitTreeSummaryRow] = [
        build_label_summary_row(labels=labels, taxa=taxa)
    ]

    if categorical_column is not None:
        coverage_rows.append(
            build_coverage_row(
                surface=categorical_column,
                source_kind="categorical_trait",
                required=True,
                taxa=taxa,
                observed_taxa=set(categorical_traits),
            )
        )
        summary_rows.append(
            build_string_summary_row(
                surface=categorical_column,
                source_kind="categorical_trait",
                values=categorical_traits,
                taxa=taxa,
            )
        )
    if continuous_column is not None:
        coverage_rows.append(
            build_coverage_row(
                surface=continuous_column,
                source_kind="continuous_trait",
                required=True,
                taxa=taxa,
                observed_taxa=set(continuous_traits),
            )
        )
        summary_rows.append(
            build_numeric_summary_row(
                surface=continuous_column,
                source_kind="continuous_trait",
                values=continuous_traits,
                taxa=taxa,
            )
        )
    for strip in metadata_strips:
        coverage_rows.append(
            build_coverage_row(
                surface=strip.name,
                source_kind="metadata_strip",
                required=True,
                taxa=taxa,
                observed_taxa=set(strip.values),
            )
        )
        summary_rows.append(
            build_string_summary_row(
                surface=strip.name,
                source_kind="metadata_strip",
                values=strip.values,
                taxa=taxa,
            )
        )
    for strip in heatmap_strips:
        coverage_rows.append(
            build_coverage_row(
                surface=strip.name,
                source_kind="heatmap",
                required=True,
                taxa=taxa,
                observed_taxa=set(strip.values),
            )
        )
        summary_rows.append(build_heatmap_summary_row(column=strip, taxa=taxa))

    required_rows = [row for row in coverage_rows if row.required]
    complete_surface_count = sum(1 for row in required_rows if row.complete)
    missing_surface_count = len(required_rows) - complete_surface_count
    limitations: list[str] = []
    for row in required_rows:
        if not row.complete:
            limitations.append(
                f"{row.surface} is incomplete for publication review because one or more tree taxa are missing annotation values"
            )
    limitations.extend(figure_package.audit.limitations)
    limitations.extend(figure_package.legibility_audit.warnings)
    reviewer_summary = [
        f"{complete_surface_count} of {len(required_rows)} required annotation surfaces cover all tree taxa",
        f"figure package rendered {figure_package.render.visible_tip_count} visible tips with {len(figure_package.legend_entries)} explicit legend entries",
        "caption, legend, and legibility audits passed"
        if figure_package.caption_draft.caption_ready
        and figure_package.audit.legend_audit.complete
        and figure_package.legibility_audit.legible
        else "one or more base publication audits still require reviewer attention",
    ]
    audit = AnnotatedTraitTreePublicationAudit(
        publication_ready=complete_surface_count == len(required_rows)
        and figure_package.audit.legend_audit.complete
        and figure_package.caption_draft.caption_ready
        and figure_package.legibility_audit.legible,
        required_surface_count=len(required_rows),
        complete_surface_count=complete_surface_count,
        missing_surface_count=missing_surface_count,
        legend_entry_count=len(figure_package.legend_entries),
        caption_ready=figure_package.caption_draft.caption_ready,
        legible=figure_package.legibility_audit.legible,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )

    write_annotation_coverage_table(coverage_path, coverage_rows)
    write_annotation_summary_table(summary_path, summary_rows)

    manifest = build_package_manifest(
        title=title,
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        label_column=label_column,
        categorical_column=categorical_column,
        continuous_column=continuous_column,
        metadata_strip_columns=metadata_strip_columns,
        heatmap_columns=heatmap_columns,
        layout=layout,
        show_support_values=show_support_values,
        figure_package=figure_package,
        coverage_path=coverage_path,
        summary_path=summary_path,
        review_path=review_path,
        audit=audit,
        coverage_rows=coverage_rows,
        summary_rows=summary_rows,
    )
    reproducibility_manifest = write_package_reproducibility_manifest(
        reproducibility_manifest_path,
        title=title,
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        label_column=label_column,
        categorical_column=categorical_column,
        continuous_column=continuous_column,
        metadata_strip_columns=metadata_strip_columns,
        heatmap_columns=heatmap_columns,
        layout=layout,
        show_support_values=show_support_values,
        figure_package=figure_package,
        coverage_path=coverage_path,
        summary_path=summary_path,
    )
    manifest["reproducibility_manifest_path"] = str(reproducibility_manifest_path)
    manifest["reproducibility_manifest_checksum"] = sha256(
        reproducibility_manifest_path
    )
    manifest["reproducibility_manifest"] = reproducibility_manifest
    write_package_manifest(manifest_path, manifest)

    write_review_report(
        title=title,
        review_path=review_path,
        manifest_path=manifest_path,
        manifest=manifest,
        audit=audit,
        figure_package=figure_package,
        coverage_path=coverage_path,
        summary_path=summary_path,
        coverage_rows=coverage_rows,
        summary_rows=summary_rows,
    )

    return AnnotatedTraitTreePackageResult(
        output_dir=out_dir,
        figure_package=figure_package,
        coverage_path=coverage_path,
        summary_path=summary_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        coverage_rows=coverage_rows,
        summary_rows=summary_rows,
        audit=audit,
    )
