from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import TaxonTable, load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.package import (
    TreeFigurePackageResult,
    build_tree_figure_package,
)
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.render.svg import AnnotationStrip
from bijux_phylogenetics.runtime.errors import MetadataJoinError


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_table(table: TaxonTable | None, *, path: Path | None, surface: str) -> TaxonTable:
    if table is None or path is None:
        raise MetadataJoinError(
            f"annotated trait tree package requires {surface} input table"
        )
    return table


def _build_full_label_map(
    *,
    taxa: list[str],
    metadata_table: TaxonTable | None,
    label_column: str | None,
) -> dict[str, str]:
    labels = {taxon: taxon for taxon in taxa}
    if label_column is None:
        return labels
    table = _require_table(
        metadata_table,
        path=None if metadata_table is None else metadata_table.path,
        surface="a metadata table for label rendering",
    )
    if label_column not in table.columns:
        raise MetadataJoinError(
            f"metadata table does not contain label column '{label_column}'"
        )
    for row in table.rows:
        taxon = row[table.taxon_column]
        if row[label_column]:
            labels[taxon] = row[label_column]
    return labels


def _build_string_map(table: TaxonTable, column: str) -> dict[str, str]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    return {row[table.taxon_column]: row[column] for row in table.rows if row[column]}


def _build_numeric_map(table: TaxonTable, column: str) -> dict[str, float]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    values: dict[str, float] = {}
    for row in table.rows:
        raw_value = row[column]
        if not raw_value:
            continue
        try:
            values[row[table.taxon_column]] = float(raw_value)
        except ValueError as error:
            raise MetadataJoinError(
                f"column '{column}' contains a non-numeric value for taxon '{row[table.taxon_column]}'"
            ) from error
    return values


def _build_annotation_strips(table: TaxonTable, columns: list[str]) -> list[AnnotationStrip]:
    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise MetadataJoinError(
            f"table does not contain columns: {', '.join(missing_columns)}"
        )
    return [
        AnnotationStrip(
            name=column,
            values={
                row[table.taxon_column]: row[column]
                for row in table.rows
                if row[column]
            },
        )
        for column in columns
    ]


def _coverage_row(
    *,
    surface: str,
    source_kind: str,
    required: bool,
    taxa: list[str],
    observed_taxa: set[str],
) -> AnnotatedTraitTreeCoverageRow:
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


def _label_summary_row(
    *,
    labels: dict[str, str],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
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


def _string_summary_row(
    *,
    surface: str,
    source_kind: str,
    values: dict[str, str],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    unique_values = sorted(set(values.values()))
    return AnnotatedTraitTreeSummaryRow(
        surface=surface,
        source_kind=source_kind,
        value_kind="categorical",
        observed_taxon_count=len(values),
        missing_taxon_count=len(taxa) - len({taxon for taxon in taxa if taxon in values}),
        distinct_value_count=len(unique_values),
        minimum_numeric_value=None,
        maximum_numeric_value=None,
        example_values=unique_values[:4],
    )


def _numeric_summary_row(
    *,
    surface: str,
    source_kind: str,
    values: dict[str, float],
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
    observed_values = list(values.values())
    return AnnotatedTraitTreeSummaryRow(
        surface=surface,
        source_kind=source_kind,
        value_kind="numeric",
        observed_taxon_count=len(observed_values),
        missing_taxon_count=len(taxa) - len({taxon for taxon in taxa if taxon in values}),
        distinct_value_count=len(set(observed_values)),
        minimum_numeric_value=min(observed_values) if observed_values else None,
        maximum_numeric_value=max(observed_values) if observed_values else None,
        example_values=[format(value, ".6g") for value in sorted(set(observed_values))[:4]],
    )


def _heatmap_summary_row(
    *,
    column: AnnotationStrip,
    taxa: list[str],
) -> AnnotatedTraitTreeSummaryRow:
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
            missing_taxon_count=len(taxa) - len(
                {taxon for taxon in taxa if taxon in column.values}
            ),
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
        missing_taxon_count=len(taxa) - len(
            {taxon for taxon in taxa if taxon in column.values}
        ),
        distinct_value_count=len(unique_values),
        minimum_numeric_value=min(numeric_values),
        maximum_numeric_value=max(numeric_values),
        example_values=[format(value, ".6g") for value in unique_values[:4]],
    )


def _coverage_lines(rows: list[AnnotatedTraitTreeCoverageRow]) -> str:
    return "\n".join(
        (
            f"{row.surface} [{row.source_kind}] complete={row.complete} "
            f"covered={row.covered_taxon_count}/{row.visible_taxon_count} "
            f"missing={','.join(row.missing_taxa) if row.missing_taxa else 'none'} "
            f"extra={','.join(row.extra_taxa) if row.extra_taxa else 'none'}"
        )
        for row in rows
    )


def _summary_lines(rows: list[AnnotatedTraitTreeSummaryRow]) -> str:
    rendered: list[str] = []
    for row in rows:
        numeric_range = (
            ""
            if row.minimum_numeric_value is None
            else f" range={format(row.minimum_numeric_value, '.6g')}..{format(row.maximum_numeric_value, '.6g')}"
        )
        examples = ",".join(row.example_values) if row.example_values else "none"
        rendered.append(
            f"{row.surface} [{row.source_kind}/{row.value_kind}] observed={row.observed_taxon_count} "
            f"missing={row.missing_taxon_count} distinct={row.distinct_value_count}{numeric_range} examples={examples}"
        )
    return "\n".join(rendered)


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

    labels = _build_full_label_map(
        taxa=taxa,
        metadata_table=metadata_table,
        label_column=label_column,
    )
    categorical_traits: dict[str, str] = {}
    if categorical_column is not None:
        categorical_traits = _build_string_map(
            _require_table(
                traits_table,
                path=traits_path,
                surface="a trait table for categorical trait rendering",
            ),
            categorical_column,
        )
    continuous_traits: dict[str, float] = {}
    if continuous_column is not None:
        continuous_traits = _build_numeric_map(
            _require_table(
                traits_table,
                path=traits_path,
                surface="a trait table for continuous trait rendering",
            ),
            continuous_column,
        )
    metadata_strips = (
        _build_annotation_strips(
            _require_table(
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
        _build_annotation_strips(
            _require_table(
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

    coverage_rows = [
        _coverage_row(
            surface="labels",
            source_kind="labels",
            required=True,
            taxa=taxa,
            observed_taxa=set(labels),
        )
    ]
    summary_rows = [_label_summary_row(labels=labels, taxa=taxa)]

    if categorical_column is not None:
        coverage_rows.append(
            _coverage_row(
                surface=categorical_column,
                source_kind="categorical_trait",
                required=True,
                taxa=taxa,
                observed_taxa=set(categorical_traits),
            )
        )
        summary_rows.append(
            _string_summary_row(
                surface=categorical_column,
                source_kind="categorical_trait",
                values=categorical_traits,
                taxa=taxa,
            )
        )
    if continuous_column is not None:
        coverage_rows.append(
            _coverage_row(
                surface=continuous_column,
                source_kind="continuous_trait",
                required=True,
                taxa=taxa,
                observed_taxa=set(continuous_traits),
            )
        )
        summary_rows.append(
            _numeric_summary_row(
                surface=continuous_column,
                source_kind="continuous_trait",
                values=continuous_traits,
                taxa=taxa,
            )
        )
    for strip in metadata_strips:
        coverage_rows.append(
            _coverage_row(
                surface=strip.name,
                source_kind="metadata_strip",
                required=True,
                taxa=taxa,
                observed_taxa=set(strip.values),
            )
        )
        summary_rows.append(
            _string_summary_row(
                surface=strip.name,
                source_kind="metadata_strip",
                values=strip.values,
                taxa=taxa,
            )
        )
    for strip in heatmap_strips:
        coverage_rows.append(
            _coverage_row(
                surface=strip.name,
                source_kind="heatmap",
                required=True,
                taxa=taxa,
                observed_taxa=set(strip.values),
            )
        )
        summary_rows.append(_heatmap_summary_row(column=strip, taxa=taxa))

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

    write_taxon_rows(
        coverage_path,
        columns=[
            "surface",
            "source_kind",
            "required",
            "visible_taxon_count",
            "observed_taxon_count",
            "covered_taxon_count",
            "complete",
            "missing_taxa",
            "extra_taxa",
        ],
        rows=[
            {
                "surface": row.surface,
                "source_kind": row.source_kind,
                "required": row.required,
                "visible_taxon_count": row.visible_taxon_count,
                "observed_taxon_count": row.observed_taxon_count,
                "covered_taxon_count": row.covered_taxon_count,
                "complete": row.complete,
                "missing_taxa": "|".join(row.missing_taxa),
                "extra_taxa": "|".join(row.extra_taxa),
            }
            for row in coverage_rows
        ],
    )
    write_taxon_rows(
        summary_path,
        columns=[
            "surface",
            "source_kind",
            "value_kind",
            "observed_taxon_count",
            "missing_taxon_count",
            "distinct_value_count",
            "minimum_numeric_value",
            "maximum_numeric_value",
            "example_values",
        ],
        rows=[
            {
                "surface": row.surface,
                "source_kind": row.source_kind,
                "value_kind": row.value_kind,
                "observed_taxon_count": row.observed_taxon_count,
                "missing_taxon_count": row.missing_taxon_count,
                "distinct_value_count": row.distinct_value_count,
                "minimum_numeric_value": row.minimum_numeric_value,
                "maximum_numeric_value": row.maximum_numeric_value,
                "example_values": "|".join(row.example_values),
            }
            for row in summary_rows
        ],
    )

    manifest = {
        "report_kind": "annotated_trait_tree_package",
        "title": title,
        "tree_path": str(tree_path),
        "metadata_path": None if metadata_path is None else str(metadata_path),
        "traits_path": None if traits_path is None else str(traits_path),
        "taxon_column": taxon_column,
        "label_column": label_column,
        "categorical_column": categorical_column,
        "continuous_column": continuous_column,
        "metadata_strip_columns": metadata_strip_columns,
        "heatmap_columns": heatmap_columns,
        "layout": layout,
        "show_support_values": show_support_values,
        "input_checksums": {
            str(path): _sha256(path)
            for path in (tree_path, metadata_path, traits_path)
            if path is not None
        },
        "figure_package_manifest_path": str(figure_package.manifest_path),
        "figure_package_manifest_checksum": _sha256(figure_package.manifest_path),
        "coverage_path": str(coverage_path),
        "coverage_checksum": _sha256(coverage_path),
        "summary_path": str(summary_path),
        "summary_checksum": _sha256(summary_path),
        "review_path": str(review_path),
        "audit": asdict(audit),
        "coverage_rows": [asdict(row) for row in coverage_rows],
        "summary_rows": [asdict(row) for row in summary_rows],
    }
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="annotated_trait_tree_package",
        input_files=[
            ("tree", tree_path),
            *([] if metadata_path is None else [("metadata", metadata_path)]),
            *([] if traits_path is None else [("traits", traits_path)]),
        ],
        generated_figures=[("annotated_trait_tree", figure_package.figure_path)],
        generated_tables=[
            ("tree_legend", figure_package.legend_path),
            ("tree_annotations", figure_package.annotations_path),
            ("annotation_coverage", coverage_path),
            ("annotation_surface_summary", summary_path),
        ],
        filters=None,
        model={
            "kind": "none",
            "name": None,
            "detail": "the annotated trait tree package overlays supplied labels and trait metadata without fitting a new statistical model",
        },
        settings={
            "title": title,
            "taxon_column": taxon_column,
            "label_column": label_column,
            "categorical_column": categorical_column,
            "continuous_column": continuous_column,
            "metadata_strip_columns": metadata_strip_columns,
            "heatmap_columns": heatmap_columns,
            "layout": layout,
            "show_support_values": show_support_values,
        },
        linked_artifacts=[
            ("tree_caption", figure_package.caption_path),
            ("tree_figure_manifest", figure_package.manifest_path),
            (
                "tree_figure_reproducibility_manifest",
                figure_package.reproducibility_manifest_path,
            ),
        ],
    )
    manifest["reproducibility_manifest_path"] = str(reproducibility_manifest_path)
    manifest["reproducibility_manifest_checksum"] = _sha256(
        reproducibility_manifest_path
    )
    manifest["reproducibility_manifest"] = reproducibility_manifest
    manifest_path.write_text(
        json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    write_html_report(
        title=title,
        out_path=review_path,
        embedded_json=manifest,
        summary_metrics=[
            ("publication_ready", audit.publication_ready),
            ("required_surface_count", audit.required_surface_count),
            ("complete_surface_count", audit.complete_surface_count),
            ("missing_surface_count", audit.missing_surface_count),
            ("legend_entry_count", audit.legend_entry_count),
            ("caption_ready", audit.caption_ready),
            ("legible", audit.legible),
        ],
        artifact_links=[
            ("figure", figure_package.figure_path.name, "annotated vector figure"),
            ("caption draft", figure_package.caption_path.name, None),
            ("legend ledger", figure_package.legend_path.name, None),
            ("tip annotations", figure_package.annotations_path.name, None),
            ("annotation coverage", coverage_path.name, None),
            ("surface summary", summary_path.name, None),
            ("package manifest", manifest_path.name, None),
        ],
        sections=[
            ("reviewer summary", "\n".join(f"- {line}" for line in audit.reviewer_summary)),
            ("publication limitations", "\n".join(f"- {line}" for line in audit.limitations) if audit.limitations else "none"),
            ("annotation coverage", _coverage_lines(coverage_rows)),
            ("annotation surface summary", _summary_lines(summary_rows)),
        ],
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
