from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.tree_figure_package import TreeFigurePackageResult

from .contracts import (
    AnnotatedTraitTreeCoverageRow,
    AnnotatedTraitTreePublicationAudit,
    AnnotatedTraitTreeSummaryRow,
)


def build_coverage_lines(rows: list[AnnotatedTraitTreeCoverageRow]) -> str:
    """Render reviewer-facing coverage ledger lines for the HTML report."""
    return "\n".join(
        (
            f"{row.surface} [{row.source_kind}] complete={row.complete} "
            f"covered={row.covered_taxon_count}/{row.visible_taxon_count} "
            f"missing={','.join(row.missing_taxa) if row.missing_taxa else 'none'} "
            f"extra={','.join(row.extra_taxa) if row.extra_taxa else 'none'}"
        )
        for row in rows
    )


def build_summary_lines(rows: list[AnnotatedTraitTreeSummaryRow]) -> str:
    """Render reviewer-facing surface-summary lines for the HTML report."""
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


def write_review_report(
    *,
    title: str,
    review_path: Path,
    manifest_path: Path,
    manifest: dict[str, Any],
    audit: AnnotatedTraitTreePublicationAudit,
    figure_package: TreeFigurePackageResult,
    coverage_path: Path,
    summary_path: Path,
    coverage_rows: list[AnnotatedTraitTreeCoverageRow],
    summary_rows: list[AnnotatedTraitTreeSummaryRow],
) -> None:
    """Write the reviewer-facing HTML report for one annotated trait tree package."""
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
            (
                "reviewer summary",
                "\n".join(f"- {line}" for line in audit.reviewer_summary),
            ),
            (
                "publication limitations",
                "\n".join(f"- {line}" for line in audit.limitations)
                if audit.limitations
                else "none",
            ),
            ("annotation coverage", build_coverage_lines(coverage_rows)),
            ("annotation surface summary", build_summary_lines(summary_rows)),
        ],
    )
