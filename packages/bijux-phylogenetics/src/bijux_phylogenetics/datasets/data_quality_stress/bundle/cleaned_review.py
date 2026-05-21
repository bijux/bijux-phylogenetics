from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import CatarrhineDataQualityStressPanelWorkflowReport
from .shared import _substantive_alignment_warnings


def _write_cleaned_linkage_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    alignment_taxa = {record.identifier for record in report.cleaned_alignment_records}
    tree_taxa = set(report.cleaned_linkage.usable_taxa)
    trait_taxa = set(report.cleaned_linkage.usable_taxa)
    taxa = sorted(alignment_taxa | tree_taxa | trait_taxa)
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "present_in_tree",
            "present_in_alignment",
            "present_in_traits",
        ],
        rows=[
            {
                "taxon": taxon,
                "present_in_tree": str(taxon in tree_taxa).lower(),
                "present_in_alignment": str(taxon in alignment_taxa).lower(),
                "present_in_traits": str(taxon in trait_taxa).lower(),
            }
            for taxon in taxa
        ],
    )


def _write_cleaned_validation_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    alignment_warning_count = len(
        _substantive_alignment_warnings(report.cleaned_alignment_validation.warnings)
    )
    rows = [
        {
            "surface": "alignment",
            "status": "pass" if alignment_warning_count == 0 else "warning",
            "warning_count": str(alignment_warning_count),
            "detail": (
                f"{report.cleaned_alignment_validation.summary.sequence_count} sequences "
                "remain in the cleaned alignment"
            ),
        },
        {
            "surface": "tree",
            "status": "pass"
            if report.cleaned_tree_validation.biologically_safe
            else "warning",
            "warning_count": str(len(report.cleaned_tree_validation.warnings)),
            "detail": report.cleaned_tree_validation.validity_decision,
        },
        {
            "surface": "traits",
            "status": "pass",
            "warning_count": "0",
            "detail": (
                f"{report.cleaned_trait_validation.row_count} cleaned trait rows keep "
                "all required comparative fields populated"
            ),
        },
        {
            "surface": "linkage",
            "status": "pass",
            "warning_count": "0",
            "detail": (
                f"{report.cleaned_linkage.linked_taxa} taxa overlap exactly across the "
                "cleaned tree and trait table"
            ),
        },
    ]
    return write_taxon_rows(
        path,
        columns=["surface", "status", "warning_count", "detail"],
        rows=rows,
    )
