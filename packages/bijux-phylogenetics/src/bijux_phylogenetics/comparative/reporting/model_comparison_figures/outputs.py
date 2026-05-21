from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import (
    ComparativeModelFigureCaptionDraft,
    ComparativeModelFigureLegendEntry,
)


def write_model_figure_table(path: Path, rows: list[object]) -> Path:
    return write_taxon_rows(
        path,
        columns=list(asdict(rows[0]).keys()),
        rows=[asdict(row) for row in rows],
    )


def write_model_figure_legend_table(
    path: Path, entries: list[ComparativeModelFigureLegendEntry]
) -> Path:
    return write_taxon_rows(
        path,
        columns=["surface", "label", "swatch", "detail"],
        rows=[
            {
                "surface": entry.surface,
                "label": entry.label,
                "swatch": entry.swatch,
                "detail": entry.detail,
            }
            for entry in entries
        ],
    )


def write_model_figure_caption(
    path: Path, draft: ComparativeModelFigureCaptionDraft
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.criteria_sentence,
                draft.likelihood_sentence,
                draft.parameter_sentence,
                draft.fit_sentence,
                draft.limitation_sentence,
                "",
                f"caption_ready: {'true' if draft.caption_ready else 'false'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path
