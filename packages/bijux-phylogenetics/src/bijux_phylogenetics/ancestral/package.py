from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.service import (
    render_ancestral_state_tree,
    write_ancestral_state_table,
)


@dataclass(slots=True)
class AncestralFigurePackageResult:
    output_dir: Path
    figure_path: Path
    node_table_path: Path
    uncertainty_table_path: Path
    legend_path: Path
    model_description_path: Path
    caption_path: Path
    manifest_path: Path


def build_ancestral_figure_package(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    out_dir: Path,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    layout: str = "phylogram",
) -> AncestralFigurePackageResult:
    """Build a publication-ready package for one ancestral-state reconstruction."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if reconstruction_kind == "continuous":
        report: ContinuousAncestralReport | DiscreteAncestralReport = (
            reconstruct_continuous_ancestral_states(
                tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                alpha=alpha,
            )
        )
    else:
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )

    figure_path = out_dir / "ancestral-figure.svg"
    node_table_path = out_dir / "node-states.tsv"
    uncertainty_table_path = out_dir / "uncertainty.tsv"
    legend_path = out_dir / "legend.md"
    model_description_path = out_dir / "model-description.md"
    caption_path = out_dir / "figure-caption.md"
    manifest_path = out_dir / "figure-manifest.json"

    render_ancestral_state_tree(tree_path, report, out_path=figure_path, layout=layout)
    write_ancestral_state_table(node_table_path, report)
    if isinstance(report, ContinuousAncestralReport):
        write_continuous_ancestral_uncertainty_table(uncertainty_table_path, report)
    else:
        _write_uncertainty_table(uncertainty_table_path, report)
    legend_path.write_text(_legend_markdown(report), encoding="utf-8")
    model_description_path.write_text(_model_description(report), encoding="utf-8")
    caption_path.write_text(_caption_markdown(report, layout=layout), encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "trait": trait,
                "reconstruction_kind": reconstruction_kind,
                "model": model,
                "layout": layout,
                "artifacts": {
                    "figure": str(figure_path),
                    "node_table": str(node_table_path),
                    "uncertainty_table": str(uncertainty_table_path),
                    "legend": str(legend_path),
                    "model_description": str(model_description_path),
                    "caption": str(caption_path),
                },
                "report": asdict(report),
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n",
        encoding="utf-8",
    )
    return AncestralFigurePackageResult(
        output_dir=out_dir,
        figure_path=figure_path,
        node_table_path=node_table_path,
        uncertainty_table_path=uncertainty_table_path,
        legend_path=legend_path,
        model_description_path=model_description_path,
        caption_path=caption_path,
        manifest_path=manifest_path,
    )


def _write_uncertainty_table(
    path: Path,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> Path:
    if isinstance(report, ContinuousAncestralReport):
        return write_ancestral_rows(
            path,
            columns=[
                "node",
                "descendant_taxa",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
                "interpretation",
            ],
            rows=[
                {
                    "node": estimate.node,
                    "descendant_taxa": ",".join(estimate.descendant_taxa),
                    "standard_error": str(estimate.standard_error),
                    "lower_95_interval": str(estimate.lower_95_interval),
                    "upper_95_interval": str(estimate.upper_95_interval),
                    "interpretation": estimate.interpretation,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
        )
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "most_likely_state",
            "confidence",
            "unstable",
            "interpretation",
        ],
        rows=[
            {
                "node": estimate.node,
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "most_likely_state": estimate.most_likely_state,
                "confidence": str(estimate.confidence),
                "unstable": str(estimate.unstable).lower(),
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def _legend_markdown(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> str:
    if isinstance(report, ContinuousAncestralReport):
        return (
            "# Legend\n\n"
            "- Internal labels: reconstructed continuous ancestral values.\n"
            "- Uncertainty table: standard error and 95% interval for each internal node.\n"
        )
    states = "\n".join(f"- `{state}`" for state in report.observed_states)
    return (
        "# Legend\n\n"
        "- Internal labels: reconstructed discrete ancestral states.\n"
        "- Observed states:\n"
        f"{states}\n"
        "- Uncertainty table: most-likely state confidence and instability flag per internal node.\n"
    )


def _model_description(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> str:
    if isinstance(report, ContinuousAncestralReport):
        return (
            "# Model Description\n\n"
            f"- Kind: continuous ancestral reconstruction\n"
            f"- Model: `{report.model}`\n"
            f"- Alpha: `{report.alpha}`\n"
            f"- Trait: `{report.trait}`\n"
        )
    return (
        "# Model Description\n\n"
        f"- Kind: discrete ancestral reconstruction\n"
        f"- Model: `{report.model}`\n"
        f"- Trait: `{report.trait}`\n"
        f"- Observed states: `{', '.join(report.observed_states)}`\n"
    )


def _caption_markdown(
    report: ContinuousAncestralReport | DiscreteAncestralReport, *, layout: str
) -> str:
    return (
        f"# Ancestral State Figure: {report.trait}\n\n"
        f"- Layout: `{layout}`\n"
        f"- Model: `{report.model}`\n"
        f"- Taxa analyzed: `{report.taxon_count}`\n"
        f"- Internal nodes summarized: `{sum(1 for estimate in report.estimates if not estimate.is_tip)}`\n"
    )
