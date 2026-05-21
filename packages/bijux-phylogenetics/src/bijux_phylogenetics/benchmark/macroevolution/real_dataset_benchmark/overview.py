from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
)

from .contracts import RealDatasetMacroevolutionBenchmarkBundle
from .shared import PROVENANCE_CITATION, PROVENANCE_DOI


def write_overview(
    path: Path,
    *,
    dataset: CentralEuropeanSeashoreFloraDataset,
    bundle: RealDatasetMacroevolutionBenchmarkBundle,
) -> Path:
    lines = [
        "# Real-Dataset Macroevolution Benchmark Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- dataset label: `{dataset.label}`",
        f"- provenance: `{PROVENANCE_CITATION}`",
        f"- DOI: `{PROVENANCE_DOI}`",
        f"- benchmark bundle directory: `{bundle.output_root.name}`",
        "",
        "Generated outputs:",
        "",
        f"- review traits input: `{bundle.review_traits_path.relative_to(bundle.output_root.parent)}`",
        f"- summary ledger: `{bundle.summary_path.relative_to(bundle.output_root.parent)}`",
        f"- native model table: `{bundle.model_table_path.relative_to(bundle.output_root.parent)}`",
        f"- alignment review ledger: `{bundle.alignment_review_path.relative_to(bundle.output_root.parent)}`",
        f"- geiger parity ledger: `{bundle.parity_table_path.relative_to(bundle.output_root.parent)}`",
        f"- stored geiger reference ledger: `{bundle.geiger_reference_path.relative_to(bundle.output_root.parent)}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
