from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
    export_central_european_seashore_flora_dataset,
    load_central_european_seashore_flora_dataset,
)
from .real_dataset_benchmark.artifact_outputs import (
    write_geiger_real_dataset_reference_payload_table,
    write_real_dataset_macroevolution_alignment_review_table,
    write_real_dataset_macroevolution_model_table,
    write_real_dataset_macroevolution_parity_table,
    write_real_dataset_macroevolution_summary_table,
)
from .real_dataset_benchmark.contracts import (
    RealDatasetMacroevolutionBenchmarkBundle,
    RealDatasetMacroevolutionBenchmarkDemoResult,
    RealDatasetMacroevolutionBenchmarkReport,
)
from .real_dataset_benchmark.review_input import (
    write_alignment_review_traits_table as _write_alignment_review_traits_table,
)
from .real_dataset_benchmark.report_assembly import (
    build_report as _build_report,
)
from .real_dataset_benchmark.shared import (
    PROVENANCE_CITATION as _PROVENANCE_CITATION,
    PROVENANCE_DOI as _PROVENANCE_DOI,
)


def benchmark_real_dataset_macroevolution() -> RealDatasetMacroevolutionBenchmarkReport:
    """Benchmark continuous and discrete comparative fits on a real published dataset."""
    dataset = load_central_european_seashore_flora_dataset()
    with TemporaryDirectory(prefix="real-dataset-macroevolution-") as temporary_root:
        review_traits_path = _write_alignment_review_traits_table(
            Path(temporary_root) / "alignment-review-traits.csv",
            dataset,
        )
        return _build_report(dataset, review_traits_path)


def write_real_dataset_macroevolution_bundle(
    output_root: Path,
) -> RealDatasetMacroevolutionBenchmarkBundle:
    """Write the benchmark ledgers and review input for the real dataset benchmark."""
    if output_root.exists():
        for path in sorted(output_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_central_european_seashore_flora_dataset()
    review_traits_path = _write_alignment_review_traits_table(
        output_root / "alignment-review-traits.csv",
        dataset,
    )
    report = _build_report(dataset, review_traits_path)
    summary_path = write_real_dataset_macroevolution_summary_table(
        output_root / "benchmark-summary.tsv",
        report,
    )
    model_table_path = write_real_dataset_macroevolution_model_table(
        output_root / "model-table.tsv",
        report,
    )
    alignment_review_path = write_real_dataset_macroevolution_alignment_review_table(
        output_root / "alignment-review.tsv",
        report,
    )
    parity_table_path = write_real_dataset_macroevolution_parity_table(
        output_root / "geiger-parity.tsv",
        report,
    )
    geiger_reference_path = write_geiger_real_dataset_reference_payload_table(
        output_root / "geiger-reference.tsv",
        report,
    )
    return RealDatasetMacroevolutionBenchmarkBundle(
        output_root=output_root,
        review_traits_path=review_traits_path,
        summary_path=summary_path,
        model_table_path=model_table_path,
        alignment_review_path=alignment_review_path,
        parity_table_path=parity_table_path,
        geiger_reference_path=geiger_reference_path,
    )


def run_real_dataset_macroevolution_benchmark_demo(
    destination: Path,
) -> RealDatasetMacroevolutionBenchmarkDemoResult:
    """Materialize the public plant dataset plus the governed benchmark bundle."""
    dataset = load_central_european_seashore_flora_dataset()
    destination.mkdir(parents=True, exist_ok=True)
    dataset_export = export_central_european_seashore_flora_dataset(
        destination / "dataset"
    )
    benchmark_bundle = write_real_dataset_macroevolution_bundle(
        destination / "benchmark"
    )
    overview_path = _write_overview(
        destination / "README.md",
        dataset=dataset,
        bundle=benchmark_bundle,
    )
    return RealDatasetMacroevolutionBenchmarkDemoResult(
        output_root=destination,
        dataset=dataset,
        dataset_export=dataset_export,
        benchmark_bundle=benchmark_bundle,
        overview_path=overview_path,
    )


def _write_overview(
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
        f"- provenance: `{_PROVENANCE_CITATION}`",
        f"- DOI: `{_PROVENANCE_DOI}`",
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
