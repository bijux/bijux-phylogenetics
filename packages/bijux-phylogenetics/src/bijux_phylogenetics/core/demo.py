from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.reports.service import render_dataset_report, render_phylo_inputs_report, render_tree_report


@dataclass(slots=True)
class DemoRunResult:
    output_root: Path
    input_root: Path
    report_root: Path
    tree_report: Path
    dataset_report: Path
    phylo_inputs_report: Path
    comparison_report: Path
    evidence_bundle: Path
    capability_summary: Path


def _fixture_root() -> Path:
    return Path(__file__).resolve().parents[3] / "tests" / "fixtures"


def _copy_demo_inputs(destination: Path) -> dict[str, Path]:
    source_root = _fixture_root()
    destination.mkdir(parents=True, exist_ok=True)
    selected = {
        "tree": source_root / "trees" / "example_tree.nwk",
        "alt_tree": source_root / "trees" / "example_tree_alt.nwk",
        "alignment": source_root / "alignments" / "example_alignment.fasta",
        "metadata": source_root / "metadata" / "example_metadata.tsv",
        "traits": source_root / "metadata" / "example_traits.tsv",
    }
    copied: dict[str, Path] = {}
    for name, source in selected.items():
        target = destination / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied[name] = target
    return copied


def _write_capability_summary(path: Path, result: DemoRunResult) -> Path:
    lines = [
        "# Bijux Phylogenetics Capability Demo",
        "",
        "This workflow demonstrates the v0.1 reporting, comparison, and evidence surfaces.",
        "",
        "## Artifacts",
        "",
        f"- `inputs/`: `{result.input_root}`",
        f"- `reports/`: `{result.report_root}`",
        f"- `tree report`: `{result.tree_report}`",
        f"- `dataset report`: `{result.dataset_report}`",
        f"- `phylo inputs report`: `{result.phylo_inputs_report}`",
        f"- `comparison report`: `{result.comparison_report}`",
        f"- `evidence bundle`: `{result.evidence_bundle}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_capability_demo(output_root: Path) -> DemoRunResult:
    """Generate a public capability demo from the repository sample inputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    input_root = output_root / "inputs"
    report_root = output_root / "reports"
    report_root.mkdir(parents=True, exist_ok=True)
    inputs = _copy_demo_inputs(input_root)

    tree_report = render_tree_report(tree_path=inputs["tree"], out_path=report_root / "tree-report.html").output_path
    dataset_report = render_dataset_report(
        tree_path=inputs["tree"],
        metadata_path=inputs["metadata"],
        traits_path=inputs["traits"],
        out_path=report_root / "dataset-report.html",
    ).output_path
    phylo_inputs_report = render_phylo_inputs_report(
        tree_path=inputs["tree"],
        alignment_path=inputs["alignment"],
        out_path=report_root / "phylo-inputs-report.html",
    ).output_path
    comparison_report = build_tree_comparison_report(
        inputs["tree"],
        inputs["alt_tree"],
        out_path=report_root / "comparison-report.html",
    ).output_path
    evidence_bundle = bundle_directory([input_root], [report_root], output_root / "evidence-pack").output_root
    capability_summary = _write_capability_summary(
        output_root / "capability-summary.md",
        DemoRunResult(
            output_root=output_root,
            input_root=input_root,
            report_root=report_root,
            tree_report=tree_report,
            dataset_report=dataset_report,
            phylo_inputs_report=phylo_inputs_report,
            comparison_report=comparison_report,
            evidence_bundle=evidence_bundle,
            capability_summary=output_root / "capability-summary.md",
        ),
    )
    return DemoRunResult(
        output_root=output_root,
        input_root=input_root,
        report_root=report_root,
        tree_report=tree_report,
        dataset_report=dataset_report,
        phylo_inputs_report=phylo_inputs_report,
        comparison_report=comparison_report,
        evidence_bundle=evidence_bundle,
        capability_summary=capability_summary,
    )
