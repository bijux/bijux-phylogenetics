# ruff: noqa: F401, F403, F405
from __future__ import annotations

import csv
from collections.abc import Callable
from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from Bio import Phylo

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    DiscreteAncestralTreeSetReport,
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.trees import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    write_clade_table,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSResult,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)
from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)
from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.ecology import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.fasta import load_permissive_fasta_records
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

from .bundle import write_rabies_cross_host_geography_panel_workflow_bundle
from .config import export_rabies_cross_host_geography_panel_dataset, load_rabies_cross_host_geography_panel_dataset
from .models import _FLAGSHIP_QUESTION
from .models import *
from .shared import _checksum, _html_list
from .workflow import run_rabies_cross_host_geography_panel_workflow


def _materialize_rabies_cross_host_geography_panel_demo(
    output_root: Path,
    *,
    config_path: Path | None,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    iqtree_seed: int | None,
    iqtree_threads: int | None,
    bootstrap_replicates: int | None,
    load_dataset: Callable[[Path | None], RabiesCrossHostGeographyPanelDataset],
    export_dataset: Callable[..., RabiesCrossHostGeographyPanelExportResult],
    run_workflow: Callable[..., RabiesCrossHostGeographyPanelWorkflowReport],
    write_workflow_bundle: Callable[
        [Path, RabiesCrossHostGeographyPanelWorkflowReport],
        RabiesCrossHostGeographyPanelWorkflowBundle,
    ],
) -> RabiesCrossHostGeographyPanelDemoResult:
    """Build one rabies demo package from the supplied workflow surfaces."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_dataset(config_path)
    dataset_export = export_dataset(
        output_root / "dataset",
        config_path=config_path,
    )
    with TemporaryDirectory(prefix="rabies-cross-host-geography-") as temporary_root:
        workflow_report = run_workflow(
            Path(temporary_root),
            config_path=config_path,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    short_answer = _build_flagship_answer_summary(workflow_bundle)
    artifact_inventory_path = output_root / "rabies-cross-host-geography-artifacts.tsv"
    reproducibility_checklist_path = (
        output_root / "rabies-cross-host-geography-reproducibility-checklist.tsv"
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        dataset=dataset,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    overview_html_path = _write_demo_overview_html(
        output_root / "rabies-cross-host-geography-overview.html",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    artifact_inventory_path, artifact_inventory_rows = _write_package_artifact_inventory(
        artifact_inventory_path,
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    reproducibility_checklist_path, checklist_rows = (
        _write_package_reproducibility_checklist(
            reproducibility_checklist_path,
            workflow_bundle=workflow_bundle,
            inventory_rows=artifact_inventory_rows,
            artifact_inventory_path=artifact_inventory_path,
        )
    )
    package_manifest_path = _write_demo_package_manifest(
        output_root / "rabies-cross-host-geography-package.manifest.json",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        artifact_inventory_rows=artifact_inventory_rows,
        reproducibility_checklist_path=reproducibility_checklist_path,
        checklist_rows=checklist_rows,
    )
    return RabiesCrossHostGeographyPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
        package_manifest_path=package_manifest_path,
    )


def run_rabies_cross_host_geography_panel_demo(
    output_root: Path,
    *,
    config_path: Path | None = None,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
) -> RabiesCrossHostGeographyPanelDemoResult:
    """Materialize the packaged integrated rabies dataset and rerun the full workflow."""
    return _materialize_rabies_cross_host_geography_panel_demo(
        output_root,
        config_path=config_path,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
        load_dataset=load_rabies_cross_host_geography_panel_dataset,
        export_dataset=export_rabies_cross_host_geography_panel_dataset,
        run_workflow=run_rabies_cross_host_geography_panel_workflow,
        write_workflow_bundle=write_rabies_cross_host_geography_panel_workflow_bundle,
    )



def _artifact_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if path.name.endswith(".manifest.json"):
        return "manifest"
    if suffix in {".html", ".htm"}:
        return "report"
    if suffix == ".md":
        return "markdown"
    if suffix == ".json":
        return "json"
    if suffix == ".tsv":
        return "table"
    if suffix == ".svg":
        return "figure"
    if suffix == ".log":
        return "log"
    if suffix in {".nwk", ".tree"}:
        return "tree"
    if suffix in {".aln", ".fasta"}:
        return "alignment"
    if suffix == ".csv":
        return "metadata"
    return "artifact"


def _workflow_artifact_section(relative_path: Path) -> str:
    parts = relative_path.parts
    if len(parts) >= 2 and parts[1] in {
        "bootstrap-review",
        "engine-artifacts",
        "biogeography",
        "comparative",
        "conclusion-stability",
    }:
        return parts[1]
    return "workflow"


def _relative_to_package_root(package_root: Path, path: Path) -> str:
    return path.relative_to(package_root).as_posix()


def _package_inventory_rows(
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    dataset_paths = [
        dataset_export.readme_path,
        dataset_export.workflow_config_path,
        dataset_export.sequences_path,
        dataset_export.metadata_path,
        dataset_export.centroids_path,
        dataset_export.accession_table_path,
    ]
    for path in dataset_paths:
        rows.append(
            {
                "section": "dataset",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    workflow_paths = sorted(
        path
        for path in workflow_bundle.output_root.rglob("*")
        if path.is_file()
    )
    for path in workflow_paths:
        rows.append(
            {
                "section": _workflow_artifact_section(
                    path.relative_to(workflow_bundle.output_root)
                ),
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    for path in (overview_path, overview_html_path):
        rows.append(
            {
                "section": "package",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    return rows


def _write_package_artifact_inventory(
    path: Path,
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> tuple[Path, list[dict[str, str]]]:
    rows = _package_inventory_rows(
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    return (
        write_taxon_rows(
            path,
            columns=["section", "kind", "relative_path", "sha256", "size_bytes"],
            rows=rows,
        ),
        rows,
    )


def _has_package_artifact(
    inventory_rows: list[dict[str, str]],
    relative_path: str,
) -> bool:
    return any(row["relative_path"] == relative_path for row in inventory_rows)


def _package_inventory_counts(
    inventory_rows: list[dict[str, str]],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in inventory_rows:
        section = row["section"]
        counts[section] = counts.get(section, 0) + 1
    return counts


def _write_package_reproducibility_checklist(
    path: Path,
    *,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    inventory_rows: list[dict[str, str]],
    artifact_inventory_path: Path,
) -> tuple[Path, list[dict[str, str]]]:
    rows = [
        {
            "section": "inputs",
            "check_id": "dataset-inputs-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "dataset/sequences.fasta",
                        "dataset/metadata.csv",
                        "dataset/region-centroids.csv",
                        "dataset/source-accessions.tsv",
                        "dataset/workflow-config.json",
                    )
                )
                else "blocked"
            ),
            "summary": "raw sequences, metadata, centroids, accession ledger, and workflow config are exported together",
            "evidence": "dataset export includes all package-level study inputs required to rerun the workflow",
            "artifact_path": "dataset/workflow-config.json",
        },
        {
            "section": "tree-inference",
            "check_id": "tree-inference-evidence-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        f"workflow/{workflow_bundle.tree_path.name}",
                        f"workflow/{workflow_bundle.model_table_path.name}",
                        f"workflow/{workflow_bundle.support_table_path.name}",
                        f"workflow/{workflow_bundle.manifest_path.name}",
                        f"workflow/{workflow_bundle.log_path.name}",
                    )
                )
                else "blocked"
            ),
            "summary": "tree inference outputs, manifest, and engine log are preserved",
            "evidence": "the package retains the rooted tree plus model, support, manifest, and log artifacts for the sequence-to-tree workflow",
            "artifact_path": f"workflow/{workflow_bundle.manifest_path.name}",
        },
        {
            "section": "uncertainty",
            "check_id": "bootstrap-uncertainty-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "workflow/bootstrap-review/bootstrap-review.summary.tsv",
                        "workflow/bootstrap-review/bootstrap-consensus.tree",
                        "workflow/bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv",
                        "workflow/bootstrap-review/rooted-tree-vs-bootstrap-consensus.report.html",
                    )
                )
                else "blocked"
            ),
            "summary": "bootstrap consensus and topology-conflict review artifacts are retained",
            "evidence": (
                "bootstrap review exported "
                f"{workflow_bundle.bootstrap_tree_count} trees over "
                f"{workflow_bundle.bootstrap_topology_count} rooted topologies"
            ),
            "artifact_path": "workflow/bootstrap-review/bootstrap-review.summary.tsv",
        },
        {
            "section": "analysis",
            "check_id": "downstream-analysis-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        f"workflow/{workflow_bundle.host_switch_summary_path.name}",
                        f"workflow/biogeography/{workflow_bundle.biogeography_report_path.name}",
                        f"workflow/comparative/{workflow_bundle.comparative_report_path.name}",
                        (
                            "workflow/conclusion-stability/"
                            f"{workflow_bundle.conclusion_stability_report_path.name}"
                        ),
                        f"workflow/{workflow_bundle.scientific_findings_path.name}",
                    )
                )
                else "blocked"
            ),
            "summary": "host-switching, biogeography, comparative, stability, and findings surfaces are preserved together",
            "evidence": (
                "the package retains one integrated downstream evidence chain from "
                "rooted tree to host/geography interpretation and stability review"
            ),
            "artifact_path": f"workflow/{workflow_bundle.scientific_findings_path.name}",
        },
        {
            "section": "package",
            "check_id": "package-navigation-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "overview.md",
                        "rabies-cross-host-geography-overview.html",
                    )
                )
                else "blocked"
            ),
            "summary": "reviewer overview surfaces are included with the package",
            "evidence": (
                "the package includes one markdown overview, one reviewer HTML overview, "
                f"and one artifact inventory at {artifact_inventory_path.name}"
            ),
            "artifact_path": "rabies-cross-host-geography-overview.html",
        },
        {
            "section": "limitations",
            "check_id": "interpretation-risks-surfaced",
            "status": (
                "risk"
                if (
                    workflow_bundle.budget_warning_count > 0
                    or workflow_bundle.conclusion_weak_count > 0
                    or workflow_bundle.conclusion_unstable_count > 0
                )
                else "pass"
            ),
            "summary": "the package records interpretation limits and stability caveats",
            "evidence": (
                "budget warnings="
                f"{workflow_bundle.budget_warning_count}; weak conclusions="
                f"{workflow_bundle.conclusion_weak_count}; unstable conclusions="
                f"{workflow_bundle.conclusion_unstable_count}"
            ),
            "artifact_path": (
                "workflow/conclusion-stability/"
                f"{workflow_bundle.conclusion_stability_report_path.name}"
            ),
        },
    ]
    return (
        write_taxon_rows(
            path,
            columns=[
                "section",
                "check_id",
                "status",
                "summary",
                "evidence",
                "artifact_path",
            ],
            rows=rows,
        ),
        rows,
    )


def _write_overview(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    lines = [
        "# Rabies Cross-Host Geography Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- workflow config: `{config.config_path.name}`",
        f"- biological question: {_FLAGSHIP_QUESTION}",
        f"- short answer: {short_answer}",
        f"- host workflow trait: `{dataset.host_trait}`",
        f"- geography workflow trait: `{dataset.geography_trait}`",
        f"- comparative formula: `{workflow_bundle.comparative_formula}`",
        "",
        "Generated outputs:",
        "",
        "- source accession ledger: `dataset/source-accessions.tsv`",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- resource observations: `{workflow_bundle.resource_observations_path.name}`",
        f"- clade table: `{workflow_bundle.clade_table_path.name}`",
        f"- bootstrap review: `bootstrap-review/{workflow_bundle.bootstrap_summary_path.name}`",
        (
            "- rooted-versus-consensus comparison: "
            f"`bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}`"
        ),
        f"- comparative report: `comparative/{workflow_bundle.comparative_report_path.name}`",
        (
            "- conclusion stability report: "
            f"`conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}`"
        ),
        f"- final report: `{workflow_bundle.final_report_path.name}`",
        f"- package artifact inventory: `{artifact_inventory_path.name}`",
        f"- package reproducibility checklist: `{reproducibility_checklist_path.name}`",
        "- package overview html: `rabies-cross-host-geography-overview.html`",
        "- package manifest: `rabies-cross-host-geography-package.manifest.json`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _build_flagship_answer_summary(
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
) -> str:
    return (
        "The rooted panel remains anchored in "
        f"{workflow_bundle.root_host} and {workflow_bundle.root_region}, and "
        "`host_group[canid]` shows a nominally supported positive longitude "
        "association under the selected comparative model, but the inference "
        "remains cautionary because the panel is intentionally compact."
    )


def _write_demo_overview_html(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Cross-Host Geography Package Overview</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #eef5ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 28px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.6; }",
            "    .panel { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 16px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    code { background: rgba(23,48,36,0.06); padding: 0 4px; border-radius: 4px; }",
            "    a { color: #23523b; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Cross-Host Geography Package</h1>",
            "  <p>This public flagship workflow starts from one real rabies nucleoprotein FASTA panel plus grouped host and geographic metadata, reruns the full owned sequence-to-tree workflow, and then connects the rooted tree to host-switching, geographic transition, migration, and comparative review surfaces.</p>",
            '  <section class="panel">',
            "    <h2>Biological Question</h2>",
            f"    <p>{escape(_FLAGSHIP_QUESTION)}</p>",
            '    <h2 style="margin-top: 16px;">Short Answer</h2>',
            f"    <p>{escape(short_answer)}</p>",
            "  </section>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">dataset id</span><span class="value">{escape(dataset.dataset_id)}</span></div>',
            f'    <div class="card"><span class="label">sequence count</span><span class="value">{dataset.sequence_count}</span></div>',
            f'    <div class="card"><span class="label">selected model</span><span class="value">{escape(workflow_bundle.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">comparative model</span><span class="value">{escape(workflow_bundle.comparative_selected_model)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Reproducibility Surface</h2>",
            _html_list(
                [
                    f'workflow config: <a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
                    f'source accession ledger: <a href="dataset/{dataset_export.accession_table_path.name}">dataset/{dataset_export.accession_table_path.name}</a>',
                    f'final workflow manifest: <a href="workflow/{workflow_bundle.final_manifest_path.name}">workflow/{workflow_bundle.final_manifest_path.name}</a>',
                    f'package artifact inventory: <a href="{artifact_inventory_path.name}">{artifact_inventory_path.name}</a>',
                    f'package reproducibility checklist: <a href="{reproducibility_checklist_path.name}">{reproducibility_checklist_path.name}</a>',
                    f'package manifest: <a href="{path.name.replace("-overview.html", "-package.manifest.json")}">{path.name.replace("-overview.html", "-package.manifest.json")}</a>',
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Key Outputs</h2>",
            _html_list(
                [
                    f'<a href="workflow/{workflow_bundle.final_report_path.name}">workflow/{workflow_bundle.final_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.workflow_summary_path.name}">workflow/{workflow_bundle.workflow_summary_path.name}</a>',
                    f'<a href="workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}">workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}</a>',
                    f'<a href="workflow/comparative/{workflow_bundle.comparative_report_path.name}">workflow/comparative/{workflow_bundle.comparative_report_path.name}</a>',
                    f'<a href="workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}">workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.scientific_findings_path.name}">workflow/{workflow_bundle.scientific_findings_path.name}</a>',
                    f'<a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _write_demo_package_manifest(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    artifact_inventory_rows: list[dict[str, str]],
    reproducibility_checklist_path: Path,
    checklist_rows: list[dict[str, str]],
) -> Path:
    inventory_counts = _package_inventory_counts(artifact_inventory_rows)
    blocked_check_count = len(
        [row for row in checklist_rows if row["status"] == "blocked"]
    )
    risk_check_count = len([row for row in checklist_rows if row["status"] == "risk"])
    payload = {
        "report_kind": "rabies_cross_host_geography_package",
        "dataset_id": dataset.dataset_id,
        "label": dataset.label,
        "biological_question": _FLAGSHIP_QUESTION,
        "short_answer": short_answer,
        "package_files": {
            "overview_markdown": {
                "path": "overview.md",
                "checksum": _checksum(path.parent / "overview.md"),
            },
            "overview_html": {
                "path": "rabies-cross-host-geography-overview.html",
                "checksum": _checksum(
                    path.parent / "rabies-cross-host-geography-overview.html"
                ),
            },
            "artifact_inventory": {
                "path": artifact_inventory_path.name,
                "checksum": _checksum(artifact_inventory_path),
                "artifact_count": len(artifact_inventory_rows),
                "section_counts": inventory_counts,
            },
            "reproducibility_checklist": {
                "path": reproducibility_checklist_path.name,
                "checksum": _checksum(reproducibility_checklist_path),
                "item_count": len(checklist_rows),
                "blocked_count": blocked_check_count,
                "risk_count": risk_check_count,
            },
        },
        "config": {
            "path": f"dataset/{dataset_export.workflow_config_path.name}",
            "checksum": _checksum(dataset_export.workflow_config_path),
            "workflow_prefix": config.workflow_prefix,
            "alignment_mode": config.alignment_mode,
            "trimming_mode": config.trimming_mode,
            "trim_gap_threshold": config.trim_gap_threshold,
            "bootstrap_consensus_threshold": config.bootstrap_consensus_threshold,
            "bootstrap_robust_support_threshold": (
                config.bootstrap_robust_support_threshold
            ),
            "comparative_formula": config.comparative_formula,
            "comparative_response": config.comparative_response,
            "comparative_branch_length_floor": (
                config.comparative_branch_length_floor
            ),
            "timeout_seconds": config.timeout_seconds,
            "max_bootstrap_tree_count": config.max_bootstrap_tree_count,
            "max_report_table_rows": config.max_report_table_rows,
            "memory_warning_threshold_bytes": config.memory_warning_threshold_bytes,
        },
        "dataset_files": {
            "readme": {
                "path": f"dataset/{dataset_export.readme_path.name}",
                "checksum": _checksum(dataset_export.readme_path),
            },
            "sequences": {
                "path": f"dataset/{dataset_export.sequences_path.name}",
                "checksum": _checksum(dataset_export.sequences_path),
            },
            "metadata": {
                "path": f"dataset/{dataset_export.metadata_path.name}",
                "checksum": _checksum(dataset_export.metadata_path),
            },
            "centroids": {
                "path": f"dataset/{dataset_export.centroids_path.name}",
                "checksum": _checksum(dataset_export.centroids_path),
            },
            "source_accessions": {
                "path": f"dataset/{dataset_export.accession_table_path.name}",
                "checksum": _checksum(dataset_export.accession_table_path),
            },
        },
        "workflow_files": {
            "final_report": {
                "path": f"workflow/{workflow_bundle.final_report_path.name}",
                "checksum": _checksum(workflow_bundle.final_report_path),
            },
            "workflow_log": {
                "path": f"workflow/{workflow_bundle.log_path.name}",
                "checksum": _checksum(workflow_bundle.log_path),
            },
            "workflow_summary": {
                "path": f"workflow/{workflow_bundle.workflow_summary_path.name}",
                "checksum": _checksum(workflow_bundle.workflow_summary_path),
            },
            "resource_observations": {
                "path": f"workflow/{workflow_bundle.resource_observations_path.name}",
                "checksum": _checksum(workflow_bundle.resource_observations_path),
            },
            "final_manifest": {
                "path": f"workflow/{workflow_bundle.final_manifest_path.name}",
                "checksum": _checksum(workflow_bundle.final_manifest_path),
            },
            "rooted_tree": {
                "path": f"workflow/{workflow_bundle.tree_path.name}",
                "checksum": _checksum(workflow_bundle.tree_path),
            },
            "rooting_report": {
                "path": f"workflow/{workflow_bundle.rooting_report_path.name}",
                "checksum": _checksum(workflow_bundle.rooting_report_path),
            },
            "model_table": {
                "path": f"workflow/{workflow_bundle.model_table_path.name}",
                "checksum": _checksum(workflow_bundle.model_table_path),
            },
            "support_table": {
                "path": f"workflow/{workflow_bundle.support_table_path.name}",
                "checksum": _checksum(workflow_bundle.support_table_path),
            },
            "bootstrap_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_summary_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.bootstrap_summary_path),
            },
            "bootstrap_tree_comparison_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_tree_comparison_summary_path.name}"
                ),
                "checksum": _checksum(
                    workflow_bundle.bootstrap_tree_comparison_summary_path
                ),
            },
            "host_switch_summary": {
                "path": f"workflow/{workflow_bundle.host_switch_summary_path.name}",
                "checksum": _checksum(workflow_bundle.host_switch_summary_path),
            },
            "biogeography_report": {
                "path": (
                    "workflow/biogeography/"
                    f"{workflow_bundle.biogeography_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.biogeography_report_path),
            },
            "comparative_report": {
                "path": (
                    "workflow/comparative/"
                    f"{workflow_bundle.comparative_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.comparative_report_path),
            },
            "conclusion_stability_report": {
                "path": (
                    "workflow/conclusion-stability/"
                    f"{workflow_bundle.conclusion_stability_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.conclusion_stability_report_path),
            },
            "scientific_findings": {
                "path": f"workflow/{workflow_bundle.scientific_findings_path.name}",
                "checksum": _checksum(workflow_bundle.scientific_findings_path),
            },
        },
        "metrics": {
            "sequence_count": dataset.sequence_count,
            "selected_model": workflow_bundle.selected_model,
            "root_host": workflow_bundle.root_host,
            "root_region": workflow_bundle.root_region,
            "bootstrap_tree_count": workflow_bundle.bootstrap_tree_count,
            "workflow_runtime_seconds": workflow_bundle.workflow_runtime_seconds,
            "bootstrap_review_runtime_seconds": (
                workflow_bundle.bootstrap_review_runtime_seconds
            ),
            "bootstrap_review_peak_memory_bytes": (
                workflow_bundle.bootstrap_review_peak_memory_bytes
            ),
            "budget_warning_count": workflow_bundle.budget_warning_count,
            "host_switch_count": workflow_bundle.host_switch_count,
            "migration_event_count": workflow_bundle.migration_event_count,
            "comparative_selected_model": workflow_bundle.comparative_selected_model,
            "comparative_pgls_lambda": workflow_bundle.comparative_pgls_lambda,
            "comparative_pgls_r_squared": workflow_bundle.comparative_pgls_r_squared,
            "comparative_branch_repair_count": (
                workflow_bundle.comparative_branch_repair_count
            ),
            "conclusion_stable_count": workflow_bundle.conclusion_stable_count,
            "conclusion_weak_count": workflow_bundle.conclusion_weak_count,
            "conclusion_unstable_count": workflow_bundle.conclusion_unstable_count,
            "scientific_finding_count": workflow_bundle.scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
