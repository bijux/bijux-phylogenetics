from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import RabiesCrossHostGeographyPanelWorkflowBundle


def _has_package_artifact(
    inventory_rows: list[dict[str, str]],
    relative_path: str,
) -> bool:
    return any(row["relative_path"] == relative_path for row in inventory_rows)


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
