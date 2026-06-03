from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from bijux_phylogenetics.datasets import (
    RabiesCrossHostGeographyPanelDemoResult,
    load_rabies_cross_host_geography_panel_dataset,
    run_rabies_cross_host_geography_panel_demo,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.demo.inventory import (
    _write_package_artifact_inventory,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.demo.manifest import (
    _write_demo_package_manifest,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.demo.overview import (
    _build_flagship_answer_summary,
    _write_overview,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.demo.presentation import (
    _write_demo_overview_html,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.demo.reproducibility import (
    _write_package_reproducibility_checklist,
)


def build_stub_rabies_cross_host_geography_package(
    output_root: Path,
    monkeypatch,
):
    """Build one lightweight rabies study package through the real demo surface."""

    dataset = load_rabies_cross_host_geography_panel_dataset()
    workflow_root = output_root / "workflow"
    workflow_root.mkdir(parents=True, exist_ok=True)
    final_report_path = workflow_root / "rabies-cross-host-geography-report.html"
    workflow_summary_path = workflow_root / "workflow-summary.tsv"
    resource_observations_path = workflow_root / "resource-observations.tsv"
    final_manifest_path = workflow_root / "rabies-cross-host-geography.manifest.json"
    scientific_findings_path = workflow_root / "scientific-findings.tsv"
    alignment_path = workflow_root / "rabies-cross-host-geography-panel.aln"
    trimmed_alignment_path = (
        workflow_root / "rabies-cross-host-geography-panel.trimmed.aln"
    )
    bootstrap_summary_path = (
        workflow_root / "bootstrap-review" / "bootstrap-review.summary.tsv"
    )
    comparative_report_path = workflow_root / "comparative" / "comparative-report.html"
    biogeography_report_path = (
        workflow_root / "biogeography" / "biogeography-report.html"
    )
    biogeography_tree_figure_path = (
        workflow_root / "biogeography" / "biogeography-tree.svg"
    )
    biogeography_map_path = workflow_root / "biogeography" / "biogeography-map.svg"
    conclusion_stability_report_path = (
        workflow_root / "conclusion-stability" / "conclusion-stability-report.html"
    )
    for path, contents in (
        (final_report_path, "<html></html>\n"),
        (workflow_summary_path, "metric\tvalue\nsequence_count\t9\n"),
        (resource_observations_path, "metric\tvalue\nworkflow_runtime_seconds\t1\n"),
        (final_manifest_path, "{}\n"),
        (
            scientific_findings_path,
            "finding_id\tquestion\tclaim\tevidence\tcaution\tsource_artifact\n",
        ),
        (
            alignment_path,
            ">taxon_a\nACGT\n>taxon_b\nACGT\n",
        ),
        (
            trimmed_alignment_path,
            ">taxon_a\nACGT\n>taxon_b\nACGT\n",
        ),
        (bootstrap_summary_path, "metric\tvalue\n"),
        (comparative_report_path, "<html></html>\n"),
        (biogeography_report_path, "<html></html>\n"),
        (biogeography_tree_figure_path, "<svg><text>tree</text></svg>\n"),
        (biogeography_map_path, "<svg><text>map</text></svg>\n"),
        (conclusion_stability_report_path, "<html></html>\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf-8")

    fake_config = SimpleNamespace(
        config_path=dataset.workflow_config_path,
        workflow_prefix=dataset.workflow_prefix,
        alignment_mode="mafft_auto",
        trimming_mode="strictplus",
        trim_gap_threshold=0.2,
        bootstrap_consensus_threshold=0.5,
        bootstrap_robust_support_threshold=0.8,
        comparative_formula=dataset.comparative_formula,
        comparative_response=dataset.comparative_response,
        comparative_branch_length_floor=dataset.comparative_branch_length_floor,
        timeout_seconds=dataset.timeout_seconds,
        max_bootstrap_tree_count=dataset.max_bootstrap_tree_count,
        max_report_table_rows=dataset.max_report_table_rows,
        memory_warning_threshold_bytes=dataset.memory_warning_threshold_bytes,
    )
    fake_workflow_report = SimpleNamespace(config=fake_config)
    fake_workflow_bundle = SimpleNamespace(
        comparative_formula=dataset.comparative_formula,
        root_host="bat",
        root_region="north_asia",
        selected_model="TPM2u+F+G4",
        comparative_selected_model="brownian",
        comparative_pgls_lambda=1.0,
        comparative_pgls_r_squared=0.833944827574,
        scientific_finding_count=6,
        sequence_type="dna",
        output_root=workflow_root,
        workflow_summary_path=workflow_summary_path,
        resource_observations_path=resource_observations_path,
        tree_path=workflow_root / "rabies-cross-host-geography-panel.rooted.tree",
        rooting_report_path=workflow_root
        / "rabies-cross-host-geography-panel.rooting.tsv",
        model_table_path=workflow_root / "rabies-cross-host-geography-panel.model.tsv",
        support_table_path=workflow_root
        / "rabies-cross-host-geography-panel.support.tsv",
        manifest_path=workflow_root / "rabies-cross-host-geography-panel.manifest.json",
        log_path=workflow_root / "rabies-cross-host-geography-panel.log",
        clade_table_path=workflow_root / "clade-table.tsv",
        bootstrap_summary_path=bootstrap_summary_path,
        bootstrap_tree_comparison_summary_path=workflow_root
        / "bootstrap-review"
        / "rooted-tree-vs-bootstrap-consensus.summary.tsv",
        host_switch_summary_path=workflow_root / "host-switch-summary.tsv",
        comparative_report_path=comparative_report_path,
        biogeography_report_path=biogeography_report_path,
        conclusion_stability_report_path=conclusion_stability_report_path,
        final_report_path=final_report_path,
        final_manifest_path=final_manifest_path,
        scientific_findings_path=scientific_findings_path,
        bootstrap_tree_count=1000,
        bootstrap_topology_count=4,
        workflow_runtime_seconds=12.5,
        bootstrap_review_runtime_seconds=0.25,
        bootstrap_review_peak_memory_bytes=2048,
        budget_warning_count=0,
        host_switch_count=2,
        migration_event_count=4,
        comparative_branch_repair_count=0,
        conclusion_stable_count=6,
        conclusion_weak_count=3,
        conclusion_unstable_count=1,
    )

    for path in (
        fake_workflow_bundle.tree_path,
        fake_workflow_bundle.rooting_report_path,
        fake_workflow_bundle.model_table_path,
        fake_workflow_bundle.support_table_path,
        fake_workflow_bundle.manifest_path,
        fake_workflow_bundle.log_path,
        fake_workflow_bundle.clade_table_path,
        fake_workflow_bundle.bootstrap_tree_comparison_summary_path,
        fake_workflow_bundle.host_switch_summary_path,
        fake_workflow_bundle.biogeography_report_path,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path == fake_workflow_bundle.tree_path:
            path.write_text(
                "(bat_chile_rv108:0.1,fox_canada_rv241:0.2)root;\n", encoding="utf-8"
            )
        else:
            path.write_text("metric\tvalue\n", encoding="utf-8")

    monkeypatch.setattr(
        "bijux_phylogenetics.datasets.rabies_host_geography.run_rabies_cross_host_geography_panel_workflow",
        lambda *args, **kwargs: fake_workflow_report,
    )

    def _fake_write_bundle(*args, **kwargs):
        for path in (
            final_report_path,
            workflow_summary_path,
            resource_observations_path,
            final_manifest_path,
            scientific_findings_path,
            alignment_path,
            trimmed_alignment_path,
            bootstrap_summary_path,
            comparative_report_path,
            biogeography_report_path,
            biogeography_tree_figure_path,
            biogeography_map_path,
            conclusion_stability_report_path,
            fake_workflow_bundle.tree_path,
            fake_workflow_bundle.rooting_report_path,
            fake_workflow_bundle.model_table_path,
            fake_workflow_bundle.support_table_path,
            fake_workflow_bundle.manifest_path,
            fake_workflow_bundle.log_path,
            fake_workflow_bundle.clade_table_path,
            fake_workflow_bundle.bootstrap_tree_comparison_summary_path,
            fake_workflow_bundle.host_switch_summary_path,
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                if path == fake_workflow_bundle.tree_path:
                    path.write_text(
                        "(bat_chile_rv108:0.1,fox_canada_rv241:0.2)root;\n",
                        encoding="utf-8",
                    )
                elif path in {alignment_path, trimmed_alignment_path}:
                    path.write_text(
                        ">taxon_a\nACGT\n>taxon_b\nACGT\n", encoding="utf-8"
                    )
                elif path == biogeography_tree_figure_path:
                    path.write_text("<svg><text>tree</text></svg>\n", encoding="utf-8")
                elif path == biogeography_map_path:
                    path.write_text("<svg><text>map</text></svg>\n", encoding="utf-8")
                else:
                    path.write_text("metric\tvalue\n", encoding="utf-8")
        return fake_workflow_bundle

    monkeypatch.setattr(
        "bijux_phylogenetics.datasets.rabies_host_geography.write_rabies_cross_host_geography_panel_workflow_bundle",
        _fake_write_bundle,
    )
    return run_rabies_cross_host_geography_panel_demo(output_root)


def _stub_workflow_config(dataset) -> SimpleNamespace:
    return SimpleNamespace(
        config_path=dataset.workflow_config_path,
        workflow_prefix=dataset.workflow_prefix,
        alignment_mode="mafft_auto",
        trimming_mode="strictplus",
        trim_gap_threshold=0.2,
        bootstrap_consensus_threshold=0.5,
        bootstrap_robust_support_threshold=0.8,
        comparative_formula=dataset.comparative_formula,
        comparative_response=dataset.comparative_response,
        comparative_branch_length_floor=dataset.comparative_branch_length_floor,
        timeout_seconds=dataset.timeout_seconds,
        max_bootstrap_tree_count=dataset.max_bootstrap_tree_count,
        max_report_table_rows=dataset.max_report_table_rows,
        memory_warning_threshold_bytes=dataset.memory_warning_threshold_bytes,
    )


def refresh_stub_rabies_cross_host_geography_package(
    result: RabiesCrossHostGeographyPanelDemoResult,
) -> RabiesCrossHostGeographyPanelDemoResult:
    """Rebuild package-control artifacts after mutating a stub rabies package."""

    config = _stub_workflow_config(result.dataset)
    short_answer = _build_flagship_answer_summary(result.workflow_bundle)  # noqa: SLF001
    _write_overview(  # noqa: SLF001
        result.overview_path,
        dataset=result.dataset,
        workflow_bundle=result.workflow_bundle,
        config=config,
        short_answer=short_answer,
        artifact_inventory_path=result.artifact_inventory_path,
        reproducibility_checklist_path=result.reproducibility_checklist_path,
    )
    _write_demo_overview_html(  # noqa: SLF001
        result.overview_html_path,
        dataset=result.dataset,
        dataset_export=result.dataset_export,
        workflow_bundle=result.workflow_bundle,
        config=config,
        short_answer=short_answer,
        artifact_inventory_path=result.artifact_inventory_path,
        reproducibility_checklist_path=result.reproducibility_checklist_path,
    )
    artifact_inventory_path, artifact_inventory_rows = (
        _write_package_artifact_inventory(  # noqa: SLF001
            result.artifact_inventory_path,
            output_root=result.output_root,
            dataset_export=result.dataset_export,
            workflow_bundle=result.workflow_bundle,
            overview_path=result.overview_path,
            overview_html_path=result.overview_html_path,
        )
    )
    reproducibility_checklist_path, checklist_rows = (
        _write_package_reproducibility_checklist(  # noqa: SLF001
            result.reproducibility_checklist_path,
            workflow_bundle=result.workflow_bundle,
            inventory_rows=artifact_inventory_rows,
            artifact_inventory_path=artifact_inventory_path,
        )
    )
    package_manifest_path = _write_demo_package_manifest(  # noqa: SLF001
        result.package_manifest_path,
        dataset=result.dataset,
        dataset_export=result.dataset_export,
        workflow_bundle=result.workflow_bundle,
        config=config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        artifact_inventory_rows=artifact_inventory_rows,
        reproducibility_checklist_path=reproducibility_checklist_path,
        checklist_rows=checklist_rows,
    )
    return RabiesCrossHostGeographyPanelDemoResult(
        output_root=result.output_root,
        dataset=result.dataset,
        dataset_export=result.dataset_export,
        workflow_bundle=result.workflow_bundle,
        overview_path=result.overview_path,
        overview_html_path=result.overview_html_path,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
        package_manifest_path=package_manifest_path,
    )
