from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from bijux_phylogenetics.datasets import (
    load_rabies_cross_host_geography_panel_dataset,
    run_rabies_cross_host_geography_panel_demo,
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
    bootstrap_summary_path = (
        workflow_root / "bootstrap-review" / "bootstrap-review.summary.tsv"
    )
    comparative_report_path = (
        workflow_root / "comparative" / "comparative-report.html"
    )
    biogeography_report_path = (
        workflow_root / "biogeography" / "biogeography-report.html"
    )
    conclusion_stability_report_path = (
        workflow_root
        / "conclusion-stability"
        / "conclusion-stability-report.html"
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
        (bootstrap_summary_path, "metric\tvalue\n"),
        (comparative_report_path, "<html></html>\n"),
        (biogeography_report_path, "<html></html>\n"),
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
        rooting_report_path=workflow_root / "rabies-cross-host-geography-panel.rooting.tsv",
        model_table_path=workflow_root / "rabies-cross-host-geography-panel.model.tsv",
        support_table_path=workflow_root / "rabies-cross-host-geography-panel.support.tsv",
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
            bootstrap_summary_path,
            comparative_report_path,
            biogeography_report_path,
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
                path.write_text("metric\tvalue\n", encoding="utf-8")
        return fake_workflow_bundle

    monkeypatch.setattr(
        "bijux_phylogenetics.datasets.rabies_host_geography.write_rabies_cross_host_geography_panel_workflow_bundle",
        _fake_write_bundle,
    )
    return run_rabies_cross_host_geography_panel_demo(output_root)
