from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.datasets import (
    RabiesCrossHostGeographyPanelDemoResult,
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowBundle,
    export_rabies_cross_host_geography_panel_dataset,
    load_rabies_cross_host_geography_panel_dataset,
    run_rabies_cross_host_geography_panel_demo,
    run_rabies_cross_host_geography_panel_workflow,
    write_rabies_cross_host_geography_panel_workflow_bundle,
)
import bijux_phylogenetics.datasets.rabies_host_geography as rabies_host_geography_api

from .support.external_engines import require_alignment_engine_executables
from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)

pytestmark = [
    pytest.mark.real_local,
    pytest.mark.evaluation,
    pytest.mark.scientific_validation,
]


def _stable_generated_outputs(bundle: object) -> dict[str, Path]:
    workflow_bundle = bundle
    excluded = {
        "rabies-cross-host-geography-panel.log",
        "rabies-cross-host-geography-panel.manifest.json",
        "biogeography/biogeography-report.manifest.json",
        "bootstrap-review/bootstrap-review.distance-matrix.tsv",
        "resource-observations.tsv",
    }
    outputs: dict[str, Path] = {}
    for path in workflow_bundle.output_root.rglob("*"):
        if not path.is_file():
            continue
        relative_name = str(path.relative_to(workflow_bundle.output_root))
        if relative_name.startswith("engine-artifacts/"):
            continue
        if relative_name in excluded:
            continue
        outputs[relative_name] = path
    return outputs


def test_load_rabies_cross_host_geography_panel_dataset_exposes_packaged_surface() -> (
    None
):
    dataset = load_rabies_cross_host_geography_panel_dataset()
    assert dataset.dataset_id == "rabies_cross_host_geography_panel"
    assert dataset.label == "Rabies cross-host geography panel"
    assert dataset.sequence_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_prefix == "rabies-cross-host-geography-panel"
    assert dataset.host_trait == "host_group"
    assert dataset.geography_trait == "region_group"
    assert dataset.host_model == "ard"
    assert dataset.geography_model == "ard"
    assert dataset.iqtree_seed == 1
    assert dataset.iqtree_threads == 1
    assert dataset.bootstrap_replicates == 1000
    assert dataset.timeout_seconds == 300.0
    assert dataset.max_bootstrap_tree_count == 1500
    assert dataset.max_report_table_rows == 25
    assert dataset.memory_warning_threshold_bytes == 67108864
    assert dataset.outgroup_taxa == ("bat_chile_rv108",)
    assert dataset.observed_host_group_count == 3
    assert dataset.observed_region_group_count == 5
    assert dataset.sequences_path.is_file()
    assert dataset.metadata_path.is_file()
    assert dataset.centroids_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert dataset.workflow_config_path.is_file()
    assert dataset.clade_metadata_columns == (
        "host_species",
        "host_group",
        "country",
        "region_group",
    )
    assert dataset.comparative_formula == "region_longitude ~ host_group"
    assert dataset.comparative_response == "region_longitude"
    assert dataset.comparative_branch_length_floor == 1e-6
    assert "MG458305" in dataset.source_accessions


def test_rabies_cross_host_geography_reference_bootstrap_summary_omits_volatile_metrics() -> (
    None
):
    dataset = load_rabies_cross_host_geography_panel_dataset()
    with (
        dataset.reference_output_root
        / "bootstrap-review"
        / "bootstrap-review.summary.tsv"
    ).open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert "runtime_seconds" not in row
    assert "peak_memory_bytes" not in row
    assert "skipped_malformed_tree_count" not in row


@pytest.mark.slow
def test_rabies_cross_host_geography_panel_workflow_bundle_writes_bootstrap_consensus_comparison(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_rabies_cross_host_geography_panel_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    assert bundle.bootstrap_tree_comparison_summary_path.is_file()
    assert bundle.bootstrap_tree_comparison_table_path.is_file()
    assert bundle.bootstrap_tree_comparison_report_path.is_file()
    with bundle.bootstrap_tree_comparison_summary_path.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    report_html = bundle.bootstrap_tree_comparison_report_path.read_text(
        encoding="utf-8"
    )
    assert int(row["rooted_rf_distance"]) == (
        bundle.bootstrap_consensus_rooted_rf_distance
    )
    assert row["same_unrooted_topology"] == (
        "true" if bundle.bootstrap_consensus_same_unrooted_topology else "false"
    )
    assert int(row["high_support_conflict_count"]) == (
        bundle.bootstrap_consensus_high_support_conflict_count
    )
    assert str(bundle.output_root) not in report_html


@pytest.mark.slow
def test_rabies_cross_host_geography_panel_workflow_bundle_writes_findings_and_config_audit(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_rabies_cross_host_geography_panel_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    assert bundle.config_audit_path.is_file()
    assert bundle.resolved_config_path.is_file()
    assert bundle.scientific_findings_path.is_file()
    with bundle.config_audit_path.open("r", encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    with bundle.scientific_findings_path.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        finding_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["status"] == "pass" for row in audit_rows)
    assert len(audit_rows) == bundle.config_check_count
    assert len(finding_rows) == bundle.scientific_finding_count
    assert {row["finding_id"] for row in finding_rows} >= {
        "bootstrap_consensus",
        "comparative_longitude",
    }


@pytest.mark.slow
def test_write_rabies_cross_host_geography_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    report = run_rabies_cross_host_geography_panel_workflow(
        tmp_path / "workflow-run",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = _stable_generated_outputs(bundle)
    expected_files = {
        str(path.relative_to(expected_root))
        for path in expected_root.rglob("*")
        if path.is_file()
    }
    assert expected_files == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_rabies_cross_host_geography_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    executables = require_alignment_engine_executables()
    result = run_rabies_cross_host_geography_panel_demo(
        tmp_path / "demo",
        mafft_executable=executables["mafft"],
        trimal_executable=executables["trimal"],
        iqtree_executable=executables["iqtree2"],
    )
    assert result.dataset.sequence_count == 9
    assert result.dataset_export.sequences_path.is_file()
    assert result.dataset_export.metadata_path.is_file()
    assert result.dataset_export.centroids_path.is_file()
    assert result.dataset_export.accession_table_path.is_file()
    assert result.dataset_export.workflow_config_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.resource_observations_path.is_file()
    assert result.workflow_bundle.config_audit_path.is_file()
    assert result.workflow_bundle.resolved_config_path.is_file()
    assert result.workflow_bundle.clade_table_path.is_file()
    assert result.workflow_bundle.bootstrap_summary_path.is_file()
    assert result.workflow_bundle.bootstrap_tree_comparison_summary_path.is_file()
    assert result.workflow_bundle.tree_path.is_file()
    assert result.workflow_bundle.host_switch_summary_path.is_file()
    assert result.workflow_bundle.biogeography_report_path.is_file()
    assert result.workflow_bundle.comparative_report_path.is_file()
    assert result.workflow_bundle.conclusion_stability_summary_path.is_file()
    assert result.workflow_bundle.key_clade_stability_path.is_file()
    assert result.workflow_bundle.support_value_stability_path.is_file()
    assert result.workflow_bundle.ancestral_state_stability_path.is_file()
    assert result.workflow_bundle.comparative_coefficient_stability_path.is_file()
    assert result.workflow_bundle.conclusion_stability_report_path.is_file()
    assert result.workflow_bundle.scientific_findings_path.is_file()
    assert result.workflow_bundle.final_report_path.is_file()
    assert result.overview_path.is_file()
    assert result.overview_html_path.is_file()
    assert result.artifact_inventory_path.is_file()
    assert result.reproducibility_checklist_path.is_file()
    assert result.package_manifest_path.is_file()
    overview = result.overview_path.read_text(encoding="utf-8")
    overview_html = result.overview_html_path.read_text(encoding="utf-8")
    package_manifest = json.loads(
        result.package_manifest_path.read_text(encoding="utf-8")
    )
    inventory_rows = list(
        csv.DictReader(
            result.artifact_inventory_path.open("r", encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    checklist_rows = list(
        csv.DictReader(
            result.reproducibility_checklist_path.open(
                "r", encoding="utf-8", newline=""
            ),
            delimiter="\t",
        )
    )
    assert "final report" in overview
    assert "comparative report" in overview
    assert "rooted-versus-consensus comparison" in overview
    assert "short answer" in overview
    assert "package artifact inventory" in overview
    assert "Biological Question" in overview_html
    assert "package reproducibility checklist" in overview_html
    assert package_manifest["biological_question"].startswith(
        "Do the host-associated rabies lineages"
    )
    assert package_manifest["workflow_files"]["final_report"]["path"] == (
        "workflow/rabies-cross-host-geography-report.html"
    )
    assert package_manifest["package_files"]["artifact_inventory"]["path"] == (
        "rabies-cross-host-geography-artifacts.tsv"
    )
    assert package_manifest["package_files"]["reproducibility_checklist"]["path"] == (
        "rabies-cross-host-geography-reproducibility-checklist.tsv"
    )
    assert package_manifest["package_files"]["artifact_inventory"][
        "artifact_count"
    ] == (len(inventory_rows))
    assert package_manifest["package_files"]["reproducibility_checklist"][
        "item_count"
    ] == len(checklist_rows)
    assert any(
        row["relative_path"] == "workflow/rabies-cross-host-geography-report.html"
        for row in inventory_rows
    )
    assert any(
        row["check_id"] == "downstream-analysis-exported" for row in checklist_rows
    )


def test_export_rabies_cross_host_geography_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_cross_host_geography_panel_dataset(tmp_path / "dataset")
    expected_files = {
        str(path.relative_to(result.expected_output_root))
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    }
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert result.metadata_path.is_file()
    assert result.centroids_path.is_file()
    assert result.accession_table_path.is_file()
    assert result.workflow_config_path.is_file()
    accession_rows = list(
        csv.DictReader(
            result.accession_table_path.open("r", encoding="utf-8", newline=""),
            delimiter="\t",
        )
    )
    assert accession_rows[0]["accession"] == "MG458305"
    assert accession_rows[0]["accession_url"].endswith("/MG458305")
    assert "rabies-cross-host-geography-report.html" in expected_files
    assert "biogeography/biogeography-report.html" in expected_files
    assert "comparative/comparative-report.html" in expected_files
    assert "bootstrap-review/bootstrap-review.summary.tsv" in expected_files
    assert (
        "bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv"
        in expected_files
    )
    assert "conclusion-stability/conclusion-stability-report.html" in expected_files
    assert "conclusion-stability/conclusion-stability-summary.tsv" in expected_files
    assert "workflow-config-audit.tsv" in expected_files
    assert "scientific-findings.tsv" in expected_files


def test_run_rabies_cross_host_geography_panel_demo_writes_flagship_package_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "demo"
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
    comparative_report_path = workflow_root / "comparative" / "comparative-report.html"
    biogeography_report_path = (
        workflow_root / "biogeography" / "biogeography-report.html"
    )
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

    result = run_rabies_cross_host_geography_panel_demo(output_root)

    assert result.dataset_export.accession_table_path.is_file()
    assert result.overview_html_path.is_file()
    assert result.package_manifest_path.is_file()
    payload = json.loads(result.package_manifest_path.read_text(encoding="utf-8"))
    assert payload["report_kind"] == "rabies_cross_host_geography_package"
    assert payload["dataset_files"]["source_accessions"]["path"] == (
        "dataset/source-accessions.tsv"
    )
    assert payload["package_files"]["artifact_inventory"]["path"] == (
        "rabies-cross-host-geography-artifacts.tsv"
    )
    assert payload["package_files"]["reproducibility_checklist"]["path"] == (
        "rabies-cross-host-geography-reproducibility-checklist.tsv"
    )
    assert payload["workflow_files"]["final_report"]["path"] == (
        "workflow/rabies-cross-host-geography-report.html"
    )


def test_cli_demo_rabies_cross_host_geography_panel_reports_flagship_package_metrics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys,
) -> None:
    output = tmp_path / "rabies-demo"
    expected_root = tmp_path / "expected"
    expected_root.mkdir(parents=True, exist_ok=True)
    for index in range(3):
        (expected_root / f"expected-{index}.txt").write_text("ok\n", encoding="utf-8")

    dataset = load_rabies_cross_host_geography_panel_dataset()
    dataset_export = RabiesCrossHostGeographyPanelExportResult(
        output_root=output / "dataset",
        readme_path=output / "dataset" / "README.md",
        workflow_config_path=output / "dataset" / "workflow-config.json",
        sequences_path=output / "dataset" / "sequences.fasta",
        metadata_path=output / "dataset" / "metadata.csv",
        centroids_path=output / "dataset" / "region-centroids.csv",
        accession_table_path=output / "dataset" / "source-accessions.tsv",
        expected_output_root=expected_root,
    )
    workflow_bundle = RabiesCrossHostGeographyPanelWorkflowBundle(
        output_root=output / "workflow",
        selected_model="TPM2u+F+G4",
        sequence_type="dna",
        inferred_sequence_type="dna",
        input_repair_applied=False,
        aligned_quality_score=90.48,
        trimmed_quality_score=90.48,
        minimum_support=84.0,
        maximum_support=100.0,
        median_support=97.0,
        weakly_supported_clade_count=0,
        clade_row_count=17,
        bootstrap_tree_count=1000,
        bootstrap_topology_count=4,
        bootstrap_unstable_branch_count=1,
        bootstrap_consensus_rooted_rf_distance=1,
        bootstrap_consensus_same_unrooted_topology=True,
        bootstrap_consensus_high_support_conflict_count=0,
        bootstrap_consensus_branch_score_distance=0.0,
        rooted_outgroup_taxa=("bat_chile_rv108",),
        root_host="bat",
        root_host_confidence=0.97,
        host_switch_count=2,
        certain_host_switch_count=1,
        uncertain_host_switch_count=1,
        root_region="north_asia",
        root_region_probability=0.98,
        changed_region_branch_count=4,
        migration_event_count=4,
        strongly_supported_migration_event_count=3,
        comparative_selected_model="brownian",
        comparative_response="region_longitude",
        comparative_formula="region_longitude ~ host_group",
        comparative_pgls_lambda=0.0,
        comparative_pgls_r_squared=0.61,
        comparative_branch_repair_count=0,
        conclusion_stable_count=6,
        conclusion_weak_count=3,
        conclusion_unstable_count=1,
        timeout_seconds=300.0,
        max_bootstrap_tree_count=1500,
        max_report_table_rows=25,
        memory_warning_threshold_bytes=67108864,
        workflow_runtime_seconds=12.5,
        bootstrap_review_runtime_seconds=0.25,
        bootstrap_review_peak_memory_bytes=2048,
        budget_warning_count=0,
        config_check_count=12,
        scientific_finding_count=6,
        workflow_summary_path=output / "workflow" / "workflow-summary.tsv",
        resource_observations_path=output / "workflow" / "resource-observations.tsv",
        config_audit_path=output / "workflow" / "workflow-config-audit.tsv",
        resolved_config_path=output / "workflow" / "workflow-config.resolved.json",
        input_validation_path=output / "workflow" / "input-validation.tsv",
        alignment_quality_path=output / "workflow" / "alignment-quality.tsv",
        alignment_sequence_ranking_path=output
        / "workflow"
        / "alignment-sequence-ranking.tsv",
        alignment_path=output / "workflow" / "rabies-cross-host-geography-panel.aln",
        trimmed_alignment_path=output
        / "workflow"
        / "rabies-cross-host-geography-panel.trimmed.aln",
        tree_path=output / "workflow" / "rabies-cross-host-geography-panel.rooted.tree",
        rooting_report_path=output
        / "workflow"
        / "rabies-cross-host-geography-panel.rooting.tsv",
        model_table_path=output
        / "workflow"
        / "rabies-cross-host-geography-panel.model.tsv",
        support_table_path=output
        / "workflow"
        / "rabies-cross-host-geography-panel.support.tsv",
        log_path=output / "workflow" / "rabies-cross-host-geography-panel.log",
        manifest_path=output
        / "workflow"
        / "rabies-cross-host-geography-panel.manifest.json",
        engine_artifact_root=output / "workflow" / "engine-artifacts",
        clade_table_path=output / "workflow" / "clade-table.tsv",
        bootstrap_output_root=output / "workflow" / "bootstrap-review",
        bootstrap_summary_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-review.summary.tsv",
        bootstrap_consensus_tree_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-consensus.tree",
        bootstrap_clade_frequencies_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-clade-frequencies.tsv",
        bootstrap_unstable_branches_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-unstable-branches.tsv",
        bootstrap_unstable_clades_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-unstable-clades.tsv",
        bootstrap_distance_matrix_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-review.distance-matrix.tsv",
        bootstrap_topology_clusters_path=output
        / "workflow"
        / "bootstrap-review"
        / "bootstrap-topology-clusters.tsv",
        bootstrap_tree_comparison_summary_path=output
        / "workflow"
        / "bootstrap-review"
        / "rooted-tree-vs-bootstrap-consensus.summary.tsv",
        bootstrap_tree_comparison_table_path=output
        / "workflow"
        / "bootstrap-review"
        / "rooted-tree-vs-bootstrap-consensus.table.tsv",
        bootstrap_tree_comparison_report_path=output
        / "workflow"
        / "bootstrap-review"
        / "rooted-tree-vs-bootstrap-consensus.report.html",
        host_switch_summary_path=output / "workflow" / "host-switch-summary.tsv",
        host_state_nodes_path=output / "workflow" / "host-state-nodes.tsv",
        host_switch_branches_path=output / "workflow" / "host-switch-branches.tsv",
        host_switch_counts_path=output / "workflow" / "host-switch-counts.tsv",
        host_switch_fits_path=output / "workflow" / "host-switch-fits.tsv",
        host_switch_unsupported_path=output
        / "workflow"
        / "host-switch-unsupported.tsv",
        host_switch_exclusions_path=output / "workflow" / "host-switch-exclusions.tsv",
        biogeography_output_root=output / "workflow" / "biogeography",
        biogeography_report_path=output
        / "workflow"
        / "biogeography"
        / "biogeography-report.html",
        biogeography_tree_figure_path=output
        / "workflow"
        / "biogeography"
        / "ancestral-region-tree.svg",
        biogeography_map_path=output
        / "workflow"
        / "biogeography"
        / "geographic-region-map.html",
        comparative_traits_path=output
        / "workflow"
        / "comparative"
        / "comparative-traits.tsv",
        comparative_tree_path=output
        / "workflow"
        / "comparative"
        / "comparative-tree.nwk",
        comparative_repairs_path=output
        / "workflow"
        / "comparative"
        / "comparative-branch-repairs.tsv",
        comparative_output_root=output / "workflow" / "comparative",
        comparative_report_path=output
        / "workflow"
        / "comparative"
        / "comparative-report.html",
        comparative_summary_path=output
        / "workflow"
        / "comparative"
        / "comparative-summary.tsv",
        comparative_coefficients_path=output
        / "workflow"
        / "comparative"
        / "comparative-coefficients.tsv",
        comparative_residuals_path=output
        / "workflow"
        / "comparative"
        / "comparative-residuals.tsv",
        comparative_signal_path=output
        / "workflow"
        / "comparative"
        / "comparative-signal.tsv",
        comparative_model_comparison_path=output
        / "workflow"
        / "comparative"
        / "comparative-model-comparison.tsv",
        comparative_interpretation_path=output
        / "workflow"
        / "comparative"
        / "comparative-interpretation.tsv",
        comparative_audit_path=output
        / "workflow"
        / "comparative"
        / "comparative-audit.tsv",
        comparative_contrasts_path=output
        / "workflow"
        / "comparative"
        / "comparative-contrasts.tsv",
        comparative_model_matrix_path=output
        / "workflow"
        / "comparative"
        / "comparative-model-matrix.tsv",
        comparative_categorical_contrasts_path=output
        / "workflow"
        / "comparative"
        / "comparative-categorical-contrasts.tsv",
        comparative_lambda_profile_path=output
        / "workflow"
        / "comparative"
        / "comparative-lambda-profile.tsv",
        comparative_manifest_path=output
        / "workflow"
        / "comparative"
        / "comparative.manifest.json",
        conclusion_stability_output_root=output / "workflow" / "conclusion-stability",
        conclusion_stability_summary_path=output
        / "workflow"
        / "conclusion-stability"
        / "conclusion-stability-summary.tsv",
        key_clade_stability_path=output
        / "workflow"
        / "conclusion-stability"
        / "key-clade-stability.tsv",
        support_value_stability_path=output
        / "workflow"
        / "conclusion-stability"
        / "support-value-stability.tsv",
        ancestral_state_stability_path=output
        / "workflow"
        / "conclusion-stability"
        / "ancestral-state-stability.tsv",
        comparative_coefficient_stability_path=output
        / "workflow"
        / "conclusion-stability"
        / "comparative-coefficient-stability.tsv",
        conclusion_stability_report_path=output
        / "workflow"
        / "conclusion-stability"
        / "conclusion-stability-report.html",
        scientific_findings_path=output / "workflow" / "scientific-findings.tsv",
        final_report_path=output
        / "workflow"
        / "rabies-cross-host-geography-report.html",
        final_manifest_path=output
        / "workflow"
        / "rabies-cross-host-geography.manifest.json",
    )
    fake_result = RabiesCrossHostGeographyPanelDemoResult(
        output_root=output,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=output / "overview.md",
        overview_html_path=output / "rabies-cross-host-geography-overview.html",
        artifact_inventory_path=output / "rabies-cross-host-geography-artifacts.tsv",
        reproducibility_checklist_path=output
        / "rabies-cross-host-geography-reproducibility-checklist.tsv",
        package_manifest_path=output
        / "rabies-cross-host-geography-package.manifest.json",
    )

    for path in [
        dataset_export.readme_path,
        dataset_export.workflow_config_path,
        dataset_export.sequences_path,
        dataset_export.metadata_path,
        dataset_export.centroids_path,
        dataset_export.accession_table_path,
        workflow_bundle.workflow_summary_path,
        workflow_bundle.resource_observations_path,
        workflow_bundle.config_audit_path,
        workflow_bundle.resolved_config_path,
        workflow_bundle.input_validation_path,
        workflow_bundle.alignment_quality_path,
        workflow_bundle.alignment_sequence_ranking_path,
        workflow_bundle.alignment_path,
        workflow_bundle.trimmed_alignment_path,
        workflow_bundle.tree_path,
        workflow_bundle.rooting_report_path,
        workflow_bundle.model_table_path,
        workflow_bundle.support_table_path,
        workflow_bundle.clade_table_path,
        workflow_bundle.bootstrap_summary_path,
        workflow_bundle.bootstrap_tree_comparison_summary_path,
        workflow_bundle.host_switch_summary_path,
        workflow_bundle.host_switch_counts_path,
        workflow_bundle.biogeography_report_path,
        workflow_bundle.biogeography_tree_figure_path,
        workflow_bundle.biogeography_map_path,
        workflow_bundle.comparative_report_path,
        workflow_bundle.comparative_summary_path,
        workflow_bundle.conclusion_stability_summary_path,
        workflow_bundle.key_clade_stability_path,
        workflow_bundle.support_value_stability_path,
        workflow_bundle.ancestral_state_stability_path,
        workflow_bundle.comparative_coefficient_stability_path,
        workflow_bundle.conclusion_stability_report_path,
        workflow_bundle.scientific_findings_path,
        workflow_bundle.final_report_path,
        workflow_bundle.final_manifest_path,
        fake_result.overview_path,
        fake_result.overview_html_path,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("ok\n", encoding="utf-8")
    fake_result.artifact_inventory_path.write_text(
        "section\tkind\trelative_path\tsha256\tsize_bytes\n"
        "package\tmarkdown\toverview.md\tabc\t3\n",
        encoding="utf-8",
    )
    fake_result.reproducibility_checklist_path.write_text(
        "section\tcheck_id\tstatus\tsummary\tevidence\tartifact_path\n"
        "package\tpackage-navigation-exported\tpass\tok\tok\toverview.md\n",
        encoding="utf-8",
    )
    fake_result.package_manifest_path.write_text("ok\n", encoding="utf-8")

    monkeypatch.setattr(
        "bijux_phylogenetics.command_line.demo.run_rabies_cross_host_geography_panel_demo",
        lambda *args, **kwargs: fake_result,
    )

    exit_code = main(
        [
            "demo",
            "rabies-cross-host-geography-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["artifact_count"] == 43
    assert payload["metrics"]["reference_output_count"] == 3
    assert payload["metrics"]["package_artifact_count"] == 1
    assert payload["metrics"]["package_checklist_item_count"] == 1
    assert payload["metrics"]["biological_question"].startswith(
        "Do the host-associated rabies lineages"
    )
    assert (
        "nominally supported positive longitude association"
        in payload["metrics"]["short_answer"]
    )
    assert payload["data"]["dataset_export"]["accession_table_path"] == str(
        output / "dataset" / "source-accessions.tsv"
    )
    assert payload["data"]["overview_html_path"] == str(
        output / "rabies-cross-host-geography-overview.html"
    )
    assert payload["data"]["package_manifest_path"] == str(
        output / "rabies-cross-host-geography-package.manifest.json"
    )
    assert payload["metrics"]["conclusion_stable_count"] == 6
    assert payload["metrics"]["conclusion_weak_count"] == 3
    assert payload["metrics"]["conclusion_unstable_count"] == 1


def test_public_runtime_exports_include_rabies_cross_host_geography_panel_surface() -> (
    None
):
    assert (
        rabies_host_geography_api.load_rabies_cross_host_geography_panel_dataset
        is load_rabies_cross_host_geography_panel_dataset
    )
    assert (
        rabies_host_geography_api.export_rabies_cross_host_geography_panel_dataset
        is export_rabies_cross_host_geography_panel_dataset
    )
    assert (
        rabies_host_geography_api.run_rabies_cross_host_geography_panel_workflow
        is run_rabies_cross_host_geography_panel_workflow
    )
    assert (
        rabies_host_geography_api.write_rabies_cross_host_geography_panel_workflow_bundle
        is write_rabies_cross_host_geography_panel_workflow_bundle
    )
    assert (
        rabies_host_geography_api.run_rabies_cross_host_geography_panel_demo
        is run_rabies_cross_host_geography_panel_demo
    )


@pytest.mark.slow
def test_cli_demo_rabies_cross_host_geography_panel_json_output_reports_integrated_workflow(
    tmp_path: Path, capsys
) -> None:
    executables = require_alignment_engine_executables()
    output = tmp_path / "rabies-integrated-demo"
    dataset = load_rabies_cross_host_geography_panel_dataset()
    exit_code = main(
        [
            "demo",
            "rabies-cross-host-geography-panel",
            "--out",
            str(output),
            "--config",
            str(dataset.workflow_config_path),
            "--mafft-executable",
            executables["mafft"],
            "--trimal-executable",
            executables["trimal"],
            "--iqtree-executable",
            executables["iqtree2"],
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 43
    assert payload["metrics"]["package_artifact_count"] > 0
    assert payload["metrics"]["package_checklist_item_count"] == 6
    assert payload["metrics"]["sequence_count"] == 9
    assert payload["metrics"]["config_path"] == str(
        output / "dataset" / "workflow-config.json"
    )
    assert payload["metrics"]["biological_question"].startswith(
        "Do the host-associated rabies lineages"
    )
    assert (
        "nominally supported positive longitude association"
        in payload["metrics"]["short_answer"]
    )
    assert payload["metrics"]["host_trait"] == "host_group"
    assert payload["metrics"]["geography_trait"] == "region_group"
    assert payload["metrics"]["selected_model"] == "TPM2u+F+G4"
    assert payload["metrics"]["aligned_quality_score"] == 90.48
    assert payload["metrics"]["trimmed_quality_score"] == 90.48
    assert payload["metrics"]["minimum_support"] == 84.0
    assert payload["metrics"]["maximum_support"] == 100.0
    assert payload["metrics"]["root_host"] == "bat"
    assert payload["metrics"]["root_region"] == "north_asia"
    assert payload["metrics"]["host_switch_count"] == 2
    assert payload["metrics"]["migration_event_count"] == 4
    assert payload["metrics"]["clade_row_count"] == 17
    assert payload["metrics"]["bootstrap_tree_count"] == 1000
    assert payload["metrics"]["bootstrap_consensus_rooted_rf_distance"] == 1
    assert payload["metrics"]["comparative_formula"] == "region_longitude ~ host_group"
    assert payload["metrics"]["comparative_selected_model"] == "brownian"
    assert payload["metrics"]["config_check_count"] == 12
    assert payload["metrics"]["scientific_finding_count"] == 6
    assert payload["metrics"]["reference_output_count"] == 66
    assert payload["metrics"]["conclusion_stable_count"] == 30
    assert payload["metrics"]["conclusion_weak_count"] == 0
    assert payload["metrics"]["conclusion_unstable_count"] == 2
    assert payload["data"]["dataset"]["dataset_id"] == (
        "rabies_cross_host_geography_panel"
    )
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
    assert payload["data"]["dataset_export"]["accession_table_path"] == str(
        output / "dataset" / "source-accessions.tsv"
    )
    assert payload["data"]["overview_html_path"] == str(
        output / "rabies-cross-host-geography-overview.html"
    )
    assert payload["data"]["package_manifest_path"] == str(
        output / "rabies-cross-host-geography-package.manifest.json"
    )
