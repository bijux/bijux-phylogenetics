from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.rabies_geography import (
    export_rabies_geographic_transition_panel_dataset,
    load_rabies_geographic_transition_panel_dataset,
    run_rabies_geographic_transition_panel_demo,
    run_rabies_geographic_transition_panel_workflow,
    write_rabies_geographic_transition_panel_workflow_bundle,
)


def test_load_rabies_geographic_transition_panel_dataset_exposes_packaged_surface(
) -> None:
    dataset = load_rabies_geographic_transition_panel_dataset()
    assert dataset.dataset_id == "rabies_geographic_transition_panel"
    assert dataset.label == "Rabies geographic transition panel"
    assert dataset.taxon_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.workflow_trait == "region_group"
    assert dataset.workflow_model == "ard"
    assert dataset.observed_region_group_count == 5
    assert dataset.sequences_path.is_file()
    assert dataset.tree_path.is_file()
    assert dataset.regions_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "MG458305" in dataset.source_accessions


def test_write_rabies_geographic_transition_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_rabies_geographic_transition_panel_workflow()
    bundle = write_rabies_geographic_transition_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.workflow_summary_path.name: bundle.workflow_summary_path,
        bundle.geographic_state_summary_path.name: bundle.geographic_state_summary_path,
        bundle.geographic_region_probability_path.name: (
            bundle.geographic_region_probability_path
        ),
        bundle.geographic_transition_rate_path.name: (
            bundle.geographic_transition_rate_path
        ),
        bundle.geographic_transition_event_path.name: (
            bundle.geographic_transition_event_path
        ),
        bundle.geographic_state_exclusion_path.name: (
            bundle.geographic_state_exclusion_path
        ),
        bundle.geographic_migration_summary_path.name: (
            bundle.geographic_migration_summary_path
        ),
        bundle.geographic_migration_event_path.name: (
            bundle.geographic_migration_event_path
        ),
        bundle.geographic_migration_exclusion_path.name: (
            bundle.geographic_migration_exclusion_path
        ),
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    for name, generated_path in generated.items():
        assert generated_path.read_text(encoding="utf-8") == (
            expected_root / name
        ).read_text(encoding="utf-8")


def test_run_rabies_geographic_transition_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_rabies_geographic_transition_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 9
    assert result.dataset_export.sequences_path.is_file()
    assert result.dataset_export.tree_path.is_file()
    assert result.dataset_export.regions_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.geographic_state_summary_path.is_file()
    assert result.workflow_bundle.geographic_migration_summary_path.is_file()
    assert result.overview_path.is_file()
    assert "geographic-state summary" in result.overview_path.read_text(
        encoding="utf-8"
    )


def test_export_rabies_geographic_transition_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_geographic_transition_panel_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.sequences_path.is_file()
    assert result.tree_path.is_file()
    assert result.regions_path.is_file()
    assert len(expected_files) == 9
    assert "geographic-migration-events.tsv" in expected_files
