from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.biogeography.migration import (
    summarize_geographic_migration_events,
    write_geographic_migration_event_table,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    build_comparative_report_package,
)
from bijux_phylogenetics.datasets.rabies_cross_host_geography.workflow.comparative_inputs import (
    _build_comparative_trait_rows,
)
from bijux_phylogenetics.datasets.rabies_host_geography import (
    load_rabies_cross_host_geography_panel_dataset,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.artifact_schema import (
    validate_artifact_schema,
    validate_comparative_report_manifest_schema,
    validate_comparative_summary_table_schema,
    validate_comparative_traits_table_schema,
    validate_geographic_event_table_schema,
)

DATASET_CONFIG = (
    Path(__file__).parent.parent
    / "src"
    / "bijux_phylogenetics"
    / "resources"
    / "datasets"
    / "pathogens"
    / "rabies_cross_host_geography_panel"
    / "workflow-config.json"
)


def _dataset():
    return load_rabies_cross_host_geography_panel_dataset(DATASET_CONFIG)


def _rooted_tree_path(dataset) -> Path:
    return dataset.reference_output_root / f"{dataset.workflow_prefix}.rooted.tree"


def _comparative_tree_path(dataset) -> Path:
    return dataset.reference_output_root / "comparative-tree.nwk"


def _write_live_comparative_traits(path: Path) -> Path:
    dataset = _dataset()
    rows = _build_comparative_trait_rows(
        metadata_path=dataset.metadata_path,
        centroids_path=dataset.centroids_path,
        host_trait=dataset.host_trait,
        geography_trait=dataset.geography_trait,
    )
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def test_comparative_traits_schema_matches_live_output(tmp_path: Path) -> None:
    traits_path = _write_live_comparative_traits(tmp_path / "comparative-traits.tsv")

    assert validate_comparative_traits_table_schema(traits_path).valid
    assert validate_artifact_schema(traits_path, "comparative_traits_tsv").valid


def test_comparative_summary_and_manifest_schemas_match_live_output(
    tmp_path: Path,
) -> None:
    dataset = _dataset()
    traits_path = _write_live_comparative_traits(tmp_path / "comparative-traits.tsv")
    package = build_comparative_report_package(
        _comparative_tree_path(dataset),
        traits_path,
        out_dir=tmp_path / "comparative",
        formula=dataset.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )

    assert validate_comparative_summary_table_schema(package.summary_table_path).valid
    assert validate_comparative_report_manifest_schema(package.manifest_path).valid
    assert validate_artifact_schema(
        package.summary_table_path,
        "comparative_summary_tsv",
    ).valid
    assert validate_artifact_schema(
        package.manifest_path,
        "comparative_report_manifest_json",
    ).valid
    manifest = json.loads(package.manifest_path.read_text(encoding="utf-8"))
    assert package.methods_summary_path.exists()
    assert manifest["outputs"]["methods_summary_path"] == str(
        package.methods_summary_path
    )


def test_geographic_event_schema_matches_live_output(tmp_path: Path) -> None:
    dataset = _dataset()
    report = summarize_geographic_migration_events(
        _rooted_tree_path(dataset),
        dataset.metadata_path,
        trait=dataset.geography_trait,
        model=dataset.geography_model,
    )
    path = write_geographic_migration_event_table(tmp_path / "events.tsv", report)

    assert validate_geographic_event_table_schema(path).valid
    assert validate_artifact_schema(path, "geographic_event_tsv").valid
