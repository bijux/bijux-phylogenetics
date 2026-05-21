from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.datasets.rabies_host_geography import (
    load_rabies_cross_host_geography_panel_dataset,
)
from bijux_phylogenetics.ecology.host_switching import (
    summarize_host_switching,
    write_host_switch_branch_table,
)
from bijux_phylogenetics.io.artifact_schema import (
    validate_artifact_schema,
    validate_clade_table_schema,
    validate_host_switch_branch_table_schema,
)
from bijux_phylogenetics.trees import extract_tree_clades, write_clade_table

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


def _rooted_tree_path(dataset) -> Path:
    return dataset.reference_output_root / f"{dataset.workflow_prefix}.rooted.tree"


def test_clade_table_schema_matches_live_output(tmp_path: Path) -> None:
    dataset = load_rabies_cross_host_geography_panel_dataset(DATASET_CONFIG)
    report = extract_tree_clades(
        _rooted_tree_path(dataset),
        metadata_path=dataset.metadata_path,
        metadata_columns=list(dataset.clade_metadata_columns),
    )
    path = write_clade_table(tmp_path / "clades.tsv", report)

    schema = validate_clade_table_schema(path)

    assert schema.valid
    assert (
        f"metadata_column_count={len(dataset.clade_metadata_columns)}" in schema.notes
    )
    assert validate_artifact_schema(path, "clade_table_tsv").valid


@pytest.mark.slow
def test_host_switch_branch_schema_matches_live_output(tmp_path: Path) -> None:
    dataset = load_rabies_cross_host_geography_panel_dataset(DATASET_CONFIG)
    report = summarize_host_switching(
        _rooted_tree_path(dataset),
        dataset.metadata_path,
        trait=dataset.host_trait,
        model=dataset.host_model,
    )
    path = write_host_switch_branch_table(tmp_path / "host-switch-branches.tsv", report)

    assert validate_host_switch_branch_table_schema(path).valid
    assert validate_artifact_schema(path, "host_switch_branch_tsv").valid
