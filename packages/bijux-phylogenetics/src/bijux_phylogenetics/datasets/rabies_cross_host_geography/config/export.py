from __future__ import annotations

import csv
from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    _WORKFLOW_CONFIG_NAME,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelExportResult,
)
from .panel_dataset import load_rabies_cross_host_geography_panel_dataset


def export_rabies_cross_host_geography_panel_dataset(
    destination: Path,
    *,
    config_path: Path | None = None,
) -> RabiesCrossHostGeographyPanelExportResult:
    """Copy the packaged integrated rabies dataset and stable expected outputs."""
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    workflow_config_path = shutil.copy2(
        dataset.workflow_config_path, destination / _WORKFLOW_CONFIG_NAME
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path, destination / "sequences.fasta"
    )
    metadata_path = shutil.copy2(dataset.metadata_path, destination / "metadata.csv")
    centroids_path = shutil.copy2(
        dataset.centroids_path, destination / "region-centroids.csv"
    )
    accession_table_path = _write_source_accession_table(
        destination / "source-accessions.tsv",
        dataset=dataset,
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesCrossHostGeographyPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        workflow_config_path=Path(workflow_config_path),
        sequences_path=Path(sequences_path),
        metadata_path=Path(metadata_path),
        centroids_path=Path(centroids_path),
        accession_table_path=accession_table_path,
        expected_output_root=expected_output_root,
    )


def _write_source_accession_table(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
) -> Path:
    metadata_rows = _read_metadata_rows(dataset.metadata_path)
    accession_index = {
        str(row["accession"]): row for row in metadata_rows if row.get("accession")
    }
    ordered_rows = []
    for accession in dataset.source_accessions:
        row = accession_index[accession]
        ordered_rows.append(
            {
                "accession": accession,
                "accession_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "taxon": str(row["taxon"]),
                "isolate": str(row["isolate"]),
                "host_species": str(row["host_species"]),
                "host_group": str(row["host_group"]),
                "country": str(row["country"]),
                "region_group": str(row["region_group"]),
                "collection_date": str(row.get("collection_date", "")),
            }
        )
    return write_taxon_rows(
        path,
        columns=[
            "accession",
            "accession_url",
            "taxon",
            "isolate",
            "host_species",
            "host_group",
            "country",
            "region_group",
            "collection_date",
        ],
        rows=ordered_rows,
    )


def _read_metadata_rows(metadata_path: Path) -> list[dict[str, str]]:
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))
