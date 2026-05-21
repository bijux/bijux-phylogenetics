from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .contracts import CoordinateMovementExclusionRow


@dataclass(slots=True)
class PreparedCoordinateDataset:
    taxon_column: str
    rows_by_taxon: dict[str, tuple[float, float]]
    included_taxa: list[str]
    exclusion_rows: list[CoordinateMovementExclusionRow]


def prepare_coordinate_dataset(
    tree_path: Path,
    table_path: Path,
    *,
    latitude_column: str,
    longitude_column: str,
    taxon_column: str | None,
) -> PreparedCoordinateDataset:
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    if latitude_column not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{latitude_column}'"
        )
    if longitude_column not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{longitude_column}'"
        )
    from bijux_phylogenetics.io.trees import load_tree

    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    rows_by_taxon: dict[str, tuple[float, float]] = {}
    exclusion_rows: list[CoordinateMovementExclusionRow] = []
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_latitude = row[latitude_column].strip()
        raw_longitude = row[longitude_column].strip()
        if taxon not in tree_taxa:
            exclusion_rows.append(
                CoordinateMovementExclusionRow(
                    taxon=taxon,
                    raw_latitude=raw_latitude,
                    raw_longitude=raw_longitude,
                    reason="taxon-not-in-tree",
                    note="taxon does not overlap the tree and is excluded before coordinate reconstruction",
                )
            )
            continue
        if not raw_latitude or not raw_longitude:
            exclusion_rows.append(
                CoordinateMovementExclusionRow(
                    taxon=taxon,
                    raw_latitude=raw_latitude,
                    raw_longitude=raw_longitude,
                    reason="missing-coordinate",
                    note="latitude or longitude is blank and cannot be used",
                )
            )
            continue
        try:
            latitude = float(raw_latitude)
            longitude = float(raw_longitude)
        except ValueError:
            exclusion_rows.append(
                CoordinateMovementExclusionRow(
                    taxon=taxon,
                    raw_latitude=raw_latitude,
                    raw_longitude=raw_longitude,
                    reason="non-numeric-coordinate",
                    note="latitude and longitude must be numeric",
                )
            )
            continue
        if latitude < -90.0 or latitude > 90.0:
            exclusion_rows.append(
                CoordinateMovementExclusionRow(
                    taxon=taxon,
                    raw_latitude=raw_latitude,
                    raw_longitude=raw_longitude,
                    reason="latitude-out-of-range",
                    note="latitude must be within [-90, 90]",
                )
            )
            continue
        if longitude < -180.0 or longitude > 180.0:
            exclusion_rows.append(
                CoordinateMovementExclusionRow(
                    taxon=taxon,
                    raw_latitude=raw_latitude,
                    raw_longitude=raw_longitude,
                    reason="longitude-out-of-range",
                    note="longitude must be within [-180, 180]",
                )
            )
            continue
        rows_by_taxon[taxon] = (latitude, longitude)
    included_taxa = sorted(rows_by_taxon)
    return PreparedCoordinateDataset(
        taxon_column=table.taxon_column,
        rows_by_taxon=rows_by_taxon,
        included_taxa=included_taxa,
        exclusion_rows=exclusion_rows,
    )


def write_filtered_coordinate_table(
    path: Path,
    prepared: PreparedCoordinateDataset,
) -> None:
    write_taxon_rows(
        path,
        columns=[prepared.taxon_column, "latitude", "longitude"],
        rows=[
            {
                prepared.taxon_column: taxon,
                "latitude": str(prepared.rows_by_taxon[taxon][0]),
                "longitude": str(prepared.rows_by_taxon[taxon][1]),
            }
            for taxon in prepared.included_taxa
        ],
    )
