from __future__ import annotations

import csv
from pathlib import Path


def _build_comparative_trait_rows(
    *,
    metadata_path: Path,
    centroids_path: Path,
    host_trait: str,
    geography_trait: str,
) -> list[dict[str, str]]:
    centroids_by_region: dict[str, dict[str, str]] = {}
    with centroids_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            centroids_by_region[row["region"].strip()] = row
    rows: list[dict[str, str]] = []
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            region = row[geography_trait].strip()
            centroid = centroids_by_region[region]
            rows.append(
                {
                    "taxon": row["taxon"].strip(),
                    "host_group": row[host_trait].strip(),
                    "region_group": region,
                    "region_latitude": centroid["latitude"].strip(),
                    "region_longitude": centroid["longitude"].strip(),
                }
            )
    return rows
