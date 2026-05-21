from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .contracts import GeographicSamplingBiasNodeRow, GeographicSamplingCountRow

DOMINANT_REGION_THRESHOLD = 0.8


def included_region_counts(audit) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in audit.rows:
        if not row.included or row.normalized_state is None:
            continue
        counts[row.normalized_state] = counts.get(row.normalized_state, 0) + 1
    if not counts:
        raise AncestralReconstructionError(
            "geographic sampling-bias review requires at least one usable observed region"
        )
    return dict(sorted(counts.items()))


def resolve_region_weights(
    counts: dict[str, int],
    *,
    weights_path: Path | None,
    region_column: str,
    weight_column: str,
) -> tuple[dict[str, float], str]:
    if weights_path is None:
        total = sum(counts.values())
        region_count = len(counts)
        return (
            {
                region: stable_value(total / max(region_count * count, 1))
                for region, count in counts.items()
            },
            "inverse-frequency",
        )
    table = load_taxon_table(weights_path, taxon_column=region_column)
    if weight_column not in table.columns:
        raise AncestralReconstructionError(
            f"region weight table does not contain column '{weight_column}'"
        )
    weights: dict[str, float] = {}
    for row in table.rows:
        region = row[table.taxon_column]
        raw_weight = row[weight_column]
        try:
            weight = float(raw_weight)
        except ValueError as error:
            raise AncestralReconstructionError(
                f"region weight for '{region}' is not numeric: {raw_weight!r}"
            ) from error
        if not math.isfinite(weight) or weight <= 0.0:
            raise AncestralReconstructionError(
                f"region weight for '{region}' must be finite and positive"
            )
        weights[region] = stable_value(weight)
    missing = sorted(region for region in counts if region not in weights)
    if missing:
        raise AncestralReconstructionError(
            "region weight table is missing observed regions: " + ", ".join(missing)
        )
    return ({region: weights[region] for region in sorted(counts)}, "explicit")


def build_count_rows(
    counts: dict[str, int],
    weights: dict[str, float],
) -> list[GeographicSamplingCountRow]:
    total = sum(counts.values())
    weighted_counts = {
        region: stable_value(count * weights[region])
        for region, count in counts.items()
    }
    weighted_total = sum(weighted_counts.values())
    dominant_fraction = max(counts.values()) / max(total, 1)
    weighted_dominant_fraction = max(weighted_counts.values()) / max(
        weighted_total, 1.0
    )
    return [
        GeographicSamplingCountRow(
            region=region,
            sample_count=count,
            sample_fraction=stable_value(count / max(total, 1)),
            applied_weight=weights[region],
            weighted_sample_count=weighted_counts[region],
            weighted_sample_fraction=stable_value(
                weighted_counts[region] / max(weighted_total, 1.0)
            ),
            dominant_unweighted=stable_value(count / max(total, 1))
            == stable_value(dominant_fraction),
            dominant_weighted=stable_value(
                weighted_counts[region] / max(weighted_total, 1.0)
            )
            == stable_value(weighted_dominant_fraction),
        )
        for region, count in sorted(counts.items())
    ]


def build_warnings(
    count_rows: list[GeographicSamplingCountRow],
    weighting_mode: str,
    node_rows: list[GeographicSamplingBiasNodeRow],
) -> list[str]:
    warnings: list[str] = []
    dominant_row = max(count_rows, key=lambda row: (row.sample_fraction, row.region))
    if dominant_row.sample_fraction >= DOMINANT_REGION_THRESHOLD:
        warnings.append(
            "observed regions are dominated by one sampled region and the baseline reconstruction may reflect sample imbalance"
        )
    weighted_dominant_row = max(
        count_rows,
        key=lambda row: (row.weighted_sample_fraction, row.region),
    )
    if weighting_mode == "inverse-frequency":
        warnings.append(
            "inverse-frequency region weights rebalance observed regions to equal weighted mass before comparing ancestral conclusions"
        )
    else:
        warnings.append(
            "explicit region weights reweight ancestral region probabilities and branchwise transitions for sampling-bias review"
        )
    if weighted_dominant_row.weighted_sample_fraction >= DOMINANT_REGION_THRESHOLD:
        warnings.append(
            "weighted region counts remain dominated by one region after correction"
        )
    if any(row.changed and row.is_root for row in node_rows):
        warnings.append(
            "the weighted correction changes the most likely root region relative to the baseline reconstruction"
        )
    return warnings
