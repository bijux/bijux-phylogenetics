from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path

from ..models import GeneticDistanceMatrix


def write_genetic_distance_matrix(path: Path, report: GeneticDistanceMatrix) -> Path:
    """Write a pairwise genetic distance matrix as a deterministic TSV."""
    rows = {
        (pair.left_identifier, pair.right_identifier): pair for pair in report.pairs
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tdistance\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            normalized_distance = (
                None
                if pair.distance is None
                else 0.0
                if math.isclose(pair.distance, 0.0, abs_tol=1e-15)
                else pair.distance
            )
            distance = (
                ""
                if normalized_distance is None
                else format(normalized_distance, ".15g")
            )
            lines.append(f"{left}\t{right}\t{distance}\t{pair.comparable_sites}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_genetic_distance_component_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    """Write one deterministic pairwise distance component table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "left_identifier",
                "right_identifier",
                "distance",
                "comparable_sites",
                "mismatch_sites",
                "transition_sites",
                "ag_transition_sites",
                "ct_transition_sites",
                "transversion_sites",
                "ambiguity_sites",
                "skipped_sites",
                "saturated",
                "saturation_reason",
            ]
        )
    ]
    for pair in report.pairs:
        normalized_distance = (
            None
            if pair.distance is None
            else 0.0
            if math.isclose(pair.distance, 0.0, abs_tol=1e-15)
            else pair.distance
        )
        lines.append(
            "\t".join(
                [
                    pair.left_identifier,
                    pair.right_identifier,
                    ""
                    if normalized_distance is None
                    else format(normalized_distance, ".15g"),
                    str(pair.comparable_sites),
                    format(pair.mismatch_sites, ".15g"),
                    format(pair.transition_sites, ".15g"),
                    format(pair.ag_transition_sites, ".15g"),
                    format(pair.ct_transition_sites, ".15g"),
                    format(pair.transversion_sites, ".15g"),
                    str(pair.ambiguity_sites),
                    str(pair.skipped_sites),
                    "true" if pair.saturated else "false",
                    "" if pair.saturation_reason is None else pair.saturation_reason,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_genetic_distance_parameter_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    """Write one deterministic alignment-wide parameter table for DNA distances."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["parameter\tvalue"]
    if report.model_parameters is None:
        lines.append(f"model\t{report.model}")
    else:
        for parameter, value in asdict(report.model_parameters).items():
            rendered = "" if value is None else format(value, ".15g")
            lines.append(f"{parameter}\t{rendered}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
