from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.runtime.errors import MetadataJoinError


def _split_csv_values(raw: str | None) -> list[str]:
    if raw is None:
        return []
    return [item.strip() for item in raw.split(",") if item.strip()]


def _parse_float_csv_row(raw: str) -> list[float]:
    values = _split_csv_values(raw)
    if not values:
        raise ValueError("matrix rows must contain one or more comma-separated numbers")
    try:
        return [float(value) for value in values]
    except ValueError as error:
        raise ValueError(
            f"matrix rows must contain only numeric values, got '{raw}'"
        ) from error


def _parse_assignment_map(raw: str | None) -> dict[str, str]:
    assignments: dict[str, str] = {}
    for item in _split_csv_values(raw):
        if "=" not in item:
            raise ValueError(f"mapping item must be KEY=VALUE, got '{item}'")
        key, value = item.split("=", 1)
        if not key.strip() or not value.strip():
            raise ValueError(
                f"mapping item must include both KEY and VALUE, got '{item}'"
            )
        assignments[key.strip()] = value.strip()
    return assignments


def _parse_transition_pairs(raw: str | None) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for item in _split_csv_values(raw):
        if "->" not in item:
            raise ValueError(f"transition item must be SOURCE->TARGET, got '{item}'")
        source_state, target_state = item.split("->", 1)
        if not source_state.strip() or not target_state.strip():
            raise ValueError(
                f"transition item must include both SOURCE and TARGET, got '{item}'"
            )
        pairs.append((source_state.strip(), target_state.strip()))
    return pairs


def _parse_rate_rows(raw_rows: list[str]) -> list[DiscreteHistoryRateRow]:
    from bijux_phylogenetics.simulation import DiscreteHistoryRateRow

    rows: list[DiscreteHistoryRateRow] = []
    for raw in raw_rows:
        if "=" not in raw or "->" not in raw:
            raise ValueError(
                f"rate item must be in SOURCE->TARGET=RATE form, got '{raw}'"
            )
        transition, raw_rate = raw.split("=", 1)
        source_state, target_state = transition.split("->", 1)
        if not source_state.strip() or not target_state.strip():
            raise ValueError(
                f"rate item must include both SOURCE and TARGET, got '{raw}'"
            )
        try:
            rate = float(raw_rate.strip())
        except ValueError as error:
            raise ValueError(
                f"rate item must end with a numeric RATE, got '{raw}'"
            ) from error
        rows.append(
            DiscreteHistoryRateRow(
                source_state=source_state.strip(),
                target_state=target_state.strip(),
                rate=rate,
            )
        )
    return rows


def _parse_probability_assignments(raw_rows: list[str]) -> dict[str, float]:
    probabilities: dict[str, float] = {}
    for raw in raw_rows:
        if "=" not in raw:
            raise ValueError(
                f"probability item must be in STATE=PROBABILITY form, got '{raw}'"
            )
        state, raw_probability = raw.split("=", 1)
        if not state.strip():
            raise ValueError(f"probability item must include a STATE, got '{raw}'")
        try:
            probability = float(raw_probability.strip())
        except ValueError as error:
            raise ValueError(
                f"probability item must end with a numeric PROBABILITY, got '{raw}'"
            ) from error
        probabilities[state.strip()] = probability
    return probabilities


def _parse_labelled_run(raw: str) -> tuple[str, Path]:
    if "=" not in raw:
        raise ValueError(f"run source must be in LABEL=PATH form, got '{raw}'")
    label, raw_path = raw.split("=", 1)
    if not label.strip() or not raw_path.strip():
        raise ValueError(f"run source must include both LABEL and PATH, got '{raw}'")
    return label.strip(), Path(raw_path.strip())


def _parse_time_bin_definition(raw: str) -> TimeBinDefinition:
    from bijux_phylogenetics.biogeography import TimeBinDefinition

    parts = [part.strip() for part in raw.split(":", 2)]
    if len(parts) != 3:
        raise ValueError(f"time bin must be in LABEL:START:END form, got '{raw}'")
    label, raw_start, raw_end = parts
    if not label:
        raise ValueError(f"time bin label must be non-empty, got '{raw}'")
    try:
        start_depth = float(raw_start)
        end_depth = float(raw_end)
    except ValueError as error:
        raise ValueError(
            f"time bin depths must be numeric in LABEL:START:END form, got '{raw}'"
        ) from error
    return TimeBinDefinition(
        label=label,
        start_depth=start_depth,
        end_depth=end_depth,
    )


def _validate_ancestral_discrete_model_arguments(
    args: Any, parser: argparse.ArgumentParser
) -> None:
    if (
        getattr(args, "kind", None) == "discrete"
        and getattr(args, "state_ordering", "unordered") == "ordered"
    ):
        resolved_model = getattr(args, "model", None) or "fitch"
        if resolved_model == "fitch":
            parser.error(
                "ordered ancestral discrete reconstruction requires a likelihood model"
            )


def _build_annotation_strips(table, columns: list[str]) -> list[AnnotationStrip]:
    from bijux_phylogenetics.render.tree_svg import AnnotationStrip

    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise MetadataJoinError(
            f"table does not contain columns: {', '.join(missing_columns)}"
        )
    return [
        AnnotationStrip(
            name=column,
            values={
                row[table.taxon_column]: row[column]
                for row in table.rows
                if row[column]
            },
        )
        for column in columns
    ]


def _build_numeric_trait_map(table, column: str) -> dict[str, float]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    values: dict[str, float] = {}
    for row in table.rows:
        raw = row[column]
        if not raw:
            continue
        try:
            values[row[table.taxon_column]] = float(raw)
        except ValueError as error:
            raise MetadataJoinError(
                f"column '{column}' contains a non-numeric value for taxon '{row[table.taxon_column]}'"
            ) from error
    return values


def _build_string_trait_map(table, column: str) -> dict[str, str]:
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    return {row[table.taxon_column]: row[column] for row in table.rows if row[column]}


def _json_requested(args: Any) -> bool:
    return bool(getattr(args, "json", False) or getattr(args, "format", "") == "json")


def _add_manifest_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Write a reproducibility manifest to this JSON path.",
    )


def _add_distance_tree_method_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--method",
        required=True,
        help=(
            "Distance-tree method. Supported: neighbor-joining, bionj, upgma, wpgma, single-linkage, complete-linkage."
        ),
    )


def _add_missing_distance_policy_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--missing-distance-policy",
        choices=("reject", "mean-impute", "nearest-valid", "triangle-bound"),
        default="reject",
        help=(
            "Policy used when one or more pairwise distances are missing before distance analysis."
        ),
    )


def _add_ultrametric_tolerance_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Absolute tolerance used when comparing the two largest distances in each taxon triple.",
    )


def _add_external_adapter_execution_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_resume: bool = True,
) -> None:
    if include_resume:
        parser.add_argument(
            "--resume",
            action="store_true",
            help="Reuse one completed governed engine run when the manifest, inputs, and outputs still match.",
        )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        help="Stop the governed engine execution if it exceeds this wall-clock budget.",
    )
    parser.add_argument(
        "--incomplete-run-policy",
        choices=("reject", "clean"),
        default="reject",
        help="Reject or clean incomplete governed engine outputs before starting a new run.",
    )


def _add_preflight_executable_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mafft-executable", type=str)
    parser.add_argument("--trimal-executable", type=str)
    parser.add_argument("--iqtree-executable", type=str)
    parser.add_argument("--fasttree-executable", type=str)
    parser.add_argument("--mrbayes-executable", type=str)
    parser.add_argument("--beast-executable", type=str)


def _adapter_version_args(engine_name: str) -> tuple[str, ...]:
    normalized = engine_name.lower()
    if normalized == "fasttree":
        return ("-help",)
    if normalized == "mrbayes":
        return ("-v",)
    return ("--version",)
