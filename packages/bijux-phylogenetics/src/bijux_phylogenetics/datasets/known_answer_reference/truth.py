from __future__ import annotations

import csv
from pathlib import Path

from .models import (
    KnownAnswerContinuousNodeTruth,
    KnownAnswerDiscreteNodeTruth,
    KnownAnswerRecoveryThreshold,
    KnownAnswerTransitionTruth,
)


def load_true_parameter_map(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["parameter"]: row["value"] for row in reader}


def load_true_continuous_nodes(path: Path) -> list[KnownAnswerContinuousNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerContinuousNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_value=float(row["true_value"]),
            )
            for row in reader
        ]


def load_true_discrete_nodes(path: Path) -> list[KnownAnswerDiscreteNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerDiscreteNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_state=row["true_state"],
            )
            for row in reader
        ]


def load_true_transition_rows(path: Path) -> list[KnownAnswerTransitionTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerTransitionTruth(
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                branch_length=float(row["branch_length"]),
                source_state=row["source_state"],
                target_state=row["target_state"],
                changed=row["changed"].strip().lower() == "true",
                event_count=int(row["event_count"]),
            )
            for row in reader
        ]


def load_recovery_thresholds(path: Path) -> list[KnownAnswerRecoveryThreshold]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerRecoveryThreshold(
                metric=row["metric"],
                comparator=row["comparator"],
                threshold=row["threshold"],
                rationale=row["rationale"],
            )
            for row in reader
        ]


def _split_descendant_taxa(value: str) -> list[str]:
    if not value:
        return []
    return value.split(",")
