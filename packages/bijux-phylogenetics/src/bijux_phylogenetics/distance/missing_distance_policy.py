from __future__ import annotations

import math
from statistics import fmean

from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .models import (
    MissingDistanceImputation,
    MissingDistancePolicy,
    MissingDistancePolicyReport,
)
from .shared import _pair_key

SUPPORTED_MISSING_DISTANCE_POLICIES = (
    "reject",
    "mean-impute",
    "nearest-valid",
    "triangle-bound",
)


def normalize_missing_distance_policy(
    policy: MissingDistancePolicy,
) -> MissingDistancePolicy:
    """Normalize and validate one missing-distance policy name."""
    normalized = str(policy).strip().lower()
    if normalized not in SUPPORTED_MISSING_DISTANCE_POLICIES:
        raise ValueError(f"unsupported missing-distance policy: {policy}")
    return normalized


def apply_missing_distance_policy(
    identifiers: list[str],
    known_pair_distances: dict[tuple[str, str], float],
    *,
    policy: MissingDistancePolicy = "reject",
) -> tuple[dict[tuple[str, str], float], MissingDistancePolicyReport]:
    """Resolve missing off-diagonal distances under one explicit policy."""
    normalized_policy = normalize_missing_distance_policy(policy)
    pair_distances = _normalized_pair_distances(known_pair_distances)
    missing_pairs = _missing_pairs(identifiers, pair_distances)
    report = MissingDistancePolicyReport(
        policy=normalized_policy,
        taxon_count=len(identifiers),
        requested_pair_count=(len(identifiers) * (len(identifiers) - 1)) // 2,
        missing_pairs=[f"{left}/{right}" for left, right in missing_pairs],
        imputed_rows=[],
        unresolved_pairs=[],
        warnings=[],
    )
    if not missing_pairs:
        return _full_distance_lookup(identifiers, pair_distances), report
    if normalized_policy == "reject":
        raise InvalidDistanceMatrixError(
            "missing-distance policy 'reject' blocks incomplete distance pairs: "
            + ", ".join(report.missing_pairs),
            code="missing_distance_policy_rejects_incomplete_pairs",
            details={
                "policy": normalized_policy,
                "missing_pairs": report.missing_pairs,
            },
        )

    resolved_distances = dict(pair_distances)
    for left_identifier, right_identifier in missing_pairs:
        pair_key = _pair_key(left_identifier, right_identifier)
        imputed_distance, rationale = _imputed_distance(
            pair_distances=resolved_distances,
            identifiers=identifiers,
            left_identifier=left_identifier,
            right_identifier=right_identifier,
            policy=normalized_policy,
        )
        if imputed_distance is None:
            report.unresolved_pairs.append(f"{left_identifier}/{right_identifier}")
            continue
        resolved_distances[pair_key] = round(imputed_distance, 15)
        report.imputed_rows.append(
            MissingDistanceImputation(
                left_identifier=left_identifier,
                right_identifier=right_identifier,
                imputed_distance=round(imputed_distance, 15),
                policy=normalized_policy,
                rationale=rationale,
            )
        )
    if report.unresolved_pairs:
        raise InvalidDistanceMatrixError(
            "missing-distance policy could not impute all incomplete distance pairs: "
            + ", ".join(report.unresolved_pairs),
            code="missing_distance_policy_cannot_impute_pairs",
            details={
                "policy": normalized_policy,
                "unresolved_pairs": report.unresolved_pairs,
            },
        )
    report.warnings.append(
        "one or more missing pairwise distances were imputed before distance analysis"
    )
    return _full_distance_lookup(identifiers, resolved_distances), report


def _normalized_pair_distances(
    known_pair_distances: dict[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    normalized: dict[tuple[str, str], float] = {}
    for (left_identifier, right_identifier), distance in known_pair_distances.items():
        if left_identifier == right_identifier:
            if not math.isclose(float(distance), 0.0, abs_tol=1e-12, rel_tol=0.0):
                raise InvalidDistanceMatrixError(
                    "distance matrix has nonzero diagonal entries",
                    code="missing_distance_policy_nonzero_diagonal",
                )
            continue
        normalized_distance = float(distance)
        if not math.isfinite(normalized_distance):
            raise InvalidDistanceMatrixError(
                f"distance matrix contains non-finite distance for {left_identifier}/{right_identifier}",
                code="missing_distance_policy_nonfinite_distance",
            )
        if normalized_distance < 0.0:
            raise InvalidDistanceMatrixError(
                "distance matrix contains negative distances",
                code="missing_distance_policy_negative_distance",
            )
        pair_key = _pair_key(left_identifier, right_identifier)
        previous = normalized.get(pair_key)
        if previous is not None and not math.isclose(
            previous,
            normalized_distance,
            abs_tol=1e-12,
            rel_tol=1e-12,
        ):
            raise InvalidDistanceMatrixError(
                "distance matrix contains asymmetric directional entries",
                code="missing_distance_policy_asymmetric_inputs",
                details={
                    "left_identifier": pair_key[0],
                    "right_identifier": pair_key[1],
                    "first_distance": previous,
                    "second_distance": normalized_distance,
                },
            )
        normalized[pair_key] = normalized_distance
    return normalized


def _missing_pairs(
    identifiers: list[str],
    pair_distances: dict[tuple[str, str], float],
) -> list[tuple[str, str]]:
    missing: list[tuple[str, str]] = []
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            if _pair_key(left_identifier, right_identifier) not in pair_distances:
                missing.append((left_identifier, right_identifier))
    return missing


def _full_distance_lookup(
    identifiers: list[str],
    pair_distances: dict[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    lookup = {(identifier, identifier): 0.0 for identifier in identifiers}
    for (left_identifier, right_identifier), distance in pair_distances.items():
        lookup[(left_identifier, right_identifier)] = distance
        lookup[(right_identifier, left_identifier)] = distance
    return lookup


def _imputed_distance(
    *,
    pair_distances: dict[tuple[str, str], float],
    identifiers: list[str],
    left_identifier: str,
    right_identifier: str,
    policy: MissingDistancePolicy,
) -> tuple[float | None, str]:
    if policy == "mean-impute":
        defined_distances = list(pair_distances.values())
        if not defined_distances:
            return (
                None,
                "no defined pairwise distances were available for mean imputation",
            )
        return (
            fmean(defined_distances),
            f"mean of {len(defined_distances)} defined pairwise distances",
        )
    if policy == "nearest-valid":
        adjacent = [
            distance
            for (pair_left, pair_right), distance in pair_distances.items()
            if left_identifier in {pair_left, pair_right}
            or right_identifier in {pair_left, pair_right}
        ]
        if not adjacent:
            return None, "no adjacent defined pairwise distances were available"
        return min(adjacent), "nearest defined distance touching either endpoint"
    triangle_bounds = [
        pair_distances[_pair_key(left_identifier, middle_identifier)]
        + pair_distances[_pair_key(middle_identifier, right_identifier)]
        for middle_identifier in identifiers
        if middle_identifier not in {left_identifier, right_identifier}
        and _pair_key(left_identifier, middle_identifier) in pair_distances
        and _pair_key(middle_identifier, right_identifier) in pair_distances
    ]
    if not triangle_bounds:
        return None, "no complete triangle path was available for imputation"
    return min(triangle_bounds), "tightest available triangle upper bound"
