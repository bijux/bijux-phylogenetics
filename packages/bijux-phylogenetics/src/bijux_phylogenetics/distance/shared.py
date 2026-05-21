from __future__ import annotations

import hashlib
from pathlib import Path

from bijux_phylogenetics.runtime.errors import UnsupportedDistanceTreeMethodError

from .models import (
    DistanceModel,
    DistanceTreeMethodPolicy,
    GeneticDistanceMatrix,
    PairwiseGeneticDistance,
    UPGMAUltrametricViolation,
)

_DISTANCE_MODEL_ALIASES = {
    "raw": "p-distance",
    "p-distance": "p-distance",
    "jc69": "jukes-cantor",
    "jukes-cantor": "jukes-cantor",
    "k80": "kimura-2-parameter",
    "kimura-2-parameter": "kimura-2-parameter",
    "f81": "felsenstein-81",
    "felsenstein-81": "felsenstein-81",
    "tn93": "tamura-nei-93",
    "tamura-nei-93": "tamura-nei-93",
    "amino-acid-p-distance": "amino-acid-p-distance",
}
_PACKAGE_ROOT = Path(__file__).resolve().parents[3]


def _normalize_distance_model(model: DistanceModel) -> DistanceModel:
    normalized = _DISTANCE_MODEL_ALIASES.get(model)
    if normalized is None:
        raise ValueError(f"unsupported distance model: {model}")
    return normalized


def list_distance_tree_method_policies() -> list[DistanceTreeMethodPolicy]:
    """Return the governed distance-tree method support policy surface."""
    return [
        DistanceTreeMethodPolicy(
            method="neighbor-joining",
            supported=True,
            reference_surface="ape::nj",
            support_scope="owned-runtime-and-live-ape-parity",
            summary="Neighbor-Joining is fully supported and covered by the governed live ape parity lane.",
            limitations=[
                "neighbor-joining remains a distance-summary method rather than a full likelihood inference",
            ],
        ),
        DistanceTreeMethodPolicy(
            method="upgma",
            supported=True,
            reference_surface=None,
            support_scope="owned-runtime",
            summary="UPGMA is supported as an owned runtime method for ultrametric reviewer workflows.",
            limitations=[
                "upgma assumes an ultrametric clock-like process and can misplace taxa when rates vary among lineages",
            ],
        ),
        DistanceTreeMethodPolicy(
            method="bionj",
            supported=False,
            reference_surface="ape::bionj",
            support_scope="explicitly-excluded",
            summary="BIONJ is explicitly out of scope for this round, and ape::bionj does not have an owned Bijux runtime or governed ape parity lane.",
            limitations=[
                "request neighbor-joining when a governed ape-parity distance-tree workflow is required in this round",
            ],
        ),
    ]


def resolve_distance_tree_method_policy(method: str) -> DistanceTreeMethodPolicy:
    """Resolve one distance-tree method name to its owned support policy."""
    normalized = method.strip().lower()
    for policy in list_distance_tree_method_policies():
        if policy.method == normalized:
            return policy
    supported = [
        policy.method
        for policy in list_distance_tree_method_policies()
        if policy.supported
    ]
    raise UnsupportedDistanceTreeMethodError(
        f"unsupported tree-building method '{method}'; supported methods are {', '.join(supported)}",
        details={
            "requested_method": method,
            "supported_methods": supported,
        },
    )


def _require_supported_distance_tree_method(method: str) -> DistanceTreeMethodPolicy:
    """Resolve one distance-tree method and reject explicitly excluded methods."""
    policy = resolve_distance_tree_method_policy(method)
    if policy.supported:
        return policy
    supported = [
        row.method for row in list_distance_tree_method_policies() if row.supported
    ]
    raise UnsupportedDistanceTreeMethodError(
        policy.summary,
        details={
            "requested_method": policy.method,
            "supported_methods": supported,
            "excluded_method": policy.method,
            "reference_surface": policy.reference_surface,
        },
    )


def _allowed_models_for_alphabet(alphabet: str) -> set[str]:
    if alphabet in {"dna", "rna"}:
        return {
            "p-distance",
            "jukes-cantor",
            "kimura-2-parameter",
            "felsenstein-81",
            "tamura-nei-93",
        }
    if alphabet == "protein":
        return {"amino-acid-p-distance"}
    return set()


def _pair_key(left_identifier: str, right_identifier: str) -> tuple[str, str]:
    return tuple(sorted((left_identifier, right_identifier)))


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _unique_genetic_distance_pairs(
    report: GeneticDistanceMatrix,
) -> list[PairwiseGeneticDistance]:
    unique: dict[tuple[str, str], PairwiseGeneticDistance] = {}
    for pair in report.pairs:
        if pair.left_identifier == pair.right_identifier:
            continue
        unique.setdefault(_pair_key(pair.left_identifier, pair.right_identifier), pair)
    return [unique[key] for key in sorted(unique)]


def _iter_ultrametric_violations(
    identifiers: list[str],
    distances: dict[tuple[str, str], float],
    *,
    tolerance: float,
) -> list[UPGMAUltrametricViolation]:
    violations: list[UPGMAUltrametricViolation] = []
    for left_index, left_identifier in enumerate(identifiers):
        for middle_index in range(left_index + 1, len(identifiers)):
            middle_identifier = identifiers[middle_index]
            for right_index in range(middle_index + 1, len(identifiers)):
                right_identifier = identifiers[right_index]
                pair_keys = [
                    _pair_key(left_identifier, middle_identifier),
                    _pair_key(left_identifier, right_identifier),
                    _pair_key(middle_identifier, right_identifier),
                ]
                if any(pair_key not in distances for pair_key in pair_keys):
                    continue
                triple = [distances[pair_key] for pair_key in pair_keys]
                ordered = sorted(triple)
                deviation = abs(ordered[2] - ordered[1])
                if deviation > tolerance:
                    violations.append(
                        UPGMAUltrametricViolation(
                            left_identifier=left_identifier,
                            middle_identifier=middle_identifier,
                            right_identifier=right_identifier,
                            smallest_distance=ordered[0],
                            middle_distance=ordered[1],
                            largest_distance=ordered[2],
                            deviation=deviation,
                        )
                    )
    return sorted(
        violations,
        key=lambda row: (
            row.left_identifier,
            row.middle_identifier,
            row.right_identifier,
        ),
    )
