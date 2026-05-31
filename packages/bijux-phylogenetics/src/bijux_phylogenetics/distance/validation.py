from __future__ import annotations

import math

from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades

from .genetic_distance_matrix import compute_pairwise_genetic_distance_matrix
from .imported import build_tree_from_imported_distance_matrix
from .models import (
    DistanceReferenceObservation,
    DistanceReferenceValidationReport,
    DistanceTreeReferenceObservation,
)
from .shared import _PACKAGE_ROOT

_REFERENCE_DISTANCE_CASES = (
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_distance.fasta",
        "case": "dna-p-distance",
        "model": "p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "B",
        "expected": 0.125,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_distance.fasta",
        "case": "dna-jukes-cantor",
        "model": "jukes-cantor",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "B",
        "expected": 0.136741167595466,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_distance.fasta",
        "case": "dna-kimura-2-parameter",
        "model": "kimura-2-parameter",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "B",
        "expected": 0.14384103622589,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_distance_gaps.fasta",
        "case": "dna-felsenstein-81",
        "model": "felsenstein-81",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "B",
        "expected": 0.189450794670797,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_duplicates.fasta",
        "case": "dna-tamura-nei-93",
        "model": "tamura-nei-93",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "C",
        "expected": 0.182459298648674,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_protein.fasta",
        "case": "protein-p-distance",
        "model": "amino-acid-p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "P1",
        "right": "P2",
        "expected": 0.125,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_distance_gaps.fasta",
        "case": "gap-complete-deletion",
        "model": "p-distance",
        "gap_handling": "complete-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "D",
        "expected": 0.0,
        "expected_ambiguity_sites": 0,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
        "case": "ambiguity-ignore",
        "model": "p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "ignore",
        "left": "A",
        "right": "B",
        "expected": 0.0,
        "expected_ambiguity_sites": 1,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
        "case": "ambiguity-partial-match",
        "model": "p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "partial-match",
        "left": "A",
        "right": "B",
        "expected": 0.15,
        "expected_ambiguity_sites": 1,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
        "case": "ambiguity-strict-mismatch",
        "model": "p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "strict-mismatch",
        "left": "A",
        "right": "B",
        "expected": 0.2,
        "expected_ambiguity_sites": 1,
    },
    {
        "path": _PACKAGE_ROOT
        / "tests/fixtures/alignments/example_alignment_ambiguity.fasta",
        "case": "ambiguity-report-only",
        "model": "p-distance",
        "gap_handling": "pairwise-deletion",
        "ambiguity_policy": "report-only",
        "left": "A",
        "right": "B",
        "expected": 0.0,
        "expected_ambiguity_sites": 1,
    },
)


def _reference_distance_observations() -> list[DistanceReferenceObservation]:
    observations: list[DistanceReferenceObservation] = []
    for example in _REFERENCE_DISTANCE_CASES:
        report = compute_pairwise_genetic_distance_matrix(
            example["path"],
            model=str(example["model"]),
            gap_handling=str(example["gap_handling"]),
            ambiguity_policy=str(example["ambiguity_policy"]),
        )
        pair = next(
            row
            for row in report.pairs
            if row.left_identifier == example["left"]
            and row.right_identifier == example["right"]
        )
        expected_distance = float(example["expected"])
        observed = pair.distance
        passed = observed is not None and math.isclose(
            observed, expected_distance, rel_tol=1e-12, abs_tol=1e-12
        )
        observations.append(
            DistanceReferenceObservation(
                case=str(example["case"]),
                model=str(example["model"]),
                gap_handling=str(example["gap_handling"]),
                ambiguity_policy=str(example["ambiguity_policy"]),
                left_identifier=str(example["left"]),
                right_identifier=str(example["right"]),
                expected_distance=expected_distance,
                observed_distance=observed,
                comparable_sites=pair.comparable_sites,
                expected_ambiguity_sites=int(example["expected_ambiguity_sites"]),
                observed_ambiguity_sites=pair.ambiguity_sites,
                passed=passed
                and pair.ambiguity_sites == int(example["expected_ambiguity_sites"]),
            )
        )
    return observations


def _reference_tree_observations() -> list[DistanceTreeReferenceObservation]:
    observations: list[DistanceTreeReferenceObservation] = []
    matrix_path = (
        _PACKAGE_ROOT
        / "tests/fixtures/metadata/example_distance_matrix_ultrametric.tsv"
    )
    nj_tree, _ = build_tree_from_imported_distance_matrix(
        matrix_path, method="neighbor-joining"
    )
    nj_observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in informative_rooted_clades(nj_tree, set(nj_tree.tip_names))
    )
    observations.append(
        DistanceTreeReferenceObservation(
            case="neighbor-joining-reference-clustering",
            method="neighbor-joining",
            matrix_path=matrix_path,
            expected_clades=["A|B"],
            observed_clades=nj_observed_clades,
            passed=nj_observed_clades == ["A|B"],
        )
    )

    upgma_tree, _ = build_tree_from_imported_distance_matrix(
        matrix_path, method="upgma"
    )
    upgma_observed_clades = sorted(
        "|".join(sorted(clade))
        for clade in informative_rooted_clades(upgma_tree, set(upgma_tree.tip_names))
    )
    observations.append(
        DistanceTreeReferenceObservation(
            case="upgma-ultrametric-clustering",
            method="upgma",
            matrix_path=matrix_path,
            expected_clades=["A|B", "C|D"],
            observed_clades=upgma_observed_clades,
            passed=upgma_observed_clades == ["A|B", "C|D"],
        )
    )
    return observations


def validate_distance_reference_examples() -> DistanceReferenceValidationReport:
    """Validate core distance examples against durable reference expectations."""
    observations = _reference_distance_observations()
    tree_observations = _reference_tree_observations()
    return DistanceReferenceValidationReport(
        observations=observations,
        tree_observations=tree_observations,
        all_passed=all(observation.passed for observation in observations)
        and all(observation.passed for observation in tree_observations),
    )
