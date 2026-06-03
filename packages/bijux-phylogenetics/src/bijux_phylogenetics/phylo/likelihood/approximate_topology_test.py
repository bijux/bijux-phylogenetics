from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import json
import math
from pathlib import Path
import random

import numpy

from bijux_phylogenetics.io.newick import loads_newick, write_newick_tree_set
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.candidate_tree_site_likelihood_matrix import (
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix,
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    ApproximateTopologyTestReport,
    ApproximateTopologyTestResamplingRow,
    ApproximateTopologyTestSummaryRow,
    CandidateTreeSiteLikelihoodMatrixReport,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_APPROXIMATE_TOPOLOGY_TEST_CAUTION_LABEL = "site-resampled p-like statistics are approximate ranking aids and must not be interpreted as AU/SH-style p-values"
_APPROXIMATE_TOPOLOGY_TEST_RESAMPLING_METHOD = "site-resampling-with-replacement"


def evaluate_nucleotide_approximate_topology_test(
    trees: list[PhyloTree] | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    resampling_replicate_count: int = 100,
    resampling_seed: int = 1,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> ApproximateTopologyTestReport:
    """Compare candidate topologies with one site-resampled approximate likelihood test."""
    matrix_report = evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
        trees,
        records,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )
    return evaluate_nucleotide_approximate_topology_test_from_matrix_report(
        matrix_report,
        resampling_replicate_count=resampling_replicate_count,
        resampling_seed=resampling_seed,
    )


def evaluate_nucleotide_approximate_topology_test_from_alignment(
    tree_set_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    resampling_replicate_count: int = 100,
    resampling_seed: int = 1,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
) -> ApproximateTopologyTestReport:
    """Compare candidate topologies from tree-set and alignment paths."""
    matrix_report = (
        evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
            tree_set_path,
            alignment_path,
            model_name=model_name,
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
        )
    )
    return evaluate_nucleotide_approximate_topology_test_from_matrix_report(
        matrix_report,
        resampling_replicate_count=resampling_replicate_count,
        resampling_seed=resampling_seed,
    )


def evaluate_nucleotide_approximate_topology_test_from_matrix_report(
    matrix_report: CandidateTreeSiteLikelihoodMatrixReport,
    *,
    resampling_replicate_count: int = 100,
    resampling_seed: int = 1,
) -> ApproximateTopologyTestReport:
    """Compare candidate topologies from one precomputed site-likelihood matrix report."""
    validated_replicate_count = validate_approximate_topology_test_replicate_count(
        resampling_replicate_count
    )
    candidate_vectors = resolve_candidate_tree_site_likelihood_vectors(matrix_report)
    observed_best = select_observed_best_candidate_tree(candidate_vectors)
    rng = random.Random(resampling_seed)  # nosec B311
    resampling_rows: list[ApproximateTopologyTestResamplingRow] = []
    rows_by_candidate_id: dict[str, list[ApproximateTopologyTestResamplingRow]] = (
        defaultdict(list)
    )
    site_count = matrix_report.site_count

    for replicate_index in range(1, validated_replicate_count + 1):
        sampled_indices = [rng.randrange(site_count) for _ in range(site_count)]
        resampled_totals = {
            candidate.candidate_tree_id: sum(
                candidate.site_log_likelihoods[index] for index in sampled_indices
            )
            for candidate in candidate_vectors
        }
        observed_best_resampled_log_likelihood = resampled_totals[
            observed_best.candidate_tree_id
        ]
        for candidate in candidate_vectors:
            resampled_log_likelihood = resampled_totals[candidate.candidate_tree_id]
            resampled_delta_log_likelihood = (
                observed_best_resampled_log_likelihood - resampled_log_likelihood
            )
            candidate_matches_or_beats_observed_best = (
                resampled_log_likelihood > observed_best_resampled_log_likelihood
                or math.isclose(
                    resampled_log_likelihood,
                    observed_best_resampled_log_likelihood,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )
            )
            row = ApproximateTopologyTestResamplingRow(
                replicate_index=replicate_index,
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                resampled_log_likelihood=resampled_log_likelihood,
                observed_best_tree_id=observed_best.candidate_tree_id,
                observed_best_tree_label=observed_best.candidate_tree_label,
                observed_best_resampled_log_likelihood=observed_best_resampled_log_likelihood,
                resampled_delta_log_likelihood=resampled_delta_log_likelihood,
                candidate_matches_or_beats_observed_best=candidate_matches_or_beats_observed_best,
            )
            resampling_rows.append(row)
            rows_by_candidate_id[candidate.candidate_tree_id].append(row)

    summary_rows = [
        ApproximateTopologyTestSummaryRow(
            candidate_tree_id=candidate.candidate_tree_id,
            candidate_tree_label=candidate.candidate_tree_label,
            tree_newick=candidate.tree_newick,
            observed_log_likelihood=candidate.log_likelihood,
            observed_delta_log_likelihood=observed_best.log_likelihood
            - candidate.log_likelihood,
            observed_best_tree=candidate.candidate_tree_id
            == observed_best.candidate_tree_id,
            resampling_win_count=sum(
                1
                for row in rows_by_candidate_id[candidate.candidate_tree_id]
                if row.candidate_matches_or_beats_observed_best
            ),
            resampling_frequency=sum(
                1.0
                for row in rows_by_candidate_id[candidate.candidate_tree_id]
                if row.candidate_matches_or_beats_observed_best
            )
            / float(validated_replicate_count),
            p_like_statistic=sum(
                1.0
                for row in rows_by_candidate_id[candidate.candidate_tree_id]
                if row.candidate_matches_or_beats_observed_best
            )
            / float(validated_replicate_count),
            resampling_mean_delta_log_likelihood=_mean(
                [
                    row.resampled_delta_log_likelihood
                    for row in rows_by_candidate_id[candidate.candidate_tree_id]
                ]
            ),
            resampling_min_delta_log_likelihood=min(
                row.resampled_delta_log_likelihood
                for row in rows_by_candidate_id[candidate.candidate_tree_id]
            ),
            resampling_max_delta_log_likelihood=max(
                row.resampled_delta_log_likelihood
                for row in rows_by_candidate_id[candidate.candidate_tree_id]
            ),
            caution_label=_APPROXIMATE_TOPOLOGY_TEST_CAUTION_LABEL,
        )
        for candidate in candidate_vectors
    ]
    return ApproximateTopologyTestReport(
        algorithm="nucleotide-approximate-topology-test",
        model_name=matrix_report.model_name,
        tree_set_path=matrix_report.tree_set_path,
        alignment_path=matrix_report.alignment_path,
        taxa=matrix_report.taxa,
        tree_count=matrix_report.tree_count,
        site_count=matrix_report.site_count,
        pattern_count=matrix_report.pattern_count,
        compression_used=matrix_report.compression_used,
        expansion_policy=matrix_report.expansion_policy,
        parameter_values=matrix_report.parameter_values,
        resampling_method=_APPROXIMATE_TOPOLOGY_TEST_RESAMPLING_METHOD,
        resampling_replicate_count=validated_replicate_count,
        resampling_seed=resampling_seed,
        observed_best_tree_id=observed_best.candidate_tree_id,
        observed_best_tree_label=observed_best.candidate_tree_label,
        caution_label=_APPROXIMATE_TOPOLOGY_TEST_CAUTION_LABEL,
        summary_rows=summary_rows,
        resampling_rows=resampling_rows,
    )


def validate_approximate_topology_test_replicate_count(
    resampling_replicate_count: int,
) -> int:
    """Require at least one site-resampling replicate for approximate topology tests."""
    if resampling_replicate_count < 1:
        raise ValueError("resampling_replicate_count must be at least one")
    return resampling_replicate_count


def resolve_candidate_tree_site_likelihood_vectors(
    matrix_report: CandidateTreeSiteLikelihoodMatrixReport,
) -> list[CandidateTreeSiteLikelihoodVector]:
    """Project one expanded candidate-tree matrix report into ordered site vectors."""
    site_vectors_by_candidate_id: dict[str, list[float | None]] = {
        row.candidate_tree_id: [None] * matrix_report.site_count
        for row in matrix_report.candidate_trees
    }
    for row in matrix_report.matrix_rows:
        site_index = row.site_position - 1
        if not (0 <= site_index < matrix_report.site_count):
            raise ValueError(
                f"site_position {row.site_position} is outside the declared site_count"
            )
        candidate_site_vector = site_vectors_by_candidate_id.get(row.candidate_tree_id)
        if candidate_site_vector is None:
            raise ValueError(
                f"matrix row references unknown candidate_tree_id '{row.candidate_tree_id}'"
            )
        if candidate_site_vector[site_index] is not None:
            raise ValueError(
                f"candidate tree '{row.candidate_tree_id}' has duplicate site_position {row.site_position}"
            )
        candidate_site_vector[site_index] = row.log_likelihood
    candidate_vectors: list[CandidateTreeSiteLikelihoodVector] = []
    for candidate in matrix_report.candidate_trees:
        candidate_site_vector = site_vectors_by_candidate_id[
            candidate.candidate_tree_id
        ]
        if any(value is None for value in candidate_site_vector):
            raise ValueError(
                f"candidate tree '{candidate.candidate_tree_id}' is missing one or more site log likelihood rows"
            )
        site_log_likelihoods = tuple(float(value) for value in candidate_site_vector)
        candidate_vectors.append(
            CandidateTreeSiteLikelihoodVector(
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                tree_newick=candidate.tree_newick,
                log_likelihood=candidate.log_likelihood,
                site_log_likelihoods=site_log_likelihoods,
            )
        )
    return candidate_vectors


def select_observed_best_candidate_tree(
    candidates: list[CandidateTreeSiteLikelihoodVector],
) -> CandidateTreeSiteLikelihoodVector:
    """Choose the observed best candidate by likelihood with deterministic tie breaks."""
    if not candidates:
        raise ValueError(
            "approximate topology test requires at least one candidate tree"
        )
    best_candidate = candidates[0]
    for candidate in candidates[1:]:
        if prefer_higher_likelihood_candidate(candidate, best_candidate):
            best_candidate = candidate
    return best_candidate


def prefer_higher_likelihood_candidate(
    left: CandidateTreeSiteLikelihoodVector,
    right: CandidateTreeSiteLikelihoodVector,
) -> bool:
    """Prefer higher likelihoods, then canonical tree text, across candidate trees."""
    if left.log_likelihood > right.log_likelihood and not math.isclose(
        left.log_likelihood,
        right.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        return True
    if right.log_likelihood > left.log_likelihood and not math.isclose(
        left.log_likelihood,
        right.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        return False
    if left.tree_newick != right.tree_newick:
        return left.tree_newick < right.tree_newick
    return left.candidate_tree_id < right.candidate_tree_id


def _mean(values: list[float]) -> float:
    if not values:
        raise ValueError("mean requires at least one value")
    return sum(values) / float(len(values))


@dataclass(frozen=True, slots=True)
class CandidateTreeSiteLikelihoodVector:
    """One candidate tree with a complete ordered vector of site log likelihoods."""

    candidate_tree_id: str
    candidate_tree_label: str
    tree_newick: str
    log_likelihood: float
    site_log_likelihoods: tuple[float, ...]


def write_approximate_topology_test_summary_table(
    path: Path,
    report: ApproximateTopologyTestReport,
) -> Path:
    """Write one candidate-tree approximate topology test summary table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "candidate_tree_id",
        "candidate_tree_label",
        "observed_log_likelihood",
        "observed_delta_log_likelihood",
        "observed_best_tree",
        "resampling_win_count",
        "resampling_frequency",
        "p_like_statistic",
        "resampling_mean_delta_log_likelihood",
        "resampling_min_delta_log_likelihood",
        "resampling_max_delta_log_likelihood",
        "caution_label",
        "tree_newick",
    ]
    rows = ["\t".join(columns)]
    for row in report.summary_rows:
        rows.append(
            "\t".join(
                [
                    row.candidate_tree_id,
                    row.candidate_tree_label,
                    repr(row.observed_log_likelihood),
                    repr(row.observed_delta_log_likelihood),
                    str(row.observed_best_tree),
                    str(row.resampling_win_count),
                    repr(row.resampling_frequency),
                    repr(row.p_like_statistic),
                    repr(row.resampling_mean_delta_log_likelihood),
                    repr(row.resampling_min_delta_log_likelihood),
                    repr(row.resampling_max_delta_log_likelihood),
                    row.caution_label,
                    row.tree_newick,
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_approximate_topology_test_resampling_table(
    path: Path,
    report: ApproximateTopologyTestReport,
) -> Path:
    """Write one resampling-distribution table for an approximate topology test."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "replicate_index",
        "candidate_tree_id",
        "candidate_tree_label",
        "resampled_log_likelihood",
        "observed_best_tree_id",
        "observed_best_tree_label",
        "observed_best_resampled_log_likelihood",
        "resampled_delta_log_likelihood",
        "candidate_matches_or_beats_observed_best",
    ]
    rows = ["\t".join(columns)]
    for row in report.resampling_rows:
        rows.append(
            "\t".join(
                [
                    str(row.replicate_index),
                    row.candidate_tree_id,
                    row.candidate_tree_label,
                    repr(row.resampled_log_likelihood),
                    row.observed_best_tree_id,
                    row.observed_best_tree_label,
                    repr(row.observed_best_resampled_log_likelihood),
                    repr(row.resampled_delta_log_likelihood),
                    str(row.candidate_matches_or_beats_observed_best),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_approximate_topology_test_run_json(
    path: Path,
    report: ApproximateTopologyTestReport,
) -> Path:
    """Write one machine-readable approximate topology test payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "model_name": report.model_name,
        "tree_set_path": report.tree_set_path,
        "alignment_path": report.alignment_path,
        "taxa": report.taxa,
        "tree_count": report.tree_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "compression_used": report.compression_used,
        "expansion_policy": report.expansion_policy,
        "parameter_values": report.parameter_values,
        "resampling_method": report.resampling_method,
        "resampling_replicate_count": report.resampling_replicate_count,
        "resampling_seed": report.resampling_seed,
        "observed_best_tree_id": report.observed_best_tree_id,
        "observed_best_tree_label": report.observed_best_tree_label,
        "caution_label": report.caution_label,
        "summary_rows": [
            {
                "candidate_tree_id": row.candidate_tree_id,
                "candidate_tree_label": row.candidate_tree_label,
                "tree_newick": row.tree_newick,
                "observed_log_likelihood": row.observed_log_likelihood,
                "observed_delta_log_likelihood": row.observed_delta_log_likelihood,
                "observed_best_tree": row.observed_best_tree,
                "resampling_win_count": row.resampling_win_count,
                "resampling_frequency": row.resampling_frequency,
                "p_like_statistic": row.p_like_statistic,
                "resampling_mean_delta_log_likelihood": row.resampling_mean_delta_log_likelihood,
                "resampling_min_delta_log_likelihood": row.resampling_min_delta_log_likelihood,
                "resampling_max_delta_log_likelihood": row.resampling_max_delta_log_likelihood,
                "caution_label": row.caution_label,
            }
            for row in report.summary_rows
        ],
        "resampling_rows": [
            {
                "replicate_index": row.replicate_index,
                "candidate_tree_id": row.candidate_tree_id,
                "candidate_tree_label": row.candidate_tree_label,
                "resampled_log_likelihood": row.resampled_log_likelihood,
                "observed_best_tree_id": row.observed_best_tree_id,
                "observed_best_tree_label": row.observed_best_tree_label,
                "observed_best_resampled_log_likelihood": row.observed_best_resampled_log_likelihood,
                "resampled_delta_log_likelihood": row.resampled_delta_log_likelihood,
                "candidate_matches_or_beats_observed_best": row.candidate_matches_or_beats_observed_best,
            }
            for row in report.resampling_rows
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_approximate_topology_test_artifacts(
    out_dir: Path,
    report: ApproximateTopologyTestReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one approximate topology test."""
    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_tree_path = write_newick_tree_set(
        out_dir / "candidate_trees.nwk",
        [loads_newick(row.tree_newick) for row in report.summary_rows],
    )
    summary_path = write_approximate_topology_test_summary_table(
        out_dir / "topology_test_summary.tsv",
        report,
    )
    resampling_path = write_approximate_topology_test_resampling_table(
        out_dir / "resampling_distribution.tsv",
        report,
    )
    run_json_path = write_approximate_topology_test_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "candidate_tree_path": candidate_tree_path,
        "summary_path": summary_path,
        "resampling_path": resampling_path,
        "run_json_path": run_json_path,
    }
