from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.approximate_topology_test import (
    CandidateTreeSiteLikelihoodVector,
    prefer_higher_likelihood_candidate,
    resolve_candidate_tree_site_likelihood_vectors,
    select_observed_best_candidate_tree,
)
from bijux_phylogenetics.phylo.likelihood.candidate_tree_site_likelihood_matrix import (
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideShLikeBranchSupportLocalTopologyRow,
    NucleotideShLikeBranchSupportReport,
    NucleotideShLikeBranchSupportResamplingRow,
    NucleotideShLikeBranchSupportRow,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    RootedNniMoveCandidate,
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
    validate_rooted_nni_tree,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SH_LIKE_BRANCH_SUPPORT_RESAMPLING_METHOD = "site-resampling-with-replacement"
_SH_LIKE_BRANCH_SUPPORT_CAUTION_LABEL = (
    "SH-like local branch support is an approximate site-resampled likelihood-ranking aid over local rearrangements and must not be interpreted as bootstrap support or an exact SH-aLRT p-value"
)


@dataclass(frozen=True, slots=True)
class _BranchLocalArrangement:
    """One local topology arrangement scored against a single rooted branch."""

    branch_id: str
    node_label: str | None
    descendant_taxa: list[str]
    candidate_tree_id: str
    candidate_tree_label: str
    local_arrangement_kind: str
    tree: PhyloTree
    topology_fingerprint: str


def validate_nucleotide_sh_like_branch_support_replicate_count(
    resampling_replicate_count: int,
) -> int:
    """Require at least one site-resampling replicate for native SH-like support."""
    if resampling_replicate_count < 1:
        raise ValueError("resampling_replicate_count must be at least one")
    return resampling_replicate_count


def evaluate_nucleotide_sh_like_branch_support(
    tree: PhyloTree | Path,
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
) -> NucleotideShLikeBranchSupportReport:
    """Estimate native SH-like branch support from local rooted-NNI likelihood comparisons."""
    resolved_tree, tree_path = _resolve_sh_like_support_tree(tree)
    resolved_records, alignment_path = _resolve_sh_like_support_alignment(records)
    validate_rooted_nni_tree(resolved_tree)
    validated_replicate_count = validate_nucleotide_sh_like_branch_support_replicate_count(
        resampling_replicate_count
    )
    branch_candidates = _group_branch_local_nni_candidates(resolved_tree)
    if not branch_candidates:
        raise ValueError(
            "SH-like branch support requires at least one informative internal branch with local rooted NNI alternatives"
        )

    branch_support_rows: list[NucleotideShLikeBranchSupportRow] = []
    local_topology_rows: list[NucleotideShLikeBranchSupportLocalTopologyRow] = []
    resampling_rows: list[NucleotideShLikeBranchSupportResamplingRow] = []
    report_matrix = None
    for branch_id, candidates in branch_candidates.items():
        arrangements = _build_branch_local_arrangements(
            resolved_tree,
            branch_id=branch_id,
            candidates=candidates,
        )
        matrix_report = evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
            [arrangement.tree for arrangement in arrangements],
            resolved_records,
            model_name=model_name,
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        report_matrix = matrix_report
        candidate_vectors = resolve_candidate_tree_site_likelihood_vectors(matrix_report)
        observed_best = select_observed_best_candidate_tree(candidate_vectors)
        reference_vector = candidate_vectors[0]
        alternative_vectors = candidate_vectors[1:]
        best_alternative = select_observed_best_candidate_tree(alternative_vectors)
        local_topology_rows.extend(
            _build_sh_like_local_topology_rows(
                arrangements,
                candidate_vectors,
                observed_best=observed_best,
            )
        )
        support_count, branch_resampling_rows = _summarize_sh_like_branch_resampling(
            arrangements,
            candidate_vectors,
            validated_replicate_count=validated_replicate_count,
            resampling_seed=resampling_seed,
        )
        resampling_rows.extend(branch_resampling_rows)
        branch_support_rows.append(
            NucleotideShLikeBranchSupportRow(
                branch_id=arrangements[0].branch_id,
                node_label=arrangements[0].node_label,
                descendant_taxa=arrangements[0].descendant_taxa,
                alternative_arrangement_count=len(alternative_vectors),
                reference_log_likelihood=reference_vector.log_likelihood,
                best_alternative_tree_id=best_alternative.candidate_tree_id,
                best_alternative_tree_label=_candidate_tree_label_for_id(
                    arrangements,
                    best_alternative.candidate_tree_id,
                ),
                best_alternative_topology_fingerprint=_topology_fingerprint_for_id(
                    arrangements,
                    best_alternative.candidate_tree_id,
                ),
                best_alternative_log_likelihood=best_alternative.log_likelihood,
                observed_delta_log_likelihood=(
                    reference_vector.log_likelihood - best_alternative.log_likelihood
                ),
                reference_is_observed_best=(
                    reference_vector.candidate_tree_id == observed_best.candidate_tree_id
                ),
                support_replicate_count=support_count,
                support_fraction=support_count / float(validated_replicate_count),
                support_percent=(
                    100.0 * support_count / float(validated_replicate_count)
                ),
                caution_label=_SH_LIKE_BRANCH_SUPPORT_CAUTION_LABEL,
            )
        )
    if report_matrix is None:
        raise AssertionError("SH-like branch support evaluation did not score any branches")
    return NucleotideShLikeBranchSupportReport(
        algorithm="nucleotide-sh-like-branch-support",
        tree_path=None if tree_path is None else str(tree_path),
        alignment_path=None if alignment_path is None else str(alignment_path),
        model_name=report_matrix.model_name,
        taxa=report_matrix.taxa,
        branch_count=len(branch_support_rows),
        site_count=report_matrix.site_count,
        pattern_count=report_matrix.pattern_count,
        compression_used=report_matrix.compression_used,
        parameter_values=report_matrix.parameter_values,
        reference_tree_newick=resolved_tree.to_newick(),
        resampling_method=_SH_LIKE_BRANCH_SUPPORT_RESAMPLING_METHOD,
        resampling_replicate_count=validated_replicate_count,
        resampling_seed=resampling_seed,
        caution_label=_SH_LIKE_BRANCH_SUPPORT_CAUTION_LABEL,
        branch_support_rows=branch_support_rows,
        local_topology_rows=local_topology_rows,
        resampling_rows=resampling_rows,
    )


def evaluate_nucleotide_sh_like_branch_support_from_alignment(
    tree_path: Path,
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
) -> NucleotideShLikeBranchSupportReport:
    """Estimate native SH-like branch support from one tree path and FASTA alignment."""
    return evaluate_nucleotide_sh_like_branch_support(
        tree_path,
        alignment_path,
        model_name=model_name,
        resampling_replicate_count=resampling_replicate_count,
        resampling_seed=resampling_seed,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def write_nucleotide_sh_like_branch_support_table(
    path: Path,
    report: NucleotideShLikeBranchSupportReport,
) -> Path:
    """Write one native SH-like branch support summary table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "node_label",
        "descendant_taxa",
        "alternative_arrangement_count",
        "reference_log_likelihood",
        "best_alternative_tree_id",
        "best_alternative_tree_label",
        "best_alternative_topology_fingerprint",
        "best_alternative_log_likelihood",
        "observed_delta_log_likelihood",
        "reference_is_observed_best",
        "support_replicate_count",
        "support_fraction",
        "support_percent",
        "caution_label",
    ]
    rows = ["\t".join(columns)]
    for row in report.branch_support_rows:
        rows.append(
            "\t".join(
                [
                    row.branch_id,
                    "" if row.node_label is None else row.node_label,
                    "|".join(row.descendant_taxa),
                    str(row.alternative_arrangement_count),
                    repr(row.reference_log_likelihood),
                    row.best_alternative_tree_id,
                    row.best_alternative_tree_label,
                    row.best_alternative_topology_fingerprint,
                    repr(row.best_alternative_log_likelihood),
                    repr(row.observed_delta_log_likelihood),
                    str(row.reference_is_observed_best).lower(),
                    str(row.support_replicate_count),
                    repr(row.support_fraction),
                    repr(row.support_percent),
                    row.caution_label,
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_sh_like_branch_support_local_topology_table(
    path: Path,
    report: NucleotideShLikeBranchSupportReport,
) -> Path:
    """Write one local-topology likelihood ledger for SH-like branch support."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "node_label",
        "descendant_taxa",
        "candidate_tree_id",
        "candidate_tree_label",
        "local_arrangement_kind",
        "topology_fingerprint",
        "observed_log_likelihood",
        "observed_delta_log_likelihood",
        "observed_best_local_arrangement",
        "tree_newick",
    ]
    rows = ["\t".join(columns)]
    for row in report.local_topology_rows:
        rows.append(
            "\t".join(
                [
                    row.branch_id,
                    "" if row.node_label is None else row.node_label,
                    "|".join(row.descendant_taxa),
                    row.candidate_tree_id,
                    row.candidate_tree_label,
                    row.local_arrangement_kind,
                    row.topology_fingerprint,
                    repr(row.observed_log_likelihood),
                    repr(row.observed_delta_log_likelihood),
                    str(row.observed_best_local_arrangement).lower(),
                    row.tree_newick,
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_sh_like_branch_support_resampling_table(
    path: Path,
    report: NucleotideShLikeBranchSupportReport,
) -> Path:
    """Write one site-resampled local branch-support comparison ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "branch_id",
        "descendant_taxa",
        "replicate_index",
        "reference_resampled_log_likelihood",
        "best_local_tree_id",
        "best_local_tree_label",
        "best_local_resampled_log_likelihood",
        "best_alternative_tree_id",
        "best_alternative_tree_label",
        "best_alternative_resampled_log_likelihood",
        "reference_delta_log_likelihood",
        "reference_matches_or_beats_alternatives",
    ]
    rows = ["\t".join(columns)]
    for row in report.resampling_rows:
        rows.append(
            "\t".join(
                [
                    row.branch_id,
                    "|".join(row.descendant_taxa),
                    str(row.replicate_index),
                    repr(row.reference_resampled_log_likelihood),
                    row.best_local_tree_id,
                    row.best_local_tree_label,
                    repr(row.best_local_resampled_log_likelihood),
                    row.best_alternative_tree_id,
                    row.best_alternative_tree_label,
                    repr(row.best_alternative_resampled_log_likelihood),
                    repr(row.reference_delta_log_likelihood),
                    str(row.reference_matches_or_beats_alternatives).lower(),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_sh_like_branch_support_run_json(
    path: Path,
    report: NucleotideShLikeBranchSupportReport,
) -> Path:
    """Write one governed JSON payload for native SH-like branch support."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(report), indent=2) + "\n", encoding="utf-8")
    return path


def write_nucleotide_sh_like_branch_support_artifacts(
    output_dir: Path,
    report: NucleotideShLikeBranchSupportReport,
) -> dict[str, Path]:
    """Materialize one governed artifact family for native SH-like branch support."""
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_tree_path = output_dir / "reference_tree.nwk"
    write_newick(reference_tree_path, loads_newick(report.reference_tree_newick))
    outputs = {
        "reference_tree_path": reference_tree_path,
        "branch_support_path": write_nucleotide_sh_like_branch_support_table(
            output_dir / "branch_support.tsv",
            report,
        ),
        "local_topology_path": write_nucleotide_sh_like_branch_support_local_topology_table(
            output_dir / "local_topologies.tsv",
            report,
        ),
        "resampling_path": write_nucleotide_sh_like_branch_support_resampling_table(
            output_dir / "resampling.tsv",
            report,
        ),
        "run_json_path": write_nucleotide_sh_like_branch_support_run_json(
            output_dir / "run.json",
            report,
        ),
    }
    return outputs


def _resolve_sh_like_support_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_tree(tree) if isinstance(tree, Path) else tree.copy().refresh()
    return resolved_tree.refresh(), resolved_tree_path


def _resolve_sh_like_support_alignment(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    resolved_alignment_path = records if isinstance(records, Path) else None
    resolved_records = (
        load_fasta_alignment(records) if isinstance(records, Path) else list(records)
    )
    return resolved_records, resolved_alignment_path


def _group_branch_local_nni_candidates(
    tree: PhyloTree,
) -> dict[str, list[RootedNniMoveCandidate]]:
    grouped_candidates: dict[str, list[RootedNniMoveCandidate]] = {}
    for candidate in iter_rooted_nni_move_candidates(tree):
        grouped_candidates.setdefault(candidate.child_node_id, []).append(candidate)
    return dict(
        sorted(
            grouped_candidates.items(),
            key=lambda item: (
                len(tree.node_by_id(item[0]).descendant_taxa),
                tuple(tree.node_by_id(item[0]).descendant_taxa),
            ),
        )
    )


def _build_branch_local_arrangements(
    reference_tree: PhyloTree,
    *,
    branch_id: str,
    candidates: list[RootedNniMoveCandidate],
) -> list[_BranchLocalArrangement]:
    reference_node = reference_tree.node_by_id(branch_id)
    arrangements = [
        _BranchLocalArrangement(
            branch_id=branch_id,
            node_label=reference_node.name,
            descendant_taxa=list(reference_node.descendant_taxa),
            candidate_tree_id="candidate-tree-1",
            candidate_tree_label="reference-local-arrangement",
            local_arrangement_kind="reference",
            tree=reference_tree.copy().refresh(),
            topology_fingerprint=rooted_topology_fingerprint(reference_tree),
        )
    ]
    alternative_fingerprints: set[str] = set()
    for alternative_index, candidate in enumerate(candidates, start=1):
        moved_tree = apply_rooted_nni_move(reference_tree, candidate)
        topology_fingerprint = rooted_topology_fingerprint(moved_tree)
        if topology_fingerprint in alternative_fingerprints:
            continue
        alternative_fingerprints.add(topology_fingerprint)
        arrangements.append(
            _BranchLocalArrangement(
                branch_id=branch_id,
                node_label=reference_node.name,
                descendant_taxa=list(reference_node.descendant_taxa),
                candidate_tree_id=f"candidate-tree-{alternative_index + 1}",
                candidate_tree_label=f"alternative-local-arrangement-{alternative_index}",
                local_arrangement_kind="alternative",
                tree=moved_tree.refresh(),
                topology_fingerprint=topology_fingerprint,
            )
        )
    if len(arrangements) != 3:
        raise ValueError(
            "SH-like branch support requires exactly two distinct local rooted NNI alternatives per informative branch"
        )
    return arrangements


def _build_sh_like_local_topology_rows(
    arrangements: list[_BranchLocalArrangement],
    candidate_vectors: list[CandidateTreeSiteLikelihoodVector],
    *,
    observed_best: CandidateTreeSiteLikelihoodVector,
) -> list[NucleotideShLikeBranchSupportLocalTopologyRow]:
    rows: list[NucleotideShLikeBranchSupportLocalTopologyRow] = []
    for arrangement, candidate in zip(arrangements, candidate_vectors, strict=True):
        rows.append(
            NucleotideShLikeBranchSupportLocalTopologyRow(
                branch_id=arrangement.branch_id,
                node_label=arrangement.node_label,
                descendant_taxa=arrangement.descendant_taxa,
                candidate_tree_id=arrangement.candidate_tree_id,
                candidate_tree_label=arrangement.candidate_tree_label,
                local_arrangement_kind=arrangement.local_arrangement_kind,
                tree_newick=candidate.tree_newick,
                topology_fingerprint=arrangement.topology_fingerprint,
                observed_log_likelihood=candidate.log_likelihood,
                observed_delta_log_likelihood=(
                    observed_best.log_likelihood - candidate.log_likelihood
                ),
                observed_best_local_arrangement=(
                    candidate.candidate_tree_id == observed_best.candidate_tree_id
                ),
            )
        )
    return rows


def _summarize_sh_like_branch_resampling(
    arrangements: list[_BranchLocalArrangement],
    candidate_vectors: list[CandidateTreeSiteLikelihoodVector],
    *,
    validated_replicate_count: int,
    resampling_seed: int,
) -> tuple[int, list[NucleotideShLikeBranchSupportResamplingRow]]:
    rng = random.Random(resampling_seed)  # nosec B311
    reference_candidate = candidate_vectors[0]
    alternative_candidates = candidate_vectors[1:]
    site_count = len(reference_candidate.site_log_likelihoods)
    support_count = 0
    rows: list[NucleotideShLikeBranchSupportResamplingRow] = []
    for replicate_index in range(1, validated_replicate_count + 1):
        sampled_indices = [rng.randrange(site_count) for _ in range(site_count)]
        resampled_totals = {
            candidate.candidate_tree_id: sum(
                candidate.site_log_likelihoods[index] for index in sampled_indices
            )
            for candidate in candidate_vectors
        }
        best_local = _select_best_resampled_candidate(
            candidate_vectors,
            resampled_totals=resampled_totals,
        )
        best_alternative = _select_best_resampled_candidate(
            alternative_candidates,
            resampled_totals=resampled_totals,
        )
        reference_total = resampled_totals[reference_candidate.candidate_tree_id]
        best_alternative_total = resampled_totals[best_alternative.candidate_tree_id]
        reference_matches_or_beats_alternatives = (
            reference_total > best_alternative_total
            or math.isclose(
                reference_total,
                best_alternative_total,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        )
        if reference_matches_or_beats_alternatives:
            support_count += 1
        rows.append(
            NucleotideShLikeBranchSupportResamplingRow(
                branch_id=arrangements[0].branch_id,
                descendant_taxa=arrangements[0].descendant_taxa,
                replicate_index=replicate_index,
                reference_resampled_log_likelihood=reference_total,
                best_local_tree_id=best_local.candidate_tree_id,
                best_local_tree_label=_candidate_tree_label_for_id(
                    arrangements,
                    best_local.candidate_tree_id,
                ),
                best_local_resampled_log_likelihood=resampled_totals[
                    best_local.candidate_tree_id
                ],
                best_alternative_tree_id=best_alternative.candidate_tree_id,
                best_alternative_tree_label=_candidate_tree_label_for_id(
                    arrangements,
                    best_alternative.candidate_tree_id,
                ),
                best_alternative_resampled_log_likelihood=best_alternative_total,
                reference_delta_log_likelihood=reference_total - best_alternative_total,
                reference_matches_or_beats_alternatives=(
                    reference_matches_or_beats_alternatives
                ),
            )
        )
    return support_count, rows


def _select_best_resampled_candidate(
    candidates: list[CandidateTreeSiteLikelihoodVector],
    *,
    resampled_totals: dict[str, float],
) -> CandidateTreeSiteLikelihoodVector:
    if not candidates:
        raise ValueError("resampled candidate selection requires at least one candidate")
    best_candidate = candidates[0]
    for candidate in candidates[1:]:
        if _prefer_higher_resampled_candidate(
            candidate,
            best_candidate,
            resampled_totals=resampled_totals,
        ):
            best_candidate = candidate
    return best_candidate


def _prefer_higher_resampled_candidate(
    left: CandidateTreeSiteLikelihoodVector,
    right: CandidateTreeSiteLikelihoodVector,
    *,
    resampled_totals: dict[str, float],
) -> bool:
    left_total = resampled_totals[left.candidate_tree_id]
    right_total = resampled_totals[right.candidate_tree_id]
    if left_total > right_total and not math.isclose(
        left_total,
        right_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        return True
    if right_total > left_total and not math.isclose(
        left_total,
        right_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        return False
    return prefer_higher_likelihood_candidate(left, right)


def _candidate_tree_label_for_id(
    arrangements: list[_BranchLocalArrangement],
    candidate_tree_id: str,
) -> str:
    arrangement = next(
        arrangement
        for arrangement in arrangements
        if arrangement.candidate_tree_id == candidate_tree_id
    )
    return arrangement.candidate_tree_label


def _topology_fingerprint_for_id(
    arrangements: list[_BranchLocalArrangement],
    candidate_tree_id: str,
) -> str:
    arrangement = next(
        arrangement
        for arrangement in arrangements
        if arrangement.candidate_tree_id == candidate_tree_id
    )
    return arrangement.topology_fingerprint

