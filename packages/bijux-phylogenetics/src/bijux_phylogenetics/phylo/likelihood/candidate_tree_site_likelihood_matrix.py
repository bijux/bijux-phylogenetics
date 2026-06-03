from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import (
    load_newick_tree_set,
    loads_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import normalize_unambiguous_dna_records
from bijux_phylogenetics.phylo.likelihood.models import (
    CandidateTreeSiteLikelihoodMatrixReport,
    CandidateTreeSiteLikelihoodRow,
    CandidateTreeSiteLikelihoodSummary,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.site_log_likelihoods import (
    evaluate_selected_dna_site_log_likelihoods_from_patterns,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError, TreeParseError

_CANDIDATE_TREE_COMPARISON_CAUTION_LABEL = "all candidate trees are rescored under one shared alignment/model surface; prior per-tree fitted-model differences are not preserved in this comparison"


@dataclass(frozen=True, slots=True)
class _CandidateTreeRecord:
    """One resolved candidate tree with a stable identifier inside matrix evaluation."""

    candidate_tree_id: str
    candidate_tree_label: str
    tree: PhyloTree


def evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
    trees: list[PhyloTree] | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> CandidateTreeSiteLikelihoodMatrixReport:
    """Evaluate one shared nucleotide likelihood model across multiple candidate trees."""
    resolved_tree_records, resolved_tree_set_path = resolve_candidate_tree_records(
        trees
    )
    resolved_records, resolved_alignment_path = (
        resolve_candidate_tree_alignment_records(records)
    )
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        resolved_records,
        model_name=normalized_model_name.upper(),
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} candidate tree site likelihood matrix",
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )

    candidate_trees: list[CandidateTreeSiteLikelihoodSummary] = []
    matrix_rows: list[CandidateTreeSiteLikelihoodRow] = []
    for candidate in resolved_tree_records:
        try:
            tree_report = evaluate_selected_dna_site_log_likelihoods_from_patterns(
                candidate.tree,
                compressed_patterns,
                model_name=specification.model_name,
                state_count=specification.state_count,
                observation_policy=specification.observation_policy,
                root_prior=specification.root_prior,
                parameter_values=specification.parameter_values,
                transition_matrix_for_child=lambda child: (
                    specification.transition_matrix_for_branch_length(
                        max(float(child.branch_length or 0.0), 0.0)
                    )
                ),
            )
        except PhylogeneticsError as error:
            raise ValueError(
                f"{candidate.candidate_tree_label} is incompatible with the shared alignment/model surface: {error}"
            ) from error
        candidate_trees.append(
            CandidateTreeSiteLikelihoodSummary(
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                topology_fingerprint=rooted_topology_fingerprint(candidate.tree),
                tree_newick=tree_report.tree_newick,
                log_likelihood=tree_report.log_likelihood,
                observed_delta_log_likelihood=0.0,
            )
        )
        matrix_rows.extend(
            CandidateTreeSiteLikelihoodRow(
                candidate_tree_id=candidate.candidate_tree_id,
                candidate_tree_label=candidate.candidate_tree_label,
                tree_newick=tree_report.tree_newick,
                pattern_id=row.pattern_id,
                pattern_weight=row.pattern_weight,
                site_position=row.site_position,
                site_states=row.site_states,
                log_likelihood=row.log_likelihood,
            )
            for row in tree_report.site_log_likelihoods
        )
    observed_best_candidate = select_observed_best_candidate_tree_summary(
        candidate_trees
    )
    candidate_trees = [
        CandidateTreeSiteLikelihoodSummary(
            candidate_tree_id=row.candidate_tree_id,
            candidate_tree_label=row.candidate_tree_label,
            topology_fingerprint=row.topology_fingerprint,
            tree_newick=row.tree_newick,
            log_likelihood=row.log_likelihood,
            observed_delta_log_likelihood=(
                observed_best_candidate.log_likelihood - row.log_likelihood
            ),
        )
        for row in candidate_trees
    ]
    return CandidateTreeSiteLikelihoodMatrixReport(
        model_name=specification.model_name,
        tree_set_path=None
        if resolved_tree_set_path is None
        else str(resolved_tree_set_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxa=compressed_patterns.taxon_order,
        tree_count=len(candidate_trees),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        expansion_policy="candidate-tree-expanded-site-rows",
        comparison_caution_label=_CANDIDATE_TREE_COMPARISON_CAUTION_LABEL,
        parameter_values=specification.parameter_values,
        candidate_trees=candidate_trees,
        matrix_rows=matrix_rows,
    )


def evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
    tree_set_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> CandidateTreeSiteLikelihoodMatrixReport:
    """Evaluate a candidate-tree site-likelihood matrix from tree-set and alignment paths."""
    return evaluate_nucleotide_candidate_tree_site_likelihood_matrix(
        tree_set_path,
        alignment_path,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def resolve_candidate_tree_records(
    trees: list[PhyloTree] | Path,
) -> tuple[list[_CandidateTreeRecord], Path | None]:
    """Resolve one candidate-tree set and assign stable candidate identifiers."""
    resolved_tree_set_path = trees if isinstance(trees, Path) else None
    resolved_trees = (
        load_newick_tree_set(trees) if isinstance(trees, Path) else list(trees)
    )
    if len(resolved_trees) < 2:
        raise TreeParseError(
            "candidate tree site likelihood matrix requires at least two candidate trees"
        )
    return [
        _CandidateTreeRecord(
            candidate_tree_id=f"candidate-tree-{index}",
            candidate_tree_label=f"candidate-tree-{index}",
            tree=tree.copy().refresh(),
        )
        for index, tree in enumerate(resolved_trees, start=1)
    ], resolved_tree_set_path


def resolve_candidate_tree_alignment_records(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    """Resolve one candidate-tree alignment input and preserve its path when present."""
    resolved_alignment_path = records if isinstance(records, Path) else None
    resolved_records = (
        load_fasta_alignment(records) if isinstance(records, Path) else list(records)
    )
    return resolved_records, resolved_alignment_path


def write_candidate_tree_likelihood_summary_table(
    path: Path,
    report: CandidateTreeSiteLikelihoodMatrixReport,
) -> Path:
    """Write one candidate-tree summary table with shared model totals."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "candidate_tree_id",
        "candidate_tree_label",
        "topology_fingerprint",
        "tree_newick",
        "log_likelihood",
        "observed_delta_log_likelihood",
    ]
    rows = ["\t".join(columns)]
    for row in report.candidate_trees:
        rows.append(
            "\t".join(
                [
                    row.candidate_tree_id,
                    row.candidate_tree_label,
                    row.topology_fingerprint,
                    row.tree_newick,
                    repr(row.log_likelihood),
                    repr(row.observed_delta_log_likelihood),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_candidate_tree_site_likelihood_matrix_table(
    path: Path,
    report: CandidateTreeSiteLikelihoodMatrixReport,
) -> Path:
    """Write one expanded candidate-tree by site likelihood matrix TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "model_name",
        "candidate_tree_id",
        "candidate_tree_label",
        "taxon_order",
        "pattern_id",
        "pattern_weight",
        "site_position",
        "site_states",
        "log_likelihood",
    ]
    rows = ["\t".join(columns)]
    taxon_order = "|".join(report.taxa)
    for row in report.matrix_rows:
        rows.append(
            "\t".join(
                [
                    report.model_name,
                    row.candidate_tree_id,
                    row.candidate_tree_label,
                    taxon_order,
                    row.pattern_id,
                    str(row.pattern_weight),
                    str(row.site_position),
                    "|".join(row.site_states),
                    repr(row.log_likelihood),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_candidate_tree_site_likelihood_matrix_run_json(
    path: Path,
    report: CandidateTreeSiteLikelihoodMatrixReport,
) -> Path:
    """Write one machine-readable candidate-tree site-likelihood matrix payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "model_name": report.model_name,
        "tree_set_path": report.tree_set_path,
        "alignment_path": report.alignment_path,
        "taxa": report.taxa,
        "tree_count": report.tree_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "compression_used": report.compression_used,
        "expansion_policy": report.expansion_policy,
        "comparison_caution_label": report.comparison_caution_label,
        "parameter_values": report.parameter_values,
        "candidate_trees": [
            {
                "candidate_tree_id": row.candidate_tree_id,
                "candidate_tree_label": row.candidate_tree_label,
                "topology_fingerprint": row.topology_fingerprint,
                "tree_newick": row.tree_newick,
                "log_likelihood": row.log_likelihood,
                "observed_delta_log_likelihood": row.observed_delta_log_likelihood,
            }
            for row in report.candidate_trees
        ],
        "matrix_rows": [
            {
                "candidate_tree_id": row.candidate_tree_id,
                "candidate_tree_label": row.candidate_tree_label,
                "tree_newick": row.tree_newick,
                "pattern_id": row.pattern_id,
                "pattern_weight": row.pattern_weight,
                "site_position": row.site_position,
                "site_states": list(row.site_states),
                "log_likelihood": row.log_likelihood,
            }
            for row in report.matrix_rows
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def select_observed_best_candidate_tree_summary(
    candidates: list[CandidateTreeSiteLikelihoodSummary],
) -> CandidateTreeSiteLikelihoodSummary:
    """Choose the observed best candidate summary with deterministic tie breaks."""
    if not candidates:
        raise ValueError(
            "candidate tree site likelihood matrix requires at least one candidate"
        )
    best_candidate = candidates[0]
    for candidate in candidates[1:]:
        if prefer_higher_likelihood_candidate_tree_summary(candidate, best_candidate):
            best_candidate = candidate
    return best_candidate


def prefer_higher_likelihood_candidate_tree_summary(
    left: CandidateTreeSiteLikelihoodSummary,
    right: CandidateTreeSiteLikelihoodSummary,
) -> bool:
    """Prefer higher likelihoods, then topology identity, across candidate trees."""
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
    if left.topology_fingerprint != right.topology_fingerprint:
        return left.topology_fingerprint < right.topology_fingerprint
    if left.tree_newick != right.tree_newick:
        return left.tree_newick < right.tree_newick
    return left.candidate_tree_id < right.candidate_tree_id


def write_candidate_tree_site_likelihood_matrix_artifacts(
    out_dir: Path,
    report: CandidateTreeSiteLikelihoodMatrixReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one candidate-tree site-likelihood matrix."""
    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_tree_path = write_newick_tree_set(
        out_dir / "candidate_trees.nwk",
        [loads_newick(row.tree_newick) for row in report.candidate_trees],
    )
    summary_path = write_candidate_tree_likelihood_summary_table(
        out_dir / "candidate_tree_summary.tsv",
        report,
    )
    matrix_path = write_candidate_tree_site_likelihood_matrix_table(
        out_dir / "site_likelihood_matrix.tsv",
        report,
    )
    run_json_path = write_candidate_tree_site_likelihood_matrix_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "candidate_tree_path": candidate_tree_path,
        "summary_path": summary_path,
        "matrix_path": matrix_path,
        "run_json_path": run_json_path,
    }
