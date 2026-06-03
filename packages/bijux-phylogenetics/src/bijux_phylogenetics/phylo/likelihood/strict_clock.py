from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.branch_likelihood_diagnostics import (
    summarize_fixed_tree_branch_likelihood_diagnostics,
    write_branch_likelihood_diagnostic_table,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
    normalize_unambiguous_dna_records,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import (
    _evaluate_jc69_tree_likelihood_from_patterns,
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.site_log_likelihoods import (
    evaluate_selected_dna_site_log_likelihoods_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.sites import write_site_log_likelihood_table
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

from .models import StrictClockBranchRow, StrictClockLikelihoodReport
from .validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)

_SUPPORTED_MODELS = {"jc69": "JC69"}


def _resolve_strict_clock_model_name(model: str) -> str:
    normalized = model.strip().casefold()
    resolved = _SUPPORTED_MODELS.get(normalized)
    if resolved is None:
        raise ValueError(
            f"strict-clock likelihood currently supports only JC69, got {model}"
        )
    return resolved


def _validate_clock_rate_bounds(
    *,
    lower_clock_rate_bound: float,
    upper_clock_rate_bound: float,
) -> None:
    if lower_clock_rate_bound <= 0.0:
        raise InvalidBranchLengthError("strict-clock rate lower bound must be positive")
    if upper_clock_rate_bound <= lower_clock_rate_bound:
        raise InvalidBranchLengthError(
            "strict-clock rate bounds must be strictly increasing"
        )


def _copy_with_scaled_branch_lengths(
    time_tree: PhyloTree,
    *,
    clock_rate: float,
) -> PhyloTree:
    scaled_tree = time_tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree


def _build_branch_rows(
    *,
    time_tree: PhyloTree,
    scaled_tree: PhyloTree,
    optimized_clock_rate: float,
) -> list[StrictClockBranchRow]:
    time_by_node_id = {
        child.node_id: float(child.branch_length or 0.0)
        for _parent, child in time_tree.iter_edges()
    }
    branch_rows: list[StrictClockBranchRow] = []
    for _parent, child in scaled_tree.iter_edges():
        if child.node_id is None:
            raise ValueError("strict-clock branch rows require stable node ids")
        branch_rows.append(
            StrictClockBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                time_duration=time_by_node_id[child.node_id],
                optimized_branch_length=float(child.branch_length or 0.0),
                optimized_clock_rate=optimized_clock_rate,
            )
        )
    return branch_rows


def _evaluate_jc69_site_log_likelihood_rows(
    tree: PhyloTree,
    *,
    compressed_patterns,
):
    return evaluate_selected_dna_site_log_likelihoods_from_patterns(
        tree,
        compressed_patterns,
        model_name="JC69",
        state_count=len(UNIFORM_DNA_ROOT_PRIOR),
        observation_policy="reject",
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
        parameter_values={},
        transition_matrix_for_child=lambda child: jc69_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0)
        ),
    )


def _fit_jc69_strict_clock_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    lower_clock_rate_bound: float,
    upper_clock_rate_bound: float,
) -> StrictClockLikelihoodReport:
    _validate_clock_rate_bounds(
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
    )
    time_tree = tree.copy()
    validate_explicit_branch_lengths(time_tree, model_name="strict-clock JC69")
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name="strict-clock JC69",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    validate_tree_taxa_against_patterns(
        time_tree,
        compressed_patterns,
        model_name="strict-clock JC69",
    )

    initial_clock_rate = min(
        max(1.0, lower_clock_rate_bound),
        upper_clock_rate_bound,
    )
    initial_scaled_tree = _copy_with_scaled_branch_lengths(
        time_tree,
        clock_rate=initial_clock_rate,
    )
    initial_report = _evaluate_jc69_tree_likelihood_from_patterns(
        initial_scaled_tree,
        compressed_patterns,
        observation_policy="reject",
        gap_state_frequency=None,
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
    )

    def evaluate_candidate(
        clock_rate: float,
    ) -> tuple[PhyloTree, float]:
        scaled_tree = _copy_with_scaled_branch_lengths(
            time_tree,
            clock_rate=clock_rate,
        )
        report = _evaluate_jc69_tree_likelihood_from_patterns(
            scaled_tree,
            compressed_patterns,
            observation_policy="reject",
            gap_state_frequency=None,
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
        )
        return scaled_tree, report.log_likelihood

    search_result = run_bounded_likelihood_search(
        lower_bound=lower_clock_rate_bound,
        upper_bound=upper_clock_rate_bound,
        evaluate=evaluate_candidate,
    )
    optimized_scaled_tree = search_result.payload
    optimized_report = _evaluate_jc69_tree_likelihood_from_patterns(
        optimized_scaled_tree,
        compressed_patterns,
        observation_policy="reject",
        gap_state_frequency=None,
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
    )
    site_log_likelihood_report = _evaluate_jc69_site_log_likelihood_rows(
        optimized_scaled_tree,
        compressed_patterns=compressed_patterns,
    )
    branch_diagnostic_report = summarize_fixed_tree_branch_likelihood_diagnostics(
        optimized_scaled_tree,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        model_name="JC69",
        baseline_log_likelihood=optimized_report.log_likelihood,
        evaluate_tree_log_likelihood=lambda candidate_tree: (
            _evaluate_jc69_tree_likelihood_from_patterns(
                candidate_tree,
                compressed_patterns,
                observation_policy="reject",
                gap_state_frequency=None,
                root_prior=UNIFORM_DNA_ROOT_PRIOR,
            ).log_likelihood
        ),
    )
    optimized_clock_rate = float(search_result.parameter_value)
    branch_rows = _build_branch_rows(
        time_tree=time_tree,
        scaled_tree=optimized_scaled_tree,
        optimized_clock_rate=optimized_clock_rate,
    )
    return StrictClockLikelihoodReport(
        model_name="JC69",
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        branch_count=len(branch_rows),
        compression_used=True,
        time_tree_newick=dumps_newick(time_tree),
        scaled_tree_newick=dumps_newick(optimized_scaled_tree),
        initial_clock_rate=float(initial_clock_rate),
        optimized_clock_rate=optimized_clock_rate,
        initial_log_likelihood=initial_report.log_likelihood,
        optimized_log_likelihood=optimized_report.log_likelihood,
        parameter_count=1,
        aic=(-2.0 * optimized_report.log_likelihood) + (2.0 * 1.0),
        function_evaluation_count=search_result.function_evaluation_count + 1,
        converged=search_result.converged,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
        branch_rows=branch_rows,
        site_log_likelihoods=site_log_likelihood_report.site_log_likelihoods,
        branch_likelihood_diagnostics=branch_diagnostic_report.branch_diagnostics,
    )


def fit_strict_clock_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
) -> StrictClockLikelihoodReport:
    """Fit one global strict-clock rate on one time-scaled tree and one alignment."""
    resolved_model = _resolve_strict_clock_model_name(model)
    if (
        resolved_model != "JC69"
    ):  # pragma: no cover - guard for future dispatch expansion
        raise ValueError(f"unsupported strict-clock likelihood model: {resolved_model}")
    return _fit_jc69_strict_clock_likelihood(
        tree,
        records,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
    )


def fit_strict_clock_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
) -> StrictClockLikelihoodReport:
    """Fit one strict-clock likelihood report from one time-tree path and one alignment."""
    return fit_strict_clock_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model=model,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
    )


def write_strict_clock_branch_table(
    path: Path,
    report: StrictClockLikelihoodReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "branch_id",
                "child_name",
                "descendant_taxa",
                "time_duration",
                "optimized_branch_length",
                "optimized_clock_rate",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.branch_rows:
            writer.writerow(
                {
                    "branch_id": row.branch_id,
                    "child_name": row.child_name or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "time_duration": format(row.time_duration, ".15g"),
                    "optimized_branch_length": format(
                        row.optimized_branch_length,
                        ".15g",
                    ),
                    "optimized_clock_rate": format(
                        row.optimized_clock_rate,
                        ".15g",
                    ),
                }
            )
    return path


def write_strict_clock_run_json(
    path: Path,
    report: StrictClockLikelihoodReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload.pop("branch_likelihood_diagnostics", None)
    payload.pop("site_log_likelihoods", None)
    path.write_text(
        json.dumps(payload, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_strict_clock_likelihood_artifacts(
    out_dir: Path,
    report: StrictClockLikelihoodReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    scaled_tree_path = write_newick(
        out_dir / "scaled_tree.nwk",
        loads_newick(report.scaled_tree_newick),
    )
    branch_table_path = write_strict_clock_branch_table(
        out_dir / "branch_rates.tsv",
        report,
    )
    site_log_likelihood_path = write_site_log_likelihood_table(
        out_dir / "site_log_likelihoods.tsv",
        report,
    )
    branch_likelihood_diagnostic_path = write_branch_likelihood_diagnostic_table(
        out_dir / "branch_likelihood_diagnostics.tsv",
        report,
    )
    run_json_path = write_strict_clock_run_json(out_dir / "run.json", report)
    return {
        "scaled_tree_path": scaled_tree_path,
        "branch_table_path": branch_table_path,
        "site_log_likelihood_path": site_log_likelihood_path,
        "branch_likelihood_diagnostic_path": branch_likelihood_diagnostic_path,
        "run_json_path": run_json_path,
    }
