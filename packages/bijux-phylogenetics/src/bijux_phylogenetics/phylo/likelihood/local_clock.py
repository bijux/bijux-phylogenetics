from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
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
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)

from .models import (
    LocalClockBranchRow,
    LocalClockLikelihoodReport,
    LocalClockRegimeRow,
)
from .sites import write_site_log_likelihood_table
from .strict_clock import (
    _evaluate_jc69_site_log_likelihood_rows,
    fit_strict_clock_likelihood,
)
from .validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)

_SUPPORTED_MODELS = {"jc69": "JC69"}
_BACKGROUND_REGIME_ID = "background"


@dataclass(frozen=True, slots=True)
class _LocalClockSelector:
    regime_id: str
    target_kind: str
    target_label: str | None
    descendant_taxa: tuple[str, ...]
    node_id: str
    matching_branch_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _BranchAssignment:
    regime_id: str
    target_kind: str


def _resolve_local_clock_model_name(model: str) -> str:
    normalized = model.strip().casefold()
    resolved = _SUPPORTED_MODELS.get(normalized)
    if resolved is None:
        raise ValueError(
            f"local-clock likelihood currently supports only JC69, got {model}"
        )
    return resolved


def _validate_clock_rate_bounds(
    *,
    lower_clock_rate_bound: float,
    upper_clock_rate_bound: float,
) -> None:
    if lower_clock_rate_bound <= 0.0:
        raise InvalidBranchLengthError("local-clock rate lower bound must be positive")
    if upper_clock_rate_bound <= lower_clock_rate_bound:
        raise InvalidBranchLengthError(
            "local-clock rate bounds must be strictly increasing"
        )


def _validate_max_coordinate_passes(max_coordinate_passes: int) -> None:
    if max_coordinate_passes < 1:
        raise PhylogeneticsError(
            "local-clock likelihood requires at least one coordinate pass",
            code="local_clock_error",
        )


def _load_local_clock_selectors(
    tree: PhyloTree,
    regime_path: Path,
) -> list[_LocalClockSelector]:
    with regime_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required_columns = {"regime_id", "target_kind", "target_label", "taxa"}
        missing_columns = sorted(required_columns - set(fieldnames))
        if missing_columns:
            raise PhylogeneticsError(
                "local clock regime table is missing required columns: "
                + ", ".join(missing_columns),
                code="local_clock_error",
            )
        rows = list(reader)

    if not rows:
        raise PhylogeneticsError(
            "local clock regime table must contain at least one regime row",
            code="local_clock_error",
        )

    node_by_descendant_taxa = {
        frozenset(node.descendant_taxa): node
        for node in tree.iter_nodes(order="preorder")
    }
    selectors: list[_LocalClockSelector] = []
    seen_regime_ids: set[str] = set()
    for row in rows:
        regime_id = (row.get("regime_id") or "").strip()
        if not regime_id:
            raise PhylogeneticsError(
                "local clock regime rows require a non-empty regime_id",
                code="local_clock_error",
            )
        if regime_id.casefold() == _BACKGROUND_REGIME_ID:
            raise PhylogeneticsError(
                "local clock regime_id 'background' is reserved for unassigned branches",
                code="local_clock_error",
            )
        if regime_id in seen_regime_ids:
            raise PhylogeneticsError(
                f"local clock regime_id '{regime_id}' is duplicated",
                code="local_clock_error",
            )
        seen_regime_ids.add(regime_id)

        target_kind = (row.get("target_kind") or "").strip().casefold()
        if target_kind not in {"branch", "clade"}:
            raise PhylogeneticsError(
                "local clock target_kind must be either 'branch' or 'clade'",
                code="local_clock_error",
            )
        target_label = (row.get("target_label") or "").strip() or None
        descendant_taxa = _parse_selector_taxa(
            row.get("taxa"),
            regime_id=regime_id,
        )
        node = node_by_descendant_taxa.get(frozenset(descendant_taxa))
        if node is None or node.node_id is None:
            raise PhylogeneticsError(
                f"local clock regime '{regime_id}' does not match one stable tree clade",
                code="local_clock_error",
            )
        if target_kind == "branch":
            if node is tree.root:
                raise PhylogeneticsError(
                    f"local clock branch regime '{regime_id}' cannot target the root",
                    code="local_clock_error",
                )
            matching_branch_ids = (node.node_id,)
        else:
            if node is tree.root:
                raise PhylogeneticsError(
                    f"local clock clade regime '{regime_id}' cannot target the whole tree; use the background rate instead",
                    code="local_clock_error",
                )
            target_taxa = frozenset(descendant_taxa)
            matching_branch_ids = tuple(
                child.node_id
                for _parent, child in tree.iter_edges()
                if child.node_id is not None
                and frozenset(child.descendant_taxa).issubset(target_taxa)
            )
        selectors.append(
            _LocalClockSelector(
                regime_id=regime_id,
                target_kind=target_kind,
                target_label=target_label,
                descendant_taxa=tuple(descendant_taxa),
                node_id=node.node_id,
                matching_branch_ids=matching_branch_ids,
            )
        )
    return selectors


def _parse_selector_taxa(
    raw_taxa: str | None,
    *,
    regime_id: str,
) -> list[str]:
    taxa = [token.strip() for token in (raw_taxa or "").split("|") if token.strip()]
    if not taxa:
        raise PhylogeneticsError(
            f"local clock regime '{regime_id}' requires at least one taxon",
            code="local_clock_error",
        )
    if len(set(taxa)) != len(taxa):
        raise PhylogeneticsError(
            f"local clock regime '{regime_id}' contains duplicate taxa",
            code="local_clock_error",
        )
    return sorted(taxa)


def _resolve_branch_assignments(
    tree: PhyloTree,
    selectors: list[_LocalClockSelector],
) -> tuple[dict[str, _BranchAssignment], dict[str, int]]:
    assignment_by_branch_id: dict[str, _BranchAssignment] = {}
    branch_count_by_regime_id: dict[str, int] = {_BACKGROUND_REGIME_ID: 0}
    for selector in selectors:
        branch_count_by_regime_id[selector.regime_id] = 0

    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                "local clock branch assignments require stable node ids",
                code="local_clock_error",
            )
        candidates = [
            selector
            for selector in selectors
            if child.node_id in selector.matching_branch_ids
        ]
        if not candidates:
            assignment = _BranchAssignment(
                regime_id=_BACKGROUND_REGIME_ID,
                target_kind="background",
            )
        else:
            ranked_candidates = sorted(
                candidates,
                key=lambda selector: (
                    0 if selector.target_kind == "branch" else 1,
                    len(selector.descendant_taxa),
                    selector.regime_id,
                ),
            )
            top_selector = ranked_candidates[0]
            top_priority = (
                0 if top_selector.target_kind == "branch" else 1,
                len(top_selector.descendant_taxa),
            )
            tied_regimes = [
                selector.regime_id
                for selector in ranked_candidates
                if (
                    (0 if selector.target_kind == "branch" else 1),
                    len(selector.descendant_taxa),
                )
                == top_priority
            ]
            if len(tied_regimes) > 1:
                raise PhylogeneticsError(
                    "local clock branch assignment is ambiguous for clade "
                    f"'{'|'.join(child.descendant_taxa)}': " + ", ".join(tied_regimes),
                    code="local_clock_error",
                )
            assignment = _BranchAssignment(
                regime_id=top_selector.regime_id,
                target_kind=top_selector.target_kind,
            )
        assignment_by_branch_id[child.node_id] = assignment
        branch_count_by_regime_id[assignment.regime_id] = (
            branch_count_by_regime_id.get(assignment.regime_id, 0) + 1
        )

    unused_selectors = [
        selector.regime_id
        for selector in selectors
        if branch_count_by_regime_id[selector.regime_id] == 0
    ]
    if unused_selectors:
        raise PhylogeneticsError(
            "local clock selector resolution left one or more regimes without branches: "
            + ", ".join(unused_selectors),
            code="local_clock_error",
        )
    return assignment_by_branch_id, branch_count_by_regime_id


def _copy_with_local_clock_branch_lengths(
    time_tree: PhyloTree,
    *,
    assignment_by_branch_id: dict[str, _BranchAssignment],
    rate_by_regime_id: dict[str, float],
) -> PhyloTree:
    scaled_tree = time_tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                "local clock branch scaling requires stable node ids",
                code="local_clock_error",
            )
        assignment = assignment_by_branch_id[child.node_id]
        clock_rate = rate_by_regime_id[assignment.regime_id]
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree


def _build_local_clock_regime_rows(
    *,
    selectors: list[_LocalClockSelector],
    branch_count_by_regime_id: dict[str, int],
    optimized_rate_by_regime_id: dict[str, float],
) -> list[LocalClockRegimeRow]:
    rows: list[LocalClockRegimeRow] = []
    if _BACKGROUND_REGIME_ID in optimized_rate_by_regime_id:
        rows.append(
            LocalClockRegimeRow(
                regime_id=_BACKGROUND_REGIME_ID,
                target_kind="background",
                target_label=None,
                descendant_taxa=[],
                node_id=None,
                branch_count=branch_count_by_regime_id[_BACKGROUND_REGIME_ID],
                optimized_clock_rate=optimized_rate_by_regime_id[_BACKGROUND_REGIME_ID],
            )
        )
    for selector in selectors:
        rows.append(
            LocalClockRegimeRow(
                regime_id=selector.regime_id,
                target_kind=selector.target_kind,
                target_label=selector.target_label,
                descendant_taxa=list(selector.descendant_taxa),
                node_id=selector.node_id,
                branch_count=branch_count_by_regime_id[selector.regime_id],
                optimized_clock_rate=optimized_rate_by_regime_id[selector.regime_id],
            )
        )
    return rows


def _build_local_clock_branch_rows(
    *,
    time_tree: PhyloTree,
    scaled_tree: PhyloTree,
    assignment_by_branch_id: dict[str, _BranchAssignment],
    optimized_rate_by_regime_id: dict[str, float],
) -> list[LocalClockBranchRow]:
    branch_rows: list[LocalClockBranchRow] = []
    for _parent, child in scaled_tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                "local clock branch rows require stable node ids",
                code="local_clock_error",
            )
        assignment = assignment_by_branch_id[child.node_id]
        branch_rows.append(
            LocalClockBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                regime_id=assignment.regime_id,
                target_kind=assignment.target_kind,
                time_duration=float(
                    time_tree.node_by_id(child.node_id).branch_length or 0.0
                ),
                optimized_branch_length=float(child.branch_length or 0.0),
                optimized_clock_rate=optimized_rate_by_regime_id[assignment.regime_id],
            )
        )
    return branch_rows


def fit_local_clock_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    regime_path: Path,
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
    max_coordinate_passes: int = 12,
) -> LocalClockLikelihoodReport:
    """Fit one local-clock JC69 likelihood model with user-defined branch regimes."""
    resolved_model = _resolve_local_clock_model_name(model)
    if resolved_model != "JC69":  # pragma: no cover
        raise ValueError(f"unsupported local-clock likelihood model: {resolved_model}")
    _validate_clock_rate_bounds(
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
    )
    _validate_max_coordinate_passes(max_coordinate_passes)

    time_tree = tree.copy()
    validate_explicit_branch_lengths(time_tree, model_name="local-clock JC69")
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name="local-clock JC69",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    validate_tree_taxa_against_patterns(
        time_tree,
        compressed_patterns,
        model_name="local-clock JC69",
    )

    strict_clock_report = fit_strict_clock_likelihood(
        time_tree.copy(),
        normalized_records,
        model=model,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
    )
    selectors = _load_local_clock_selectors(time_tree, regime_path)
    assignment_by_branch_id, branch_count_by_regime_id = _resolve_branch_assignments(
        time_tree,
        selectors,
    )

    initial_clock_rate = strict_clock_report.optimized_clock_rate
    parameter_names = []
    if branch_count_by_regime_id[_BACKGROUND_REGIME_ID] > 0:
        parameter_names.append(_BACKGROUND_REGIME_ID)
    parameter_names.extend(selector.regime_id for selector in selectors)
    initial_values = dict.fromkeys(parameter_names, initial_clock_rate)
    bounds_by_name = dict.fromkeys(
        parameter_names,
        (lower_clock_rate_bound, upper_clock_rate_bound),
    )

    def evaluate_candidate(
        rate_by_regime_id: dict[str, float],
    ) -> tuple[PhyloTree, float]:
        scaled_tree = _copy_with_local_clock_branch_lengths(
            time_tree,
            assignment_by_branch_id=assignment_by_branch_id,
            rate_by_regime_id=rate_by_regime_id,
        )
        report = _evaluate_jc69_tree_likelihood_from_patterns(
            scaled_tree,
            compressed_patterns,
            observation_policy="reject",
            gap_state_frequency=None,
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
        )
        return scaled_tree, report.log_likelihood

    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate_candidate,
        max_coordinate_passes=max_coordinate_passes,
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
        model_name="JC69 local-clock",
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
    optimized_rate_by_regime_id = dict(search_result.parameter_values)
    regime_rows = _build_local_clock_regime_rows(
        selectors=selectors,
        branch_count_by_regime_id=branch_count_by_regime_id,
        optimized_rate_by_regime_id=optimized_rate_by_regime_id,
    )
    branch_rows = _build_local_clock_branch_rows(
        time_tree=time_tree,
        scaled_tree=optimized_scaled_tree,
        assignment_by_branch_id=assignment_by_branch_id,
        optimized_rate_by_regime_id=optimized_rate_by_regime_id,
    )
    parameter_count = len(optimized_rate_by_regime_id)
    aic = (-2.0 * optimized_report.log_likelihood) + (2.0 * float(parameter_count))
    aic_delta_vs_strict_clock = aic - strict_clock_report.aic
    preferred_model_by_aic = (
        "local-clock" if aic < (strict_clock_report.aic - 1e-12) else "strict-clock"
    )
    return LocalClockLikelihoodReport(
        model_name="JC69 local-clock",
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        branch_count=len(branch_rows),
        regime_count=len(regime_rows),
        compression_used=True,
        time_tree_newick=dumps_newick(time_tree),
        scaled_tree_newick=dumps_newick(optimized_scaled_tree),
        strict_clock_rate=strict_clock_report.optimized_clock_rate,
        strict_clock_log_likelihood=strict_clock_report.optimized_log_likelihood,
        strict_clock_aic=strict_clock_report.aic,
        initial_clock_rate=initial_clock_rate,
        optimized_log_likelihood=optimized_report.log_likelihood,
        parameter_count=parameter_count,
        aic=aic,
        aic_delta_vs_strict_clock=aic_delta_vs_strict_clock,
        preferred_model_by_aic=preferred_model_by_aic,
        function_evaluation_count=search_result.function_evaluation_count + 1,
        optimization_pass_count=search_result.optimization_pass_count,
        converged=search_result.converged,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
        branch_rows=branch_rows,
        regime_rows=regime_rows,
        site_log_likelihoods=site_log_likelihood_report.site_log_likelihoods,
        branch_likelihood_diagnostics=branch_diagnostic_report.branch_diagnostics,
    )


def fit_local_clock_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    regime_path: Path,
    *,
    model: str = "jc69",
    lower_clock_rate_bound: float = 1e-6,
    upper_clock_rate_bound: float = 5.0,
    max_coordinate_passes: int = 12,
) -> LocalClockLikelihoodReport:
    """Fit one local-clock likelihood report from one time tree, alignment, and regime table."""
    return fit_local_clock_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        regime_path,
        model=model,
        lower_clock_rate_bound=lower_clock_rate_bound,
        upper_clock_rate_bound=upper_clock_rate_bound,
        max_coordinate_passes=max_coordinate_passes,
    )


def write_local_clock_branch_table(
    path: Path,
    report: LocalClockLikelihoodReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "branch_id",
                "child_name",
                "descendant_taxa",
                "regime_id",
                "target_kind",
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
                    "regime_id": row.regime_id,
                    "target_kind": row.target_kind,
                    "time_duration": format(row.time_duration, ".15g"),
                    "optimized_branch_length": format(
                        row.optimized_branch_length,
                        ".15g",
                    ),
                    "optimized_clock_rate": format(row.optimized_clock_rate, ".15g"),
                }
            )
    return path


def write_local_clock_regime_table(
    path: Path,
    report: LocalClockLikelihoodReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "regime_id",
                "target_kind",
                "target_label",
                "descendant_taxa",
                "node_id",
                "branch_count",
                "optimized_clock_rate",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.regime_rows:
            writer.writerow(
                {
                    "regime_id": row.regime_id,
                    "target_kind": row.target_kind,
                    "target_label": row.target_label or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "node_id": row.node_id or "",
                    "branch_count": row.branch_count,
                    "optimized_clock_rate": format(row.optimized_clock_rate, ".15g"),
                }
            )
    return path


def write_local_clock_run_json(
    path: Path,
    report: LocalClockLikelihoodReport,
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


def write_local_clock_likelihood_artifacts(
    out_dir: Path,
    report: LocalClockLikelihoodReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    scaled_tree_path = write_newick(
        out_dir / "scaled_tree.nwk",
        loads_newick(report.scaled_tree_newick),
    )
    branch_table_path = write_local_clock_branch_table(
        out_dir / "branch_rates.tsv",
        report,
    )
    regime_table_path = write_local_clock_regime_table(
        out_dir / "regimes.tsv",
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
    run_json_path = write_local_clock_run_json(out_dir / "run.json", report)
    return {
        "scaled_tree_path": scaled_tree_path,
        "branch_table_path": branch_table_path,
        "regime_table_path": regime_table_path,
        "site_log_likelihood_path": site_log_likelihood_path,
        "branch_likelihood_diagnostic_path": branch_likelihood_diagnostic_path,
        "run_json_path": run_json_path,
    }
