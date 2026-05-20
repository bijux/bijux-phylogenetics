from __future__ import annotations

import json
import math
from pathlib import Path
import statistics
import tempfile

from Bio import Phylo

from bijux_phylogenetics.ancestral.common import (
    load_continuous_dataset,
    load_discrete_dataset,
    stable_value,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AncestralReconstructionError,
    InvalidAlignmentError,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.trees import load_tree_set
from .models import (
    AncestralTreeSetExclusion,
    AncestralTreeSetTreeRow,
    ContinuousAncestralTreeSetCladeSummaryRow,
    ContinuousAncestralTreeSetNodeRow,
    ContinuousAncestralTreeSetReport,
    ContinuousAncestralTreeSetSummary,
    DiscreteAncestralTreeSetCladeSummaryRow,
    DiscreteAncestralTreeSetNodeRow,
    DiscreteAncestralTreeSetReport,
    DiscreteAncestralTreeSetSummary,
)


def summarize_continuous_ancestral_tree_set(
    tree_set_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    burnin_fraction: float = 0.0,
) -> ContinuousAncestralTreeSetReport:
    """Run continuous ancestral reconstruction across a retained posterior or bootstrap tree set."""
    if model not in {"brownian", "ou"}:
        raise AncestralReconstructionError(
            f"unsupported continuous ancestral tree-set model: {model}"
        )
    if alpha <= 0:
        raise AncestralReconstructionError(
            "continuous ancestral tree-set alpha must be positive"
        )
    _validate_burnin_fraction(burnin_fraction)
    _source_format, trees = _load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = math.floor(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise AncestralReconstructionError(
            "ancestral tree-set analysis retains no trees after burn-in removal"
        )
    kept_trees = [tree for _, tree in kept_tree_entries]
    shared_tree_taxa = sorted(_shared_taxa(kept_trees))
    warnings: list[str] = []
    if any(set(tree.tip_names) != set(shared_tree_taxa) for tree in kept_trees):
        warnings.append(
            "retained trees do not share identical tip sets and were reduced to their shared taxa"
        )
    (
        analysis_trees,
        topology_summary,
        analysis_taxa,
        exclusions,
        dataset_warnings,
        resolved_taxon_column,
    ) = _prepare_analysis_tree_set(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        trait=trait,
        kept_tree_entries=kept_tree_entries,
        shared_tree_taxa=shared_tree_taxa,
        dataset_kind="continuous",
    )
    warnings.extend(dataset_warnings)
    tree_rows = [
        AncestralTreeSetTreeRow(
            source_tree_index=source_tree_index,
            post_burnin_index=record.index,
            rooted_topology_id=record.rooted_topology_id,
            unrooted_topology_id=record.unrooted_topology_id,
            internal_clade_count=max(len(analysis_tree.tip_names) - 1, 0),
        )
        for (source_tree_index, analysis_tree), record in zip(
            analysis_trees, topology_summary.records, strict=True
        )
    ]
    node_rows: list[ContinuousAncestralTreeSetNodeRow] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-tree-set-continuous-"
    ) as tmp_dir:
        current_tree_path = Path(tmp_dir) / "ancestral-tree-set-current-tree.nwk"
        for tree_row, (_source_tree_index, analysis_tree) in zip(
            tree_rows,
            analysis_trees,
            strict=True,
        ):
            current_tree_path.write_text(
                dumps_newick(analysis_tree) + "\n",
                encoding="utf-8",
            )
            report = reconstruct_continuous_ancestral_states(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                alpha=alpha,
            )
            for estimate in report.estimates:
                if estimate.is_tip:
                    continue
                node_rows.append(
                    ContinuousAncestralTreeSetNodeRow(
                        source_tree_index=tree_row.source_tree_index,
                        post_burnin_index=tree_row.post_burnin_index,
                        rooted_topology_id=tree_row.rooted_topology_id,
                        unrooted_topology_id=tree_row.unrooted_topology_id,
                        clade_id=estimate.node,
                        clade_taxa=estimate.descendant_taxa,
                        estimate=estimate.estimate,
                        standard_error=estimate.standard_error,
                        lower_95_interval=estimate.lower_95_interval,
                        upper_95_interval=estimate.upper_95_interval,
                        confidence=estimate.confidence,
                        unstable=estimate.unstable,
                    )
                )
    clade_summaries = _summarize_continuous_clades(
        node_rows, kept_tree_count=len(tree_rows)
    )
    if any(row.tree_presence_fraction < 1.0 for row in clade_summaries):
        warnings.append(
            "one or more comparable ancestral clades are absent from some retained trees"
        )
    if any(row.stability_class != "stable" for row in clade_summaries):
        warnings.append(
            "one or more continuous ancestral clades show topology-sensitive or dispersed values across retained trees"
        )
    return ContinuousAncestralTreeSetReport(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=resolved_taxon_column,
        model=model,
        alpha=alpha,
        burnin_fraction=burnin_fraction,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(tree_rows),
        shared_tree_taxa=shared_tree_taxa,
        analysis_taxa=analysis_taxa,
        rooted_topology_count=topology_summary.rooted_topology_count,
        unrooted_topology_count=topology_summary.unrooted_topology_count,
        tree_rows=tree_rows,
        node_rows=node_rows,
        clade_summaries=clade_summaries,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_discrete_ancestral_tree_set(
    tree_set_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    burnin_fraction: float = 0.0,
) -> DiscreteAncestralTreeSetReport:
    """Run discrete ancestral reconstruction across a retained posterior or bootstrap tree set."""
    if model == "meristic":
        _resolve_discrete_model_name(model)
    if model not in {"fitch", "equal-rates", "symmetric", "all-rates-different"}:
        raise AncestralReconstructionError(
            f"unsupported discrete ancestral tree-set model: {model}"
        )
    if model == "fitch" and state_ordering != "unordered":
        raise AncestralReconstructionError(
            "ordered discrete ancestral tree-set reconstruction requires a likelihood model"
        )
    _validate_burnin_fraction(burnin_fraction)
    _source_format, trees = _load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = math.floor(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise AncestralReconstructionError(
            "ancestral tree-set analysis retains no trees after burn-in removal"
        )
    kept_trees = [tree for _, tree in kept_tree_entries]
    shared_tree_taxa = sorted(_shared_taxa(kept_trees))
    warnings: list[str] = []
    if any(set(tree.tip_names) != set(shared_tree_taxa) for tree in kept_trees):
        warnings.append(
            "retained trees do not share identical tip sets and were reduced to their shared taxa"
        )
    (
        analysis_trees,
        topology_summary,
        analysis_taxa,
        exclusions,
        dataset_warnings,
        resolved_taxon_column,
    ) = _prepare_analysis_tree_set(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        trait=trait,
        kept_tree_entries=kept_tree_entries,
        shared_tree_taxa=shared_tree_taxa,
        dataset_kind="discrete",
    )
    warnings.extend(dataset_warnings)
    tree_rows = [
        AncestralTreeSetTreeRow(
            source_tree_index=source_tree_index,
            post_burnin_index=record.index,
            rooted_topology_id=record.rooted_topology_id,
            unrooted_topology_id=record.unrooted_topology_id,
            internal_clade_count=max(len(analysis_tree.tip_names) - 1, 0),
        )
        for (source_tree_index, analysis_tree), record in zip(
            analysis_trees, topology_summary.records, strict=True
        )
    ]
    node_rows: list[DiscreteAncestralTreeSetNodeRow] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-tree-set-discrete-"
    ) as tmp_dir:
        current_tree_path = Path(tmp_dir) / "ancestral-tree-set-current-tree.nwk"
        for tree_row, (_source_tree_index, analysis_tree) in zip(
            tree_rows,
            analysis_trees,
            strict=True,
        ):
            current_tree_path.write_text(
                dumps_newick(analysis_tree) + "\n",
                encoding="utf-8",
            )
            report = reconstruct_discrete_ancestral_states(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            for estimate in report.estimates:
                if estimate.is_tip:
                    continue
                node_rows.append(
                    DiscreteAncestralTreeSetNodeRow(
                        source_tree_index=tree_row.source_tree_index,
                        post_burnin_index=tree_row.post_burnin_index,
                        rooted_topology_id=tree_row.rooted_topology_id,
                        unrooted_topology_id=tree_row.unrooted_topology_id,
                        clade_id=estimate.node,
                        clade_taxa=estimate.descendant_taxa,
                        most_likely_state=estimate.most_likely_state,
                        state_set=estimate.state_set,
                        confidence=estimate.confidence,
                        ambiguous=estimate.ambiguous,
                        unstable=estimate.unstable,
                    )
                )
    clade_summaries = _summarize_discrete_clades(
        node_rows, kept_tree_count=len(tree_rows)
    )
    if any(row.tree_presence_fraction < 1.0 for row in clade_summaries):
        warnings.append(
            "one or more comparable ancestral clades are absent from some retained trees"
        )
    if any(row.stability_class != "stable" for row in clade_summaries):
        warnings.append(
            "one or more discrete ancestral clades change state or support profile across retained trees"
        )
    return DiscreteAncestralTreeSetReport(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=resolved_taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        burnin_fraction=burnin_fraction,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(tree_rows),
        shared_tree_taxa=shared_tree_taxa,
        analysis_taxa=analysis_taxa,
        rooted_topology_count=topology_summary.rooted_topology_count,
        unrooted_topology_count=topology_summary.unrooted_topology_count,
        tree_rows=tree_rows,
        node_rows=node_rows,
        clade_summaries=clade_summaries,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_continuous_ancestral_tree_set_report(
    report: ContinuousAncestralTreeSetReport,
) -> ContinuousAncestralTreeSetSummary:
    """Summarize the main review facts for one continuous ancestral tree-set report."""
    unstable_clades = [
        row for row in report.clade_summaries if row.stability_class != "stable"
    ]
    top_unstable = (
        max(
            unstable_clades,
            key=lambda row: (row.instability_score, row.clade_id),
        ).clade_id
        if unstable_clades
        else None
    )
    return ContinuousAncestralTreeSetSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        alpha=report.alpha,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        kept_tree_count=report.kept_tree_count,
        shared_tree_taxon_count=len(report.shared_tree_taxa),
        analysis_taxon_count=len(report.analysis_taxa),
        rooted_topology_count=report.rooted_topology_count,
        unrooted_topology_count=report.unrooted_topology_count,
        clade_summary_count=len(report.clade_summaries),
        unstable_clade_count=len(unstable_clades),
        top_unstable_clade=top_unstable,
        warning_count=len(report.warnings),
    )


def summarize_discrete_ancestral_tree_set_report(
    report: DiscreteAncestralTreeSetReport,
) -> DiscreteAncestralTreeSetSummary:
    """Summarize the main review facts for one discrete ancestral tree-set report."""
    unstable_clades = [
        row for row in report.clade_summaries if row.stability_class != "stable"
    ]
    top_unstable = (
        max(
            unstable_clades,
            key=lambda row: (row.instability_score, row.clade_id),
        ).clade_id
        if unstable_clades
        else None
    )
    return DiscreteAncestralTreeSetSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        kept_tree_count=report.kept_tree_count,
        shared_tree_taxon_count=len(report.shared_tree_taxa),
        analysis_taxon_count=len(report.analysis_taxa),
        rooted_topology_count=report.rooted_topology_count,
        unrooted_topology_count=report.unrooted_topology_count,
        clade_summary_count=len(report.clade_summaries),
        unstable_clade_count=len(unstable_clades),
        top_unstable_clade=top_unstable,
        warning_count=len(report.warnings),
    )


def write_ancestral_tree_set_tree_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport | DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one retained-tree ledger for an ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "internal_clade_count",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "internal_clade_count": str(row.internal_clade_count),
            }
            for row in report.tree_rows
        ],
    )


def write_ancestral_tree_set_exclusion_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport | DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one excluded-taxa ledger for an ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.exclusions
        ],
    )


def write_continuous_ancestral_tree_set_summary_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one overall summary ledger for continuous ancestral tree-set analysis."""
    summary = summarize_continuous_ancestral_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "alpha",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "clade_summary_count",
            "unstable_clade_count",
            "top_unstable_clade",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "alpha": str(summary.alpha),
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "clade_summary_count": str(summary.clade_summary_count),
                "unstable_clade_count": str(summary.unstable_clade_count),
                "top_unstable_clade": summary.top_unstable_clade or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_tree_set_node_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one per-tree internal-node ledger for continuous ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "clade_id",
            "clade_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "confidence",
            "unstable",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "estimate": str(row.estimate),
                "standard_error": str(row.standard_error),
                "lower_95_interval": str(row.lower_95_interval),
                "upper_95_interval": str(row.upper_95_interval),
                "confidence": str(row.confidence),
                "unstable": str(row.unstable).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_continuous_ancestral_tree_set_clade_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one comparable-clade summary ledger for continuous ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "mean_estimate",
            "median_estimate",
            "standard_deviation",
            "minimum_estimate",
            "maximum_estimate",
            "lower_95_empirical_estimate",
            "upper_95_empirical_estimate",
            "empirical_interval_width",
            "mean_standard_error",
            "unstable_tree_count",
            "unstable_tree_fraction",
            "instability_score",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "mean_estimate": str(row.mean_estimate),
                "median_estimate": str(row.median_estimate),
                "standard_deviation": str(row.standard_deviation),
                "minimum_estimate": str(row.minimum_estimate),
                "maximum_estimate": str(row.maximum_estimate),
                "lower_95_empirical_estimate": str(row.lower_95_empirical_estimate),
                "upper_95_empirical_estimate": str(row.upper_95_empirical_estimate),
                "empirical_interval_width": str(row.empirical_interval_width),
                "mean_standard_error": str(row.mean_standard_error),
                "unstable_tree_count": str(row.unstable_tree_count),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "instability_score": str(row.instability_score),
                "stability_class": row.stability_class,
            }
            for row in report.clade_summaries
        ],
    )


def write_discrete_ancestral_tree_set_summary_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one overall summary ledger for discrete ancestral tree-set analysis."""
    summary = summarize_discrete_ancestral_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "clade_summary_count",
            "unstable_clade_count",
            "top_unstable_clade",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "clade_summary_count": str(summary.clade_summary_count),
                "unstable_clade_count": str(summary.unstable_clade_count),
                "top_unstable_clade": summary.top_unstable_clade or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_discrete_ancestral_tree_set_node_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one per-tree internal-node ledger for discrete ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "clade_id",
            "clade_taxa",
            "most_likely_state",
            "state_set",
            "confidence",
            "ambiguous",
            "unstable",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "most_likely_state": row.most_likely_state,
                "state_set": ",".join(row.state_set),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "unstable": str(row.unstable).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_discrete_ancestral_tree_set_clade_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one comparable-clade summary ledger for discrete ancestral tree-set analysis."""
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "unique_state_count",
            "dominant_state",
            "dominant_state_tree_count",
            "dominant_state_fraction",
            "ambiguous_tree_count",
            "ambiguous_tree_fraction",
            "unstable_tree_count",
            "unstable_tree_fraction",
            "state_distribution",
            "instability_score",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "unique_state_count": str(row.unique_state_count),
                "dominant_state": row.dominant_state,
                "dominant_state_tree_count": str(row.dominant_state_tree_count),
                "dominant_state_fraction": str(row.dominant_state_fraction),
                "ambiguous_tree_count": str(row.ambiguous_tree_count),
                "ambiguous_tree_fraction": str(row.ambiguous_tree_fraction),
                "unstable_tree_count": str(row.unstable_tree_count),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "state_distribution": json.dumps(
                    row.state_distribution,
                    sort_keys=True,
                ),
                "instability_score": str(row.instability_score),
                "stability_class": row.stability_class,
            }
            for row in report.clade_summaries
        ],
    )


def _prepare_analysis_tree_set(
    *,
    tree_set_path: Path,
    traits_path: Path,
    taxon_column: str | None,
    trait: str,
    kept_tree_entries: list[tuple[int, PhyloTree]],
    shared_tree_taxa: list[str],
    dataset_kind: str,
) -> tuple[
    list[tuple[int, PhyloTree]],
    object,
    list[str],
    list[AncestralTreeSetExclusion],
    list[str],
    str,
]:
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-tree-set-prepare-"
    ) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        reference_tree_path = tmp_dir_path / "ancestral-tree-set-reference.nwk"
        reference_tree = _prune_tree_to_taxa(
            kept_tree_entries[0][1],
            shared_tree_taxa,
            scratch_path=reference_tree_path,
        )
        reference_tree_path.write_text(
            dumps_newick(reference_tree) + "\n",
            encoding="utf-8",
        )
        if dataset_kind == "continuous":
            dataset = load_continuous_dataset(
                reference_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
            )
            exclusions = [
                *[
                    AncestralTreeSetExclusion(
                        taxon=taxon,
                        reason="missing_trait_value",
                    )
                    for taxon in dataset.dropped_missing_taxa
                ],
                *[
                    AncestralTreeSetExclusion(
                        taxon=taxon,
                        reason="non_numeric_trait_value",
                    )
                    for taxon in dataset.dropped_non_numeric_taxa
                ],
            ]
        else:
            dataset = load_discrete_dataset(
                reference_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
            )
            exclusions = [
                AncestralTreeSetExclusion(
                    taxon=taxon,
                    reason="missing_discrete_trait_state",
                )
                for taxon in dataset.dropped_missing_taxa
            ]
        analysis_taxa = list(dataset.taxa)
        analysis_tree_set_path = tmp_dir_path / "ancestral-tree-set-analysis.nwk"
        analysis_trees = [
            (
                source_tree_index,
                _prune_tree_to_taxa(
                    tree,
                    analysis_taxa,
                    scratch_path=tmp_dir_path
                    / f"ancestral-tree-set-pruned-{source_tree_index}.nwk",
                ),
            )
            for source_tree_index, tree in kept_tree_entries
        ]
        analysis_tree_set_path.write_text(
            "".join(dumps_newick(tree) + "\n" for _, tree in analysis_trees),
            encoding="utf-8",
        )
        topology_summary = load_tree_set(analysis_tree_set_path)
        return (
            analysis_trees,
            topology_summary,
            analysis_taxa,
            sorted(exclusions, key=lambda row: (row.taxon, row.reason)),
            list(dataset.warnings),
            dataset.taxon_column,
        )


def _summarize_continuous_clades(
    rows: list[ContinuousAncestralTreeSetNodeRow],
    *,
    kept_tree_count: int,
) -> list[ContinuousAncestralTreeSetCladeSummaryRow]:
    grouped: dict[str, list[ContinuousAncestralTreeSetNodeRow]] = {}
    for row in rows:
        grouped.setdefault(row.clade_id, []).append(row)
    all_estimates = [row.estimate for row in rows]
    global_range = (
        max(all_estimates) - min(all_estimates) if len(all_estimates) > 1 else 0.0
    )
    scale = global_range if global_range > 0 else 1.0
    summaries: list[ContinuousAncestralTreeSetCladeSummaryRow] = []
    for clade_id, clade_rows in sorted(grouped.items()):
        estimates = [row.estimate for row in clade_rows]
        standard_errors = [row.standard_error for row in clade_rows]
        presence_fraction = stable_value(len(clade_rows) / kept_tree_count)
        unstable_tree_count = sum(row.unstable for row in clade_rows)
        unstable_tree_fraction = stable_value(unstable_tree_count / len(clade_rows))
        empirical_low = _empirical_quantile(estimates, 0.025)
        empirical_high = _empirical_quantile(estimates, 0.975)
        empirical_width = stable_value(empirical_high - empirical_low)
        mean_standard_error = stable_value(statistics.fmean(standard_errors))
        normalized_dispersion = empirical_width / scale
        instability_score = stable_value(
            (1.0 - presence_fraction) + unstable_tree_fraction + normalized_dispersion
        )
        if presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif unstable_tree_count > 0:
            stability_class = "within_tree_uncertainty"
        elif normalized_dispersion > 0.5:
            stability_class = "value_dispersion"
        else:
            stability_class = "stable"
        summaries.append(
            ContinuousAncestralTreeSetCladeSummaryRow(
                clade_id=clade_id,
                clade_taxa=clade_rows[0].clade_taxa,
                tree_presence_count=len(clade_rows),
                tree_presence_fraction=presence_fraction,
                mean_estimate=stable_value(statistics.fmean(estimates)),
                median_estimate=stable_value(statistics.median(estimates)),
                standard_deviation=stable_value(_sample_standard_deviation(estimates)),
                minimum_estimate=stable_value(min(estimates)),
                maximum_estimate=stable_value(max(estimates)),
                lower_95_empirical_estimate=stable_value(empirical_low),
                upper_95_empirical_estimate=stable_value(empirical_high),
                empirical_interval_width=empirical_width,
                mean_standard_error=mean_standard_error,
                unstable_tree_count=unstable_tree_count,
                unstable_tree_fraction=unstable_tree_fraction,
                instability_score=instability_score,
                stability_class=stability_class,
            )
        )
    return summaries


def _summarize_discrete_clades(
    rows: list[DiscreteAncestralTreeSetNodeRow],
    *,
    kept_tree_count: int,
) -> list[DiscreteAncestralTreeSetCladeSummaryRow]:
    grouped: dict[str, list[DiscreteAncestralTreeSetNodeRow]] = {}
    for row in rows:
        grouped.setdefault(row.clade_id, []).append(row)
    summaries: list[DiscreteAncestralTreeSetCladeSummaryRow] = []
    for clade_id, clade_rows in sorted(grouped.items()):
        presence_fraction = stable_value(len(clade_rows) / kept_tree_count)
        state_counts: dict[str, int] = {}
        for row in clade_rows:
            state_counts[row.most_likely_state] = (
                state_counts.get(row.most_likely_state, 0) + 1
            )
        dominant_state = max(
            sorted(state_counts),
            key=lambda state: (state_counts[state], state),
        )
        dominant_state_tree_count = state_counts[dominant_state]
        dominant_state_fraction = stable_value(
            dominant_state_tree_count / len(clade_rows)
        )
        ambiguous_tree_count = sum(row.ambiguous for row in clade_rows)
        unstable_tree_count = sum(row.unstable for row in clade_rows)
        ambiguous_tree_fraction = stable_value(ambiguous_tree_count / len(clade_rows))
        unstable_tree_fraction = stable_value(unstable_tree_count / len(clade_rows))
        instability_score = stable_value(
            (1.0 - presence_fraction)
            + (1.0 - dominant_state_fraction)
            + ambiguous_tree_fraction
            + unstable_tree_fraction
        )
        if presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif len(state_counts) > 1:
            stability_class = "state_conflict"
        elif ambiguous_tree_count > 0 or unstable_tree_count > 0:
            stability_class = "low_confidence"
        else:
            stability_class = "stable"
        summaries.append(
            DiscreteAncestralTreeSetCladeSummaryRow(
                clade_id=clade_id,
                clade_taxa=clade_rows[0].clade_taxa,
                tree_presence_count=len(clade_rows),
                tree_presence_fraction=presence_fraction,
                unique_state_count=len(state_counts),
                dominant_state=dominant_state,
                dominant_state_tree_count=dominant_state_tree_count,
                dominant_state_fraction=dominant_state_fraction,
                ambiguous_tree_count=ambiguous_tree_count,
                ambiguous_tree_fraction=ambiguous_tree_fraction,
                unstable_tree_count=unstable_tree_count,
                unstable_tree_fraction=unstable_tree_fraction,
                state_distribution=dict(sorted(state_counts.items())),
                instability_score=instability_score,
                stability_class=stability_class,
            )
        )
    return summaries


def _validate_burnin_fraction(burnin_fraction: float) -> None:
    if not 0.0 <= burnin_fraction < 1.0:
        raise AncestralReconstructionError(
            "ancestral tree-set burnin fraction must be between 0 inclusive and 1 exclusive"
        )


def _load_tree_set_trees(path: Path) -> tuple[str, list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    return source_format, [
        tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees
    ]


def _shared_taxa(trees: list[PhyloTree]) -> set[str]:
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def _prune_tree_to_taxa(
    tree: PhyloTree,
    requested_taxa: list[str],
    *,
    scratch_path: Path,
) -> PhyloTree:
    scratch_path.write_text(dumps_newick(tree) + "\n", encoding="utf-8")
    pruned_tree, _report = prune_tree_to_requested_taxa(scratch_path, requested_taxa)
    return pruned_tree


def _empirical_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return stable_value(ordered[0])
    index = (len(ordered) - 1) * probability
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return stable_value(ordered[lower])
    fraction = index - lower
    return stable_value(ordered[lower] + (ordered[upper] - ordered[lower]) * fraction)


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return stable_value(statistics.stdev(values))
