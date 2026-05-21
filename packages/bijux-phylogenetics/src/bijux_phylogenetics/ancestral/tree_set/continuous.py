from __future__ import annotations

import math
from pathlib import Path
import statistics
import tempfile

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .models import (
    AncestralTreeSetTreeRow,
    ContinuousAncestralTreeSetCladeSummaryRow,
    ContinuousAncestralTreeSetNodeRow,
    ContinuousAncestralTreeSetReport,
)
from .preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
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
    validate_burnin_fraction(burnin_fraction)
    _source_format, trees = load_tree_set_trees(tree_set_path)
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
    shared_tree_taxa = sorted(shared_taxa(kept_trees))
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
    ) = prepare_analysis_tree_set(
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
    clade_summaries = summarize_continuous_clades(
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


def summarize_continuous_clades(
    rows: list[ContinuousAncestralTreeSetNodeRow],
    *,
    kept_tree_count: int,
) -> list[ContinuousAncestralTreeSetCladeSummaryRow]:
    """Summarize comparable continuous clades across retained trees."""
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
        empirical_low = empirical_quantile(estimates, 0.025)
        empirical_high = empirical_quantile(estimates, 0.975)
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
                standard_deviation=stable_value(sample_standard_deviation(estimates)),
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


def empirical_quantile(values: list[float], probability: float) -> float:
    """Interpolate one empirical quantile from a stable ordered sample."""
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


def sample_standard_deviation(values: list[float]) -> float:
    """Return the sample standard deviation for one empirical distribution."""
    if len(values) < 2:
        return 0.0
    return stable_value(statistics.stdev(values))
