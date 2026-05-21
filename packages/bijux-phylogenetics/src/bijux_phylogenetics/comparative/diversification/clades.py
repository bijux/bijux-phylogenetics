from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import TreeNode

from .models import (
    CladeDiversificationObservation,
    CladeDiversificationScanReport,
)
from .rates import estimate_diversification_rate
from .trees import descendant_taxa, node_age, node_depths


def detect_diversification_outlier_clades(
    tree_path: Path,
    *,
    min_tip_count: int = 2,
    model: str = "birth-death",
) -> CladeDiversificationScanReport:
    """Flag clades whose diversification rate is high or low relative to the tree-wide baseline."""
    global_report = estimate_diversification_rate(tree_path, model=model)
    tree = load_tree(tree_path)
    depths = node_depths(tree)
    observations: list[CladeDiversificationObservation] = []
    raw_rows: list[tuple[TreeNode, list[str], float, float]] = []
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        clade_taxa = descendant_taxa(node)
        if len(clade_taxa) < min_tip_count:
            continue
        crown_age = node_age(tree, depths, node)
        if crown_age <= 0.0:
            continue
        diversification_rate = float(
            format(math.log(len(clade_taxa)) / crown_age, ".15g")
        )
        raw_rows.append((node, clade_taxa, crown_age, diversification_rate))
    rates = [row[3] for row in raw_rows]
    mean_rate = sum(rates) / max(len(rates), 1)
    variance = sum((rate - mean_rate) ** 2 for rate in rates) / max(len(rates), 1)
    standard_deviation = math.sqrt(variance)
    for node, clade_taxa, crown_age, diversification_rate in raw_rows:
        z_score = (
            float(
                format((diversification_rate - mean_rate) / standard_deviation, ".15g")
            )
            if standard_deviation > 0.0
            else 0.0
        )
        if z_score >= 1.0:
            classification = "high"
        elif z_score <= -1.0:
            classification = "low"
        else:
            classification = "baseline"
        observations.append(
            CladeDiversificationObservation(
                node=node_signature(node),
                node_name=node.name,
                descendant_taxa=clade_taxa,
                tip_count=len(clade_taxa),
                crown_age=crown_age,
                diversification_rate=diversification_rate,
                z_score=z_score,
                classification=classification,
            )
        )
    high = [row for row in observations if row.classification == "high"]
    low = [row for row in observations if row.classification == "low"]
    return CladeDiversificationScanReport(
        tree_path=tree_path,
        model=model,
        global_rate=global_report.net_diversification_rate,
        observations=observations,
        high_diversification_clades=high,
        low_diversification_clades=low,
        warnings=list(global_report.warnings),
    )


def write_clade_diversification_table(
    path: Path, report: CladeDiversificationScanReport
) -> Path:
    """Export clade diversification summaries as a deterministic TSV."""
    rows = [
        {
            "node": row.node,
            "node_name": row.node_name or "",
            "descendant_taxa": ",".join(row.descendant_taxa),
            "tip_count": str(row.tip_count),
            "crown_age": format(row.crown_age, ".15g"),
            "diversification_rate": format(row.diversification_rate, ".15g"),
            "z_score": format(row.z_score, ".15g"),
            "classification": row.classification,
        }
        for row in report.observations
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "tip_count",
            "crown_age",
            "diversification_rate",
            "z_score",
            "classification",
        ],
        rows=rows,
    )
