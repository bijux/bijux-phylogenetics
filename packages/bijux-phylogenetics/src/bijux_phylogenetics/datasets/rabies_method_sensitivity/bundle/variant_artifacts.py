from __future__ import annotations

from pathlib import Path
import shutil

from ..models import (
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariantRun,
)
from .shared import _format_optional_bool, _write_tsv


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _copy_task_logs(
    output_root: Path,
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...],
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for task in task_records:
        _copy_output(task.log_path, output_root / task.log_path.name)
    return output_root


def _write_variant_outputs(
    output_root: Path, variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for variant in variant_runs:
        variant_root = output_root / variant.config.variant_id
        variant_root.mkdir(parents=True, exist_ok=True)
        _copy_output(
            variant.alignment_workflow.output_paths["alignment"],
            variant_root / f"{variant.config.variant_id}.aln",
        )
        _copy_output(
            variant.trimming_workflow.output_paths["trimmed_alignment"],
            variant_root / f"{variant.config.variant_id}.trimmed.aln",
        )
        _copy_output(
            variant.inference_comparison.output_paths["fasttree_tree"],
            variant_root / "fasttree.nwk",
        )
        _copy_output(
            variant.inference_comparison.output_paths["iqtree_support_tree"],
            variant_root / "iqtree-support.nwk",
        )
        _copy_output(
            variant.rooted_fasttree_path,
            variant_root / "rooted-fasttree.nwk",
        )
        _copy_output(
            variant.rooted_iqtree_path,
            variant_root / "rooted-iqtree-support.nwk",
        )
        _write_rooting_summary_table(
            variant_root / "rooting-summary.tsv",
            variant,
        )
        _copy_output(
            variant.inference_comparison.output_paths["stability_summary"],
            variant_root / "unrooted-stability-summary.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conclusion_table"],
            variant_root / "unrooted-conclusions.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["support_weighted_conflicts"],
            variant_root / "unrooted-support-weighted-conflicts.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["shared_clades"],
            variant_root / "unrooted-shared-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conflicting_clades"],
            variant_root / "unrooted-conflicting-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["comparison_table"],
            variant_root / "unrooted-comparison.tsv",
        )
        _copy_output(
            variant.rooted_engine_comparison_table_path,
            variant_root / "rooted-engine-comparison.tsv",
        )
    return output_root


def _write_rooting_summary_table(
    path: Path, variant: RabiesMethodSensitivityVariantRun
) -> Path:
    rows = [
        [
            "engine_name",
            "requested_taxa",
            "matched_taxa",
            "outgroup_monophyletic",
            "rooted_outgroup_taxa",
            "warning_count",
        ],
        [
            "fasttree",
            ",".join(variant.fasttree_rooting.requested_taxa),
            ",".join(variant.fasttree_rooting.matched_taxa),
            _format_optional_bool(variant.fasttree_rooting.outgroup_monophyletic),
            ",".join(variant.fasttree_rooting.rooted_outgroup_taxa),
            str(len(variant.fasttree_rooting.warnings)),
        ],
        [
            "iqtree",
            ",".join(variant.iqtree_rooting.requested_taxa),
            ",".join(variant.iqtree_rooting.matched_taxa),
            _format_optional_bool(variant.iqtree_rooting.outgroup_monophyletic),
            ",".join(variant.iqtree_rooting.rooted_outgroup_taxa),
            str(len(variant.iqtree_rooting.warnings)),
        ],
    ]
    return _write_tsv(path, rows)
