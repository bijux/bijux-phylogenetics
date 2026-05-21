from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    compare_tree_paths,
    write_tree_comparison_table,
)
from bijux_phylogenetics.engines.inference import run_tree_inference_comparison
from bijux_phylogenetics.engines.workflows.alignment import (
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.topology import root_tree_on_outgroup

from ..models import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityVariant,
    RabiesMethodSensitivityVariantRun,
)


def _run_variant_workflow(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    variant: RabiesMethodSensitivityVariant,
    variant_root: Path,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    iqtree_seed: int,
    iqtree_threads: int,
    bootstrap_replicates: int,
) -> RabiesMethodSensitivityVariantRun:
    alignment_path = variant_root / f"{variant.variant_id}.aln"
    trimmed_alignment_path = variant_root / f"{variant.variant_id}.trimmed.aln"
    alignment_workflow = run_multiple_sequence_alignment(
        dataset.sequences_path,
        alignment_path,
        executable=mafft_executable,
        mode=variant.alignment_mode,
    )
    trimming_workflow = run_alignment_trimming(
        alignment_path,
        trimmed_alignment_path,
        executable=trimal_executable,
        mode=variant.trimming_mode,
        gap_threshold=variant.trim_gap_threshold,
    )
    inference_comparison = run_tree_inference_comparison(
        trimmed_alignment_path,
        out_dir=variant_root / "engine-comparison",
        prefix=variant.variant_id,
        sequence_type=dataset.sequence_type,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
    )
    fasttree_rooted, fasttree_rooting = root_tree_on_outgroup(
        inference_comparison.output_paths["fasttree_tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    iqtree_rooted, iqtree_rooting = root_tree_on_outgroup(
        inference_comparison.output_paths["iqtree_support_tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    rooted_fasttree_path = variant_root / "rooted-fasttree.nwk"
    rooted_iqtree_path = variant_root / "rooted-iqtree-support.nwk"
    write_newick(rooted_fasttree_path, fasttree_rooted)
    write_newick(rooted_iqtree_path, iqtree_rooted)
    rooted_engine_comparison = compare_tree_paths(
        rooted_fasttree_path, rooted_iqtree_path
    )
    rooted_engine_comparison_table_path = write_tree_comparison_table(
        variant_root / "rooted-engine-comparison.tsv",
        rooted_fasttree_path,
        rooted_iqtree_path,
    )
    aligned_records = load_fasta_alignment(alignment_workflow.output_paths["alignment"])
    trimmed_records = load_fasta_alignment(
        trimming_workflow.output_paths["trimmed_alignment"]
    )
    return RabiesMethodSensitivityVariantRun(
        config=variant,
        alignment_workflow=alignment_workflow,
        trimming_workflow=trimming_workflow,
        inference_comparison=inference_comparison,
        rooted_fasttree_path=rooted_fasttree_path,
        rooted_iqtree_path=rooted_iqtree_path,
        fasttree_rooting=fasttree_rooting,
        iqtree_rooting=iqtree_rooting,
        rooted_engine_comparison=rooted_engine_comparison,
        rooted_engine_comparison_table_path=rooted_engine_comparison_table_path,
        alignment_length=len(aligned_records[0].sequence),
        trimmed_alignment_length=len(trimmed_records[0].sequence),
    )
