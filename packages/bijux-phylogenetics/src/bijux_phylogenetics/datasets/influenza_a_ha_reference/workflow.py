from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.inference import run_fasta_to_tree_workflow

from .models import InfluenzaAHAReferenceWorkflowReport
from .panel import (
    BOOTSTRAP_REPLICATES,
    IQTREE_SEED,
    IQTREE_THREADS,
    load_influenza_a_ha_reference_dataset,
)


def run_influenza_a_ha_reference_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = IQTREE_SEED,
    iqtree_threads: int = IQTREE_THREADS,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> InfluenzaAHAReferenceWorkflowReport:
    """Run the owned FASTA-to-tree workflow over the packaged influenza A HA panel."""
    dataset = load_influenza_a_ha_reference_dataset()
    workflow = run_fasta_to_tree_workflow(
        dataset.sequences_path,
        out_dir=out_dir,
        prefix=dataset.workflow_prefix,
        sequence_type=dataset.sequence_type,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
    )
    return InfluenzaAHAReferenceWorkflowReport(dataset=dataset, workflow=workflow)
