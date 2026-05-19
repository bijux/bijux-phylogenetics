from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.trees import load_tree

from .models import (
    DATASET_ID,
    DATASET_LABEL,
    DISTANCE_METHOD,
    DISTANCE_MODEL,
    SEQUENCE_TYPE,
    KnownAnswerReferenceDataset,
)


def load_known_answer_reference_dataset() -> KnownAnswerReferenceDataset:
    """Expose the packaged deterministic simulation panel as a first-class surface."""
    dataset_root = _resource_root()
    true_tree_path = dataset_root / "true-tree.nwk"
    alignment_path = dataset_root / "simulated-alignment.fasta"
    continuous_traits_path = dataset_root / "continuous-traits.tsv"
    ou_traits_path = dataset_root / "ou-traits.tsv"
    discrete_traits_path = dataset_root / "discrete-traits.tsv"
    host_traits_path = dataset_root / "host-traits.tsv"
    geographic_traits_path = dataset_root / "geographic-traits.tsv"
    true_parameters_path = dataset_root / "true-parameters.tsv"
    true_continuous_nodes_path = dataset_root / "true-continuous-nodes.tsv"
    true_ou_nodes_path = dataset_root / "true-ou-nodes.tsv"
    true_discrete_nodes_path = dataset_root / "true-discrete-nodes.tsv"
    true_host_nodes_path = dataset_root / "true-host-nodes.tsv"
    true_geographic_nodes_path = dataset_root / "true-geographic-nodes.tsv"
    true_host_switch_events_path = dataset_root / "true-host-switch-events.tsv"
    true_geographic_transition_events_path = (
        dataset_root / "true-geographic-transition-events.tsv"
    )
    recovery_thresholds_path = dataset_root / "recovery-thresholds.tsv"
    validate_fasta_input(alignment_path, sequence_type=SEQUENCE_TYPE)
    records = load_fasta_alignment(alignment_path)
    tree = load_tree(true_tree_path)
    return KnownAnswerReferenceDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        true_tree_path=true_tree_path,
        alignment_path=alignment_path,
        continuous_traits_path=continuous_traits_path,
        ou_traits_path=ou_traits_path,
        discrete_traits_path=discrete_traits_path,
        host_traits_path=host_traits_path,
        geographic_traits_path=geographic_traits_path,
        true_parameters_path=true_parameters_path,
        true_continuous_nodes_path=true_continuous_nodes_path,
        true_ou_nodes_path=true_ou_nodes_path,
        true_discrete_nodes_path=true_discrete_nodes_path,
        true_host_nodes_path=true_host_nodes_path,
        true_geographic_nodes_path=true_geographic_nodes_path,
        true_host_switch_events_path=true_host_switch_events_path,
        true_geographic_transition_events_path=true_geographic_transition_events_path,
        recovery_thresholds_path=recovery_thresholds_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=tree.tip_count,
        sequence_length=len(records[0].sequence),
        sequence_type=SEQUENCE_TYPE,
        distance_method=DISTANCE_METHOD,
        distance_model=DISTANCE_MODEL,
        source_summary=(
            "Deterministic owned simulation panel with one birth-death tree, one "
            "JC-like DNA alignment, one Brownian continuous trait, one OU "
            "continuous trait, one generic discrete trait, one host-state trait, "
            "and one geographic-state trait, packaged with node-level truth, "
            "branch-event truth, and explicit recovery thresholds."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "simulation"
        / DATASET_ID
    )
