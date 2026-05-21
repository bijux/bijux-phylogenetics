from __future__ import annotations

from pathlib import Path
import shutil

from .models import KnownAnswerReferenceExportResult
from .panel import load_known_answer_reference_dataset


def export_known_answer_reference_dataset(
    destination: Path,
) -> KnownAnswerReferenceExportResult:
    """Copy the packaged known-answer simulation dataset and reference outputs."""
    dataset = load_known_answer_reference_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    true_tree_path = shutil.copy2(dataset.true_tree_path, destination / "true-tree.nwk")
    alignment_path = shutil.copy2(
        dataset.alignment_path,
        destination / "simulated-alignment.fasta",
    )
    continuous_traits_path = shutil.copy2(
        dataset.continuous_traits_path,
        destination / "continuous-traits.tsv",
    )
    ou_traits_path = shutil.copy2(
        dataset.ou_traits_path,
        destination / "ou-traits.tsv",
    )
    discrete_traits_path = shutil.copy2(
        dataset.discrete_traits_path,
        destination / "discrete-traits.tsv",
    )
    host_traits_path = shutil.copy2(
        dataset.host_traits_path,
        destination / "host-traits.tsv",
    )
    geographic_traits_path = shutil.copy2(
        dataset.geographic_traits_path,
        destination / "geographic-traits.tsv",
    )
    true_parameters_path = shutil.copy2(
        dataset.true_parameters_path,
        destination / "true-parameters.tsv",
    )
    true_continuous_nodes_path = shutil.copy2(
        dataset.true_continuous_nodes_path,
        destination / "true-continuous-nodes.tsv",
    )
    true_ou_nodes_path = shutil.copy2(
        dataset.true_ou_nodes_path,
        destination / "true-ou-nodes.tsv",
    )
    true_discrete_nodes_path = shutil.copy2(
        dataset.true_discrete_nodes_path,
        destination / "true-discrete-nodes.tsv",
    )
    true_host_nodes_path = shutil.copy2(
        dataset.true_host_nodes_path,
        destination / "true-host-nodes.tsv",
    )
    true_geographic_nodes_path = shutil.copy2(
        dataset.true_geographic_nodes_path,
        destination / "true-geographic-nodes.tsv",
    )
    true_host_switch_events_path = shutil.copy2(
        dataset.true_host_switch_events_path,
        destination / "true-host-switch-events.tsv",
    )
    true_geographic_transition_events_path = shutil.copy2(
        dataset.true_geographic_transition_events_path,
        destination / "true-geographic-transition-events.tsv",
    )
    recovery_thresholds_path = shutil.copy2(
        dataset.recovery_thresholds_path,
        destination / "recovery-thresholds.tsv",
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return KnownAnswerReferenceExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        true_tree_path=Path(true_tree_path),
        alignment_path=Path(alignment_path),
        continuous_traits_path=Path(continuous_traits_path),
        ou_traits_path=Path(ou_traits_path),
        discrete_traits_path=Path(discrete_traits_path),
        host_traits_path=Path(host_traits_path),
        geographic_traits_path=Path(geographic_traits_path),
        true_parameters_path=Path(true_parameters_path),
        true_continuous_nodes_path=Path(true_continuous_nodes_path),
        true_ou_nodes_path=Path(true_ou_nodes_path),
        true_discrete_nodes_path=Path(true_discrete_nodes_path),
        true_host_nodes_path=Path(true_host_nodes_path),
        true_geographic_nodes_path=Path(true_geographic_nodes_path),
        true_host_switch_events_path=Path(true_host_switch_events_path),
        true_geographic_transition_events_path=Path(
            true_geographic_transition_events_path
        ),
        recovery_thresholds_path=Path(recovery_thresholds_path),
        expected_output_root=expected_output_root,
    )
