from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree

from .models import ContinuousModeRecoveryPanelDataset
from .scenarios import count_continuous_mode_recovery_panel_scenarios

DATASET_ID = "continuous_mode_recovery_panel"
DATASET_LABEL = "Continuous trait-model recovery panel"
DEFAULT_TREE_FILE = "trees/reference-tree-twelve-taxa.nwk"


def load_continuous_mode_recovery_panel_dataset() -> ContinuousModeRecoveryPanelDataset:
    """Expose the packaged continuous-mode recovery panel as a first-class surface."""
    dataset_root = _resource_root()
    default_tree_path = dataset_root / DEFAULT_TREE_FILE
    reference_tree_paths = sorted((dataset_root / "trees").glob("*.nwk"))
    simulation_cases_path = dataset_root / "simulation-cases.tsv"
    taxon_count = max(load_tree(path).tip_count for path in reference_tree_paths)
    return ContinuousModeRecoveryPanelDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        default_tree_path=default_tree_path,
        reference_tree_paths=reference_tree_paths,
        simulation_cases_path=simulation_cases_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        tree_count=len(reference_tree_paths),
        case_count=count_continuous_mode_recovery_panel_scenarios(
            simulation_cases_path,
            dataset_root,
        ),
        source_summary=(
            "Deterministic continuous-trait recovery panel with one governed "
            "twelve-taxon review tree for BM, OU, EB, and weak-OU identifiability, "
            "plus one governed twenty-four-taxon ultrametric review tree for "
            "Pagel-lambda, Pagel-kappa, and Pagel-delta transformed-branch "
            "recovery comparisons against stored local geiger references."
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
