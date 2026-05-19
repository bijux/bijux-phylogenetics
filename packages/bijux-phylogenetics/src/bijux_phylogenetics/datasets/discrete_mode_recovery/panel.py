from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree

from .models import DiscreteModeRecoveryPanelDataset
from .scenarios import count_discrete_mode_recovery_panel_scenarios

DATASET_ID = "discrete_mode_recovery_panel"
DATASET_LABEL = "Discrete trait-model recovery panel"
DEFAULT_TREE_FILE = "trees/reference-tree-twelve-taxa.nwk"


def load_discrete_mode_recovery_panel_dataset() -> DiscreteModeRecoveryPanelDataset:
    """Expose the packaged discrete-mode recovery panel as a first-class surface."""
    dataset_root = _resource_root()
    default_tree_path = dataset_root / DEFAULT_TREE_FILE
    reference_tree_paths = sorted((dataset_root / "trees").glob("*.nwk"))
    simulation_cases_path = dataset_root / "simulation-cases.tsv"
    taxon_count = max(load_tree(path).tip_count for path in reference_tree_paths)
    return DiscreteModeRecoveryPanelDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        default_tree_path=default_tree_path,
        reference_tree_paths=reference_tree_paths,
        simulation_cases_path=simulation_cases_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        tree_count=len(reference_tree_paths),
        case_count=count_discrete_mode_recovery_panel_scenarios(
            simulation_cases_path,
            dataset_root,
        ),
        source_summary=(
            "Deterministic discrete-trait recovery panel with one governed "
            "twelve-taxon overparameterized ARD failure tree plus one governed "
            "twenty-four-taxon rooted review tree for stable ER and SYM "
            "selection cases and one weak-identification ARD review surface "
            "compared against stored local geiger references."
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
