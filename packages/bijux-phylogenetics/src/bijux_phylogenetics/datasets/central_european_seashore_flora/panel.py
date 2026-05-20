from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table

from .models import CentralEuropeanSeashoreFloraDataset

DATASET_ID = "central_european_seashore_flora"
DATASET_LABEL = "Central European seashore flora dataset"
TAXON_COLUMN = "species"
CONTINUOUS_TRAIT = "seed_mass"
PGLS_PREDICTOR = "plant_height"
DISCRETE_TRAIT = "lifeform"
CLADE_TRAIT = "lifeform"
SIGNAL_PERMUTATIONS = 11
SIGNAL_SEED = 7


def load_central_european_seashore_flora_dataset() -> (
    CentralEuropeanSeashoreFloraDataset
):
    """Expose the packaged plant dataset as a first-class runtime surface."""
    dataset_root = _resource_root()
    traits_path = dataset_root / "traits.csv"
    table = load_taxon_table(traits_path, taxon_column=TAXON_COLUMN)
    return CentralEuropeanSeashoreFloraDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        tree_path=dataset_root / "tree.nwk",
        traits_path=traits_path,
        reference_output_root=dataset_root / "expected",
        taxon_column=TAXON_COLUMN,
        taxon_count=table.row_count,
        continuous_traits=(
            "seed_mass",
            "plant_height",
            "light",
            "temperature",
            "continentality",
            "humidity",
            "reaction",
            "nitrogen",
            "salt",
        ),
        categorical_traits=("lifeform", "longevity", "habitat"),
        workflow_continuous_trait=CONTINUOUS_TRAIT,
        workflow_pgls_predictor=PGLS_PREDICTOR,
        workflow_discrete_trait=DISCRETE_TRAIT,
        workflow_clade_trait=CLADE_TRAIT,
        source_summary=(
            "Published Central European flora dataset restricted to the source saltwater and seashore habitat class."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "plants"
        / DATASET_ID
    )
