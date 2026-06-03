from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table

from .models import PrimateComparativeDataset

DATASET_ID = "primate_comparative"
DATASET_LABEL = "Primate comparative mammal dataset"
TAXON_COLUMN = "species"
CONTINUOUS_TRAIT = "longevity"
PGLS_PREDICTOR = "social_group_size"
DISCRETE_TRAIT = "mating_system"
SIGNAL_PERMUTATIONS = 11
SIGNAL_SEED = 7


def load_primate_comparative_dataset() -> PrimateComparativeDataset:
    """Expose the packaged real primate dataset as a first-class runtime surface."""
    dataset_root = _resource_root()
    traits_path = dataset_root / "traits.csv"
    table = load_taxon_table(traits_path, taxon_column=TAXON_COLUMN)
    return PrimateComparativeDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        tree_path=dataset_root / "tree.nwk",
        traits_path=traits_path,
        reference_output_root=dataset_root / "expected",
        taxon_column=TAXON_COLUMN,
        taxon_count=table.row_count,
        continuous_traits=(
            "body_mass",
            "gestation",
            "home_range",
            "longevity",
            "social_group_size",
        ),
        categorical_traits=("family", "sex_dimorphism", "mating_system"),
        workflow_continuous_trait=CONTINUOUS_TRAIT,
        workflow_pgls_predictor=PGLS_PREDICTOR,
        workflow_discrete_trait=DISCRETE_TRAIT,
        source_locator=(
            "evidence-book/studies/primate-longevity-signal/datasets/"
            "reference_primate.csv"
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "mammals"
        / DATASET_ID
    )
