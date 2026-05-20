from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table

from .models import AvianReproductiveTraitDataset

DATASET_ID = "avian_reproductive_traits"
DATASET_LABEL = "Avian reproductive trait dataset"
TAXON_COLUMN = "species"
CONTINUOUS_TRAIT = "testes_mass"
PGLS_PREDICTOR = "body_mass"
DISCRETE_TRAIT = "mating_system"
CLADE_TRAIT = "mating_system"
SIGNAL_PERMUTATIONS = 11
SIGNAL_SEED = 7


def load_avian_reproductive_trait_dataset() -> AvianReproductiveTraitDataset:
    """Expose the packaged bird dataset as a first-class runtime surface."""
    dataset_root = _resource_root()
    traits_path = dataset_root / "traits.csv"
    table = load_taxon_table(traits_path, taxon_column=TAXON_COLUMN)
    return AvianReproductiveTraitDataset(
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
            "testes_mass",
            "multiple_paternity_percentage",
        ),
        categorical_traits=("mating_system", "development_mode"),
        workflow_continuous_trait=CONTINUOUS_TRAIT,
        workflow_pgls_predictor=PGLS_PREDICTOR,
        workflow_discrete_trait=DISCRETE_TRAIT,
        workflow_clade_trait=CLADE_TRAIT,
        source_summary=(
            "Curated avian reproductive trait table paired with a pruned rooted bird phylogeny."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "birds"
        / DATASET_ID
    )
