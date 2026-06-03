from __future__ import annotations

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative import (
    run_pgls,
    summarize_brownian_trait_evolution,
    summarize_clade_traits,
    summarize_ou_trait_evolution,
    summarize_phylogenetic_signal,
)

from .models import CentralEuropeanSeashoreFloraWorkflowReport
from .panel import (
    SIGNAL_PERMUTATIONS,
    SIGNAL_SEED,
    load_central_european_seashore_flora_dataset,
)


def run_central_european_seashore_flora_workflow() -> (
    CentralEuropeanSeashoreFloraWorkflowReport
):
    """Run the owned comparative workflow over the packaged plant dataset."""
    dataset = load_central_european_seashore_flora_dataset()
    pgls = run_pgls(
        dataset.tree_path,
        dataset.traits_path,
        response=dataset.workflow_continuous_trait,
        predictors=[dataset.workflow_pgls_predictor],
        taxon_column=dataset.taxon_column,
        lambda_value="estimate",
    )
    brownian = summarize_brownian_trait_evolution(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
    )
    ou = summarize_ou_trait_evolution(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
    )
    signal = summarize_phylogenetic_signal(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        permutations=SIGNAL_PERMUTATIONS,
        seed=SIGNAL_SEED,
    )
    continuous_ancestral = reconstruct_continuous_ancestral_states(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        model="brownian",
    )
    discrete_ancestral = reconstruct_discrete_ancestral_states(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_discrete_trait,
        taxon_column=dataset.taxon_column,
        model="equal-rates",
    )
    clade_traits = summarize_clade_traits(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_clade_trait,
        taxon_column=dataset.taxon_column,
        trait_kind="categorical",
    )
    return CentralEuropeanSeashoreFloraWorkflowReport(
        dataset=dataset,
        pgls=pgls,
        brownian=brownian,
        ou=ou,
        signal=signal,
        continuous_ancestral=continuous_ancestral,
        discrete_ancestral=discrete_ancestral,
        clade_traits=clade_traits,
    )
