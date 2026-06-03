from __future__ import annotations

from bijux_phylogenetics.ancestral import (
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative import (
    run_pgls,
    summarize_brownian_trait_evolution,
    summarize_ou_trait_evolution,
    summarize_phylogenetic_signal,
)

from .models import PrimateComparativeWorkflowReport
from .panel import (
    SIGNAL_PERMUTATIONS,
    SIGNAL_SEED,
    load_primate_comparative_dataset,
)


def run_primate_comparative_workflow() -> PrimateComparativeWorkflowReport:
    """Run the owned comparative workflow over the packaged real primate dataset."""
    dataset = load_primate_comparative_dataset()
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
    return PrimateComparativeWorkflowReport(
        dataset=dataset,
        pgls=pgls,
        brownian=brownian,
        ou=ou,
        signal=signal,
        continuous_ancestral=continuous_ancestral,
        discrete_ancestral=discrete_ancestral,
    )
