---
title: Native Bayesian Workflows
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Native Bayesian Workflows

Use the native Bayesian dispatcher when you want one owned posterior run
through `run_bayesian_inference(...)` and you are prepared to define the model
surface explicitly.

This public surface is Python-first today. The stable contract is one
dispatcher plus validated model-definition and proposal-schedule builders. The
dispatcher rejects unsupported or mismatched inputs explicitly instead of
quietly routing through external engines.

## Run One Fixed-Topology DNA Chain

Fixed-topology DNA workflows keep the supplied tree constant while sampling
branch lengths and supported substitution parameters.

```python
from bijux_phylogenetics.bayesian import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_branch_length_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    build_substitution_parameter_prior_bundle,
    run_bayesian_inference,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree("dataset/tree.nwk")
records = load_fasta_alignment("dataset/alignment.fasta")

model_definition = build_fixed_topology_dna_model_definition(
    substitution_model_name="HKY85",
    branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
    substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
        kappa_prior=build_exponential_positive_substitution_parameter_prior(
            rate=1.5
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
        ),
    ),
)

proposal_schedule = build_fixed_topology_dna_proposal_schedule(
    model_definition=model_definition,
    branch_length_move_weight=1.0,
    branch_length_log_scale_standard_deviation=0.25,
    kappa_move_weight=1.0,
    kappa_log_scale_standard_deviation=0.35,
    base_frequency_move_weight=1.0,
    base_frequency_coordinate_standard_deviation=0.45,
)

report = run_bayesian_inference(
    tree=tree,
    records=records,
    model_definition=model_definition,
    proposal_schedule=proposal_schedule,
    iteration_count=2_000,
    sample_every=10,
    seed=11,
)

print(report.chain_report.accepted_count)
print(len(report.posterior_rows))
```

The supported fixed-topology DNA substitution families are the public
parameterized models, not `JC69`. Use `K80`, `HKY85`, or `GTR` when building
the model definition.

## Run One Joint-Topology DNA Chain

Joint-topology DNA workflows combine the same sequence-model surface with one
topology prior and explicit topology-move weights.

```python
from bijux_phylogenetics.bayesian import (
    build_exponential_branch_length_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
    build_substitution_parameter_prior_bundle,
    build_uniform_rooted_tree_topology_prior,
    run_bayesian_inference,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree

tree = load_tree("dataset/start-tree.nwk")
records = load_fasta_alignment("dataset/alignment.fasta")

sequence_model_definition = build_fixed_topology_dna_model_definition(
    substitution_model_name="K80",
    branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
    substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
        kappa_prior=build_exponential_positive_substitution_parameter_prior(
            rate=1.5
        )
    ),
)

model_definition = build_joint_topology_dna_model_definition(
    sequence_model_definition=sequence_model_definition,
    topology_prior=build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"]),
)

proposal_schedule = build_joint_topology_dna_proposal_schedule(
    sequence_proposal_schedule=build_fixed_topology_dna_proposal_schedule(
        model_definition=sequence_model_definition,
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.35,
    ),
    nni_move_weight=1.0,
)

report = run_bayesian_inference(
    tree=tree,
    records=records,
    model_definition=model_definition,
    proposal_schedule=proposal_schedule,
    iteration_count=2_000,
    sample_every=10,
    seed=7,
)

print(report.distinct_topology_count)
print(len(report.posterior_rows))
```

Use this surface when topology uncertainty is part of the review question.
Keep the topology prior and topology-move weights explicit, because they are
part of the public contract rather than hidden defaults.

## What The Dispatcher Protects

The dispatcher validates that:

- the model definition is one supported native Bayesian definition
- the proposal schedule is one supported native Bayesian schedule
- fixed-topology definitions only pair with fixed-topology schedules
- joint-topology definitions only pair with joint-topology schedules

That validation is part of the serious public story. Mismatches fail
explicitly instead of producing one ambiguous runtime result.

## How To Read The Run Reports

Treat the returned run report as the stable output contract for the chosen
workflow family.

The most important fields to inspect are:

- `model_definition`
- `proposal_schedule`
- `observation_policy`
- `chain_report`
- `posterior_rows`

Joint-topology runs also keep topology diversity explicit through
`distinct_topology_count` and `distinct_topology_ids`.

## Current Boundaries

- this public page covers the native dispatcher for fixed-topology and
  joint-topology DNA workflows
- this page does not claim that every Bayesian submodule has the same public
  maturity as those two dispatcher-routed contracts
- this page does not rename BEAST or MrBayes wrapper surfaces as native
  inference ownership
- this page does not replace the evidence-book or broader validation studies
  when the question is scientific trust rather than API usage
