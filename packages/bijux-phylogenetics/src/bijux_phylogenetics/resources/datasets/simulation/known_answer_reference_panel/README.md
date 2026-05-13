# Known-Answer Reference Panel

This packaged dataset is a deterministic owned simulation panel built for
recovery review. It stores one true tree, one simulated DNA alignment, one
simulated continuous trait, one simulated discrete trait, and explicit node-
level truth ledgers so downstream workflows can be checked against known
answers instead of only against other inferred outputs.

## Truth Surface

The packaged truth artifacts include:

- `true-tree.nwk`: the rooted simulated tree
- `simulated-alignment.fasta`: the simulated DNA alignment generated on that tree
- `continuous-traits.tsv`: the simulated Brownian tip trait table
- `discrete-traits.tsv`: the simulated discrete tip-state table
- `true-continuous-nodes.tsv`: true continuous values for every node
- `true-discrete-nodes.tsv`: true discrete states for every node
- `true-parameters.tsv`: the exact simulation parameters and seeds

## Simulation Contract

This panel is intentionally deterministic and fully owned by the runtime:

- tree model: birth-death
- alignment model: symmetric JC-like DNA substitution simulation
- continuous trait model: Brownian motion
- discrete trait model: symmetric discrete-state transition simulation

The purpose of the panel is not to claim a realistic biological data-generating
history. Its purpose is to give the runtime one durable known-answer reference
surface where topology recovery, continuous parameter recovery, and ancestral
state recovery can be checked honestly against stored truth.

## Governed Recovery Workflow

The packaged workflow reruns the following reviewer-facing recovery surfaces:

- distance-tree recovery from the simulated alignment
- Brownian trait-evolution fitting on the true tree
- continuous ancestral reconstruction on the true tree
- discrete ancestral reconstruction on the true tree
- explicit truth-versus-estimate recovery ledgers for parameters and internal nodes

## Expected Outputs

The `expected/` directory contains the governed recovery ledgers regenerated
from the owned runtime surface. These files are used as exact packaged
references for the dataset tests and the public demo surface.
