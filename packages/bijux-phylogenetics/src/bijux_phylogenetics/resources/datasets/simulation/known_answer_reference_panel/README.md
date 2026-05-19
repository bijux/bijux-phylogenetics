# Known-Answer Reference Panel

This packaged dataset is a deterministic owned simulation panel built for
recovery review. It stores one true tree, one simulated DNA alignment, one
simulated Brownian continuous trait, one simulated OU continuous trait, one
generic discrete trait, one host-state trait, one geographic-state trait, and
explicit node- and branch-level truth ledgers so downstream workflows can be
checked against known answers instead of only against other inferred outputs.

## Truth Surface

The packaged truth artifacts include:

- `true-tree.nwk`: the rooted simulated tree
- `simulated-alignment.fasta`: the simulated DNA alignment generated on that tree
- `continuous-traits.tsv`: the simulated Brownian tip trait table
- `ou-traits.tsv`: the simulated OU tip trait table
- `discrete-traits.tsv`: the simulated generic discrete tip-state table
- `host-traits.tsv`: the simulated host-state tip table
- `geographic-traits.tsv`: the simulated geographic-state tip table
- `true-continuous-nodes.tsv`: true continuous values for every node
- `true-ou-nodes.tsv`: true OU continuous values for every node
- `true-discrete-nodes.tsv`: true discrete states for every node
- `true-host-nodes.tsv`: true host states for every node
- `true-geographic-nodes.tsv`: true geographic states for every node
- `true-host-switch-events.tsv`: true branchwise simulated host-state event ledger
- `true-geographic-transition-events.tsv`: true branchwise simulated geographic event ledger
- `true-parameters.tsv`: the exact simulation parameters and seeds
- `recovery-thresholds.tsv`: explicit pass or fail thresholds for the governed recovery bundle

## Simulation Contract

This panel is intentionally deterministic and fully owned by the runtime:

- tree model: birth-death
- alignment model: symmetric JC-like DNA substitution simulation
- continuous trait models: Brownian motion and Ornstein-Uhlenbeck
- discrete trait models: symmetric discrete-state transition simulation for
  generic state, host-state, and geographic-state surfaces

The purpose of the panel is not to claim a realistic biological data-generating
history. Its purpose is to give the runtime one durable known-answer reference
surface where topology recovery, continuous parameter recovery, ancestral-state
recovery, host-switch recovery, and geographic-transition recovery can be
checked honestly against stored truth.

## Governed Recovery Workflow

The packaged workflow reruns the following reviewer-facing recovery surfaces:

- distance-tree recovery from the simulated alignment
- Brownian and OU trait-evolution fitting on the true tree
- continuous ancestral reconstruction on the true tree
- generic discrete ancestral reconstruction on the true tree
- host-switch reconstruction against stored host-state and branch-event truth
- geographic ancestral-state and transition reconstruction against stored geographic truth
- explicit truth-versus-estimate recovery ledgers for parameters, internal nodes, branch changes, and threshold checks

## Expected Outputs

The `expected/` directory contains the governed recovery ledgers regenerated
from the owned runtime surface. These files are used as exact packaged
references for the dataset tests and the public demo surface.

The governed expected bundle contains:

- `workflow-summary.tsv`: one integrated recovery metrics row for the panel
- `recovered-distance-tree.nwk`: the tree rebuilt from the simulated alignment
- `tree-recovery.tsv`: topology recovery metrics against the true tree
- `parameter-recovery.tsv`: continuous truth-versus-estimate parameter ledger
- `brownian-fit-summary.tsv` and `ou-fit-summary.tsv`: continuous model-fit summaries on the true tree
- `continuous-ancestral-summary.tsv` and `continuous-ancestral-uncertainty.tsv`: continuous ancestral reconstruction summaries
- `continuous-node-recovery.tsv`: continuous internal-node truth-versus-estimate ledger
- `discrete-ancestral-summary.tsv` and `discrete-ancestral-probabilities.tsv`: generic discrete ancestral summaries
- `discrete-node-recovery.tsv`: generic discrete internal-node recovery ledger
- `host-switch-summary.tsv`, `host-state-nodes.tsv`, and `host-switch-branches.tsv`: host-state recovery and branch-change review outputs
- `host-node-recovery.tsv` and `host-event-recovery.tsv`: host-state truth-versus-estimate ledgers
- `geographic-ancestral-summary.tsv`, `geographic-state-probabilities.tsv`, and `geographic-transition-summary.tsv`: geographic-state recovery outputs
- `geographic-node-recovery.tsv` and `geographic-event-recovery.tsv`: geographic truth-versus-estimate ledgers
- `recovery-threshold-evaluation.tsv`: declared threshold evaluation ledger for the full recovery suite
