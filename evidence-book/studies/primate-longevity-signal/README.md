# Primate Longevity Signal

This study turns the Lund `PCM1_plots_signal` lecture into a governed set of
small evidence slices covering primate data intake, tree preparation,
phylogenetic signal, and ancestral reconstruction.

It is the introductory comparative study in the evidence book, but it is not
lightweight. The repository already carries the full evidence breakdown rather
than one monolithic summary claim.

## Study Metadata

- study id: `primate-longevity-signal`
- lecture family: `PCM1_plots_signal`
- deployment state: `governed_parity_deployed_with_plotting_boundary`
- current runtime posture:
  `primate_signal_and_ancestral_workflow_governed_with_explicit_plot_boundary`
- summary evidence id: `evidence-001`
- component evidence ids: `evidence-002` through `evidence-009`
- current reviewer status: `analytical_alignment_complete`
- evidence rule: report reproduction, deviation, and unclaimed scope explicitly

## Study Story

PCM1 starts with the practical foundations that later lectures depend on:

- reading and normalizing the primate workbook
- importing and pruning the primate tree
- checking species correspondence between the tree and the trait table
- estimating phylogenetic signal on the cleaned primate response
- projecting ancestral values back onto the tree

Because those steps are easy to collapse into one invisible preprocessing
block, the study keeps them separated as explicit evidence IDs.

## Evidence Map

| Evidence ID | Focus | Current repository state |
| --- | --- | --- |
| `primate-longevity-signal/evidence-001` | End-to-end summary bundle across the full lecture workflow | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-002` | Workbook-loading parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-003` | Type-repair parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-004` | Missing-data accounting parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-005` | Duplicate-species aggregation parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-006` | Tree import parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-007` | Tree diagnostics parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-008` | Tree-data correspondence parity | Governed parity bundle deployed |
| `primate-longevity-signal/evidence-009` | Processed-export parity | Governed parity bundle deployed |

## Ideal Comparative Checks

- Keep data normalization separate from tree operations so reviewers can see
  exactly where the comparative dataset became trustworthy.
- Show which outputs are numerical parity claims and which are only workflow or
  plotting surfaces.
- Preserve the route from raw workbook and raw tree to the governed primate
  derivatives reused by later studies.

## Elite Comparative Checks

- Keep reviewer checklists, maturity registries, scientific debt registers,
  and reproducibility manifests attached to the main summary bundle.
- Track every lecture block that matters to later studies so PCM2 and PCM3 do
  not need to rediscover baseline primate assumptions.

## Study Surfaces

- [Summary evidence bundle](./evidence-001/README.md)
- [Shared datasets](./datasets/)
- [Study provenance](./provenance/)
- [runtime reference mapping](./reference/README.md)
- checked-in lecture source copies:
  `primate_raw.xlsx`, `primate.RData`, `primate.csv`, `primatetree.nex`, and
  `trimmed_primatetree.nex`
