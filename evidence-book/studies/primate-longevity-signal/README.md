# Primate Longevity Signal

This study holds reviewer-facing parity evidence for the Lund `PCM1_plots_signal`
teaching workflow and the corresponding `bijux-phylogenetics` reproductions.

Study metadata:

- study id: `primate-longevity-signal`
- summary evidence id: `evidence-001`
- component evidence ids: `evidence-002` through `evidence-009`
- current reviewer status: `analytical_alignment_complete`
- evidence rule: report reproduction, deviation, and unclaimed scope explicitly

The checked-in reference sources for this study live here:

- [R reference checks](./reference/primate_lifespan_signal_reference_r.R)
- [Python `bijux-phylogenetics` checks](./reference/primate_lifespan_signal_reference_bijux.py)

The governed study surfaces live here:

- [Summary evidence bundle](./evidence-001/README.md)
- [Shared datasets](./datasets/)
- [Study provenance](./provenance/)

Current outcome at a glance:

- analytical and structural reproduction is in place for the registered example, with all executable script lines tracked
- preprocessing, tree import, diagnostics, correspondence, and processed-export contracts now have dedicated Evidence IDs instead of living only inside one summary bundle
- block verdicts are explicit: `verified`, `verified_with_tolerance`, `plot_only`, `artifact_only`, `seeded_input_only`, or `workflow_only`
- rendered figure surfaces are tracked honestly as `plot_only`; they are not falsely advertised as identical to the R plotting stack
- reviewer-facing bundle extras now include a checklist, reproducibility manifest, maturity registry, and scientific debt register
