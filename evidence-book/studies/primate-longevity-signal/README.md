# Primate Longevity Signal

This study holds reviewer-facing parity evidence for the Lund `PCM1_plots_signal`
teaching workflow and the corresponding `bijux-phylogenetics` reproductions.

Study metadata:

- study id: `primate-longevity-signal`
- current evidence id: `evidence-001`
- current reviewer status: `analytical_alignment_complete`
- evidence rule: report reproduction, deviation, and unclaimed scope explicitly

The checked-in reference sources for this study live here:

- [R reference checks](./reference/primate_lifespan_signal_reference_r.R)
- [Python `bijux-phylogenetics` checks](./reference/primate_lifespan_signal_reference_bijux.py)
- [Evidence builder](./build_evidence.py)

The governed study surfaces live here:

- [Study manifest](./study.json)
- [Current evidence bundle](./evidence-001/README.md)

Current outcome at a glance:

- analytical and structural reproduction is in place for the registered example, with all executable script lines tracked
- block verdicts are explicit: `verified`, `verified_with_tolerance`, `plot_only`, `artifact_only`, `seeded_input_only`, or `workflow_only`
- rendered figure surfaces are tracked honestly as `plot_only`; they are not falsely advertised as identical to the R plotting stack
- reviewer-facing bundle extras now include a checklist, reproducibility manifest, maturity registry, and scientific debt register
