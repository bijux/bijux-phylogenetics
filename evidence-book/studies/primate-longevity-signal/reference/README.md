# PCM1 Runtime Reference Mapping

PCM1 already owns a strong governed parity surface for the data-preparation,
tree-alignment, signal, and ancestral workflow used in the lecture family.

Checked-in reference artifacts:

- `primate_lifespan_signal_reference_r.R`
- `primate_lifespan_signal_reference_bijux.py`
- `lund_pcm1_plots_signal_source.R`

Checked-in lecture source copies now live beside the study:

- `../datasets/primate_raw.xlsx`
- `../datasets/primate.RData`
- `../datasets/primate.csv`
- `../datasets/primatetree.nex`
- `../datasets/trimmed_primatetree.nex`

Relevant owned surfaces include:

- workbook loading and type-repair utilities used to produce governed primate
  derivatives
- tree import, pruning, and correspondence checks
- phylogenetic signal estimation and ancestral-state reporting

What remains outside the checked-in parity claim:

- lecture plotting choices that depend on presentation-specific R graphics
  surfaces
- plot-for-plot identity with the lecture visualization stack

This directory exists so reviewers can distinguish the strong numerical PCM1
parity record from the intentionally narrower plotting boundary.
