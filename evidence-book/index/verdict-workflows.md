# Verdict Workflows

This index explains how `mismatch_explained`, `mismatch_unexplained`, and
`not_comparable` evidence states are supposed to stay visible.

## mismatch_explained

Keep the row visible, record the explanation kind explicitly, and do not promote it to a full match unless the reference source becomes more precise.

- entries: `4`

- `primate-pgls-and-signal-evidence-001-baseline-intercept`
  Path: `studies/primate-pgls-and-signal/evidence-001`
  Explanation: The R reference stores the baseline intercept rounded to four decimal places; the Bijux value rounds back to the same published scalar.
- `primate-pgls-and-signal-evidence-001-baseline-log-likelihood`
  Path: `studies/primate-pgls-and-signal/evidence-001`
  Explanation: The R reference stores the baseline log likelihood rounded to four decimal places; the Bijux value rounds back to the same published scalar.
- `primate-pgls-and-signal-evidence-001-baseline-r-squared`
  Path: `studies/primate-pgls-and-signal/evidence-001`
  Explanation: The R reference stores the baseline R-squared rounded to four decimal places; the Bijux value rounds back to the same published scalar.
- `primate-pgls-and-signal-evidence-001-baseline-slope`
  Path: `studies/primate-pgls-and-signal/evidence-001`
  Explanation: The R reference stores the baseline slope rounded to four decimal places; the Bijux value rounds back to the same published scalar.

## mismatch_unexplained

Keep the row visible, treat it as open scientific debt, and require explicit closure rather than silent tolerance inflation.

- entries: `1`

- `primate-pgls-and-signal-evidence-001-estimated-diagnostic-fitted-correlation`
  Path: `studies/primate-pgls-and-signal/evidence-001`

## not_comparable

Keep the boundary explicit, point to the governing claim and evidence bundle, and do not restate it as a pass/fail parity surface until the runtime owns the comparison.

- entries: `4`

- `primate-longevity-signal-evidence-001-pcm1-artifact-provenance-tracking`
  Path: `studies/primate-longevity-signal/evidence-001`
  Claim: PCM1 saved artifacts remain governed provenance surfaces
- `primate-longevity-signal-evidence-001-pcm1-reproducibility-contract-tracked`
  Path: `studies/primate-longevity-signal/evidence-001`
  Claim: PCM1 reproducibility contract is tracked explicitly
- `primate-longevity-signal-evidence-001-pcm1-visual-surface-tracking`
  Path: `studies/primate-longevity-signal/evidence-001`
  Claim: PCM1 visual surfaces are tracked without false equivalence claims
- `primate-pgls-and-signal-evidence-006-pcm2-coverage-boundary-explicit`
  Path: `studies/primate-pgls-and-signal/evidence-006`
  Claim: PCM2 uncovered evolutionary-mode and ancestral fragments remain explicit
