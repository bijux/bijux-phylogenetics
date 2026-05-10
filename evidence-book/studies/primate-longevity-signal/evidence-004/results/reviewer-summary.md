# Primate missing-data accounting parity bundle Reviewer Summary

- study: `primate-longevity-signal`
- evidence id: `evidence-004`
- comparison mode: `direct_parity`
- verdict: `matched`
- claims: `1`

## Summary

- Governed evidence for how missing values are handled before the processed primate table is used downstream.
- Comparison mode: direct_parity. Bundle verdict: matched.
- Claims stay explicit and reviewer-readable through governed claim rows, portable source locators, and tracked limitations.

## Claims

- PCM1 missing-data accounting is visible before downstream inference

## Source Fragments

- `primate-data-preprocessing`

## Limitations

- The raw workbook remains external, so the bundle records the governed missing-value contract at the processed artifact boundary.
