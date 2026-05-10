# Primate missing-data accounting parity bundle

Governed evidence for how missing values are handled before the processed primate table is used downstream.

- evidence id: `evidence-004`
- claim id: `pcm1-missing-data-accounting-parity`
- verdict: `matched`
- source fragments: `primate-data-preprocessing`

Governed files:

- `missing-data-accounting-parity.json`
- `claims.json`
- `manifest.json`

Reference script locators:

- `evidence-book/studies/primate-longevity-signal/reference/primate_lifespan_signal_reference_r.R#L41-L47`

Limitations:

- The raw workbook remains external, so the bundle records the governed missing-value contract at the processed artifact boundary.
