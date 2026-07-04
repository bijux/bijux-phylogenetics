# BEAST Posterior Reference Bundle

This bundle records the governed reference summaries derived from the checked
BEAST XML, posterior log, and posterior tree-set fixtures.

It is the reviewer-facing evidence surface for BEAST posterior summarization
when tests need stable checked outputs rather than a live engine run.

- `consensus-10pct*` files capture the consensus summary after discarding the
  first 10% of posterior samples as burn-in.
- `mcc-10pct*` files capture the maximum-clade-credibility summary under the
  same burn-in policy.
