---
title: PCM2 Diagnostics And Open Boundaries
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# PCM2 Diagnostics And Open Boundaries

PCM2 retains diagnostic, ancestral, and incomplete coverage as separate
evidence units. Strong agreement in the main regression path therefore does
not erase an unresolved comparison elsewhere in the lecture.

## Evidence State

| Bundle | Surface | Verdict state |
| --- | --- | --- |
| `evidence-004` | phylogenetic signal | governed parity deployed |
| `evidence-005` | residual diagnostics | governed parity deployed |
| `evidence-009` | Brownian and early-burst ancestral reconstruction | `matched_with_tolerance` |
| `evidence-010` | intercept-mode likelihood sweep | `not_comparable`; boundary register |

## Residual Diagnostics

Residual review qualifies coefficient and model-fit claims. Inspect residual
shape, fitted-value relationships, influential taxa, variance behavior, and
any remaining phylogenetic pattern under the same row exclusions and fitted
covariance. Diagnostic agreement is not merely plot resemblance; it must be
tied to comparable quantities and observations.

## Signal

Signal estimates depend on response, taxa, tree, branch lengths, transformation
bounds, and optimization. A signal result should be interpreted with the model
comparison and residual record, not as an isolated certification that
phylogenetic correction succeeded.

## Ancestral Modes

`evidence-009` compares node-level ancestral estimates for Brownian and
early-burst modes within governed tolerance. Rendered node-size equivalence is
not claimed. Stable node identity, rather than display order, anchors the
comparison.

## Intercept-Sweep Boundary

`evidence-010` records the remaining `corBlomberg` intercept-mode likelihood
sweep. Its verdict is intentionally `not_comparable`: the bundle identifies
the uncovered lecture fragment and closure location but does not manufacture a
parity result.

Closing the boundary requires a comparable implementation or reference path,
matching parameter and intercept conventions, recorded likelihood values,
explicit tolerance, and adjudicated checks. Rewording the claim or borrowing a
neighboring PGLS verdict would not close it.

## Supported Conclusion

PCM2 has governed agreement for baseline and λ PGLS, signal, diagnostics,
transformed trees, mode fits, likelihood-ratio comparison, and bounded
ancestral outputs. The intercept-mode sweep remains outside numerical parity.
Both facts must accompany any summary of study coverage.
