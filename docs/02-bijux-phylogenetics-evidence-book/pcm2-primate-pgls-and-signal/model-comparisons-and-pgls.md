---
title: PCM2 Model Comparisons And PGLS
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# PCM2 Model Comparisons And PGLS

PCM2 compares ordinary residual independence with phylogenetically structured
covariance and evaluates continuous evolutionary modes on the same prepared
primate data. Shared response, predictors, taxa, and transformation conventions
make changes attributable to the covariance/model rather than data drift.

```mermaid
flowchart TB
    data["Governed primate data"]
    gls["Baseline GLS"]
    lambda["Pagel-λ PGLS"]
    transforms["Transformed trees"]
    modes["Brownian · OU<br/>early burst"]
    compare["Coefficients · likelihoods<br/>diagnostics · LRT"]

    data --> gls --> compare
    data --> lambda --> compare
    data --> transforms --> modes --> compare
```

## Baseline And Phylogenetic Regression

Baseline GLS provides the non-phylogenetic residual model. PGLS changes the
covariance matrix according to tree structure and its fitted transformation.
Coefficient differences can arise from changed weighting and correlation even
when the formula and observations are unchanged.

The governed bundles record matched baseline GLS (`evidence-002`) and
Pagel-λ regression (`evidence-003`) claims.

## Tree Transformations And Evolutionary Modes

`evidence-006` covers transformed-tree behavior. `evidence-007` covers the
continuous mode fits. Brownian, OU, and early-burst families encode different
expectations about trait covariance or rate through time; their parameters and
likelihoods are meaningful only under matching conventions.

## Likelihood-Ratio Comparison

`evidence-008` records the likelihood-ratio parity surface. A valid comparison
requires matching data, likelihood definition, optimization treatment, and a
defensible nesting and degrees-of-freedom rule. Information criteria and
likelihood-ratio tests answer related but different selection questions.

## Comparison Ledger

For each fitted model retain:

- response transformation, predictor formula, and included taxa;
- covariance or tree transformation and parameter bounds;
- coefficient estimates and uncertainty;
- log likelihood, parameter count, and optimization status;
- residual diagnostics;
- comparison statistic, degrees of freedom, and decision rule.

This ledger prevents a change in taxa or transformation from being mistaken
for a model effect.

## Interpretation Limits

- PGLS accounts for the declared covariance; it does not remove confounding or
  establish causality.
- Better fit does not prove the evolutionary process is literally true.
- Boundary parameters can invalidate naive asymptotic test interpretations.
- Parity establishes agreement under the study contract, not universal model
  recovery.

Continue with
[signal diagnostics and open boundaries](signal-diagnostics-and-open-boundaries.md).
