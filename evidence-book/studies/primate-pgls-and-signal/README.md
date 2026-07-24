# Primate PGLS And Signal

This directory is the governed PCM2 record for baseline regression,
Pagel-lambda PGLS, phylogenetic signal, tree transformations, evolutionary-mode
fits, likelihood-ratio comparisons, residual review, and ancestral estimates
on the prepared primate population.

PCM2 deliberately carries two status layers. Ten focused bundle manifests
adjudicate bounded claims; the aggregate scalar table retains a wider set of
cross-implementation observations. A reader must inspect both before making a
study-wide statement.

## Current Bundle Ledger

| Record | Scientific surface | Verdict |
| --- | --- | --- |
| `evidence-001` | workspace reload semantics | `matched` |
| `evidence-002` | baseline GLS | `matched` |
| `evidence-003` | estimated-lambda regression | `matched_with_tolerance` |
| `evidence-004` | phylogenetic signal and lambda-zero test | `matched_with_tolerance` |
| `evidence-005` | residual diagnostics | `matched_with_tolerance` |
| `evidence-006` | OU, early-burst, and late-burst tree transforms | `matched` |
| `evidence-007` | Brownian, OU, and early-burst fits | `matched_with_tolerance` |
| `evidence-008` | evolutionary-mode likelihood-ratio tests | `matched_with_tolerance` |
| `evidence-009` | named-node ancestral estimates | `matched_with_tolerance` |
| `evidence-010` | `corBlomberg` intercept sweep | `not_comparable` |

The scalar parity table retained with `evidence-001` also contains 12
`mismatch_unexplained` observations across transformed-tree totals,
early-burst fitting, one likelihood-ratio comparison, ancestral vectors, and
one residual diagnostic. Those rows remain scientific debt; favorable focused
bundle verdicts do not remove them.

## Select Evidence From The Claim

| Intended statement | Required focused bundles | Additional record that must remain visible |
| --- | --- | --- |
| baseline GLS correspondence | `evidence-001` population plus `evidence-002` | baseline model identity and diagnostics |
| effect of fitted phylogenetic covariance on the regression | `evidence-002`, `evidence-003`, and `evidence-005` | both coefficient/likelihood records and the unresolved residual observation |
| longevity phylogenetic signal | `evidence-004` | response/tree identity, λ bounds, zero-signal rule, and tolerance |
| Brownian or OU mode comparison | `evidence-006` through `evidence-008` | corresponding aggregate rows and complete candidate denominator |
| early-burst correspondence | `evidence-006` through `evidence-008` | unresolved transform, fitted-mode, and likelihood-ratio scalar rows |
| ancestral correspondence | `evidence-009` | model-specific clade keys plus broader unresolved vector rows |
| complete source-workflow coverage | all ten bundles | `evidence-010` remains `not_comparable` and the 12 mismatches remain in scope |

The table prevents a favorable claim-scoped verdict from being reused for a
broader model family. In particular, evidence for the baseline-to-PGLS slope
change does not establish early-burst transformation or ancestral parity.

## Keep The Two Ledgers Separate

| Ledger | Unit | Status vocabulary | Appropriate use |
| --- | --- | --- | --- |
| focused bundle manifests | one registered claim contract | `matched`, `matched_with_tolerance`, `not_comparable` | decide whether that bounded claim is supportable |
| aggregate scalar parity table | one named cross-implementation observation | row-level outcome and mismatch reason | inspect broader correspondence and unresolved numerical debt |

A bundle can match while a neighboring scalar row remains unexplained because
the two records can have different observation keys, scopes, or comparison
contracts. Quote the focused verdict only with its claim identity; quote a
study-wide mismatch count only with the complete scalar denominator.

## Decision Structure

```mermaid
flowchart LR
    population["prepared primate<br/>population"] --> gls["baseline GLS"]
    population --> covariance["tree transforms and<br/>phylogenetic covariance"]
    covariance --> pgls["Pagel-lambda PGLS"]
    covariance --> modes["BM, OU, early-burst<br/>mode fits"]
    pgls --> diagnostics["residual and signal<br/>diagnostics"]
    modes --> tests["likelihood-ratio and<br/>ancestral comparisons"]
    tests --> open["intercept sweep<br/>not comparable"]
```

The baseline and phylogenetic models share response and taxon identity but do
not share residual covariance. A coefficient is therefore meaningful only
with its model, population, transform, and diagnostic record.

The dependency direction matters: tree transformations feed covariance and
mode fits; fitted models feed likelihood-ratio and ancestral comparisons.
Later matches cannot validate an earlier unexplained transformation row.

## Record Authority

| Record | What it owns |
| --- | --- |
| bundle `manifest.json` | claim identity, comparison mode, freshness basis, limitations, and verdict |
| `inputs.manifest.json` | material input identity and dependency linkage |
| `checks.json` and machine results | observation rules, tolerances, values, and outcomes |
| `evidence-001/results/scalar-parity-table.json` | broad PCM2 scalar comparison, including unexplained mismatches |
| `reference/` | lecture source and registered R reference observations |
| public PCM2 guide | reader route and interpretation boundary |

The generated reviewer summary is a projection of these records. It cannot
turn `not_comparable` into a match or drop mismatch rows from the denominator.

## What Can Be Said

The registered baseline GLS, estimated-lambda regression, signal, and several
focused model comparisons agree directly or within declared tolerances. In
particular, the retained slope change between baseline GLS and Pagel-lambda
PGLS demonstrates the consequence of changing residual covariance on this
population.

A PCM2-wide parity claim is not supported. The intercept sweep lacks a
registered comparable observation, and the scalar ledger retains unexplained
disagreements. These are different failure classes: absence of a comparison
contract is not a numerical mismatch.

For regression claims, retain the response, design columns, factor coding,
admitted taxa, covariance/transform identity, fitted lambda or mode parameters,
likelihood convention, and residual diagnostics. Similar coefficient labels
without this identity are not comparable PCM2 observations.

## Reviewer Route

1. Confirm the governed primate derivatives and source provenance.
2. Choose the focused claim whose model and output match the intended
   statement; read its manifest, inputs, checks, and results.
3. For evolutionary-mode or study-wide language, include the aggregate scalar
   table and its complete mismatch denominator.
4. Inspect `evidence-010` before discussing lecture coverage.
5. Verify freshness after any change to runtime definitions, reference
   observations, derivatives, or tolerance rules.

Continue with the public
[PCM2 study guide](../../../docs/03-bijux-phylogenetics-evidence-book/pcm2-primate-pgls-and-signal/index.md)
for the result tables and citation boundaries.
