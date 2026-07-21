---
title: phytools Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# phytools Correspondence

phytools is the live reference for phylogenetic signal, Mk fitting,
continuous-trait utilities, stochastic histories, and ancestral-state
summaries. The governed registry currently contains 65 cases across 21 shards.
Deterministic estimates and stochastic summaries use different acceptance
rules.

## Governed Surface

| Territory | phytools references | Correspondence level |
| --- | --- | --- |
| signal | `phylosig` with K and lambda | statistics, fitted lambda and named likelihood quantities |
| discrete models | `fitMk` with ER, SYM and ARD | rate rows, likelihood and supported state ordering |
| ancestral states | `fastAnc`, `anc.ML`, `rerootingMethod` | node estimates, uncertainty or state probabilities under comparable models |
| regression and residuals | `pgls.SEy`, `phyl.resid` with BM and lambda | bounded coefficient, covariance and residual summaries |
| continuous simulation | `fastBM`, `sim.corrs` | moment, variance, covariance and correlation summaries |
| stochastic histories | `make.simmap`, `countSimmap`, `describe.simmap`, `densityMap`, `sim.history` | transition-count, dwell-time, state-density and history-summary envelopes |
| comparative tests | `phylANOVA` | declared group and test summaries under governed fixtures |

## Deterministic And Stochastic Evidence

For deterministic cases, the registry can compare named scalars, vectors,
rates, node values, or probabilities directly within a tolerance. For
simulation and stochastic mapping, exact histories are not expected to match
across random-number generators. Those cases compare declared distributional
or aggregate summaries across governed replication policies.

`pgls.SEy` is intentionally narrow. It validates the comparable
measurement-error PGLS surface; it is not presented as a general
`phytools::pgls` correspondence claim. Broader coefficient comparisons use the
checked APE plus `nlme` reference where that is the actual source.

## What A Passing Stochastic Case Means

A stochastic case passes when its registered replication policy completes and
every declared aggregate or distribution rule holds. It does not require the R
and Bijux histories to contain the same random events. The record therefore
retains replicate count, seed policy, summary definition, tolerance or envelope,
failed replicates, and observed distributional quantities.

Do not remove non-finite, degenerate, or failed replicates before computing a
headline summary unless the registry declares that exclusion. Selected,
executed, passed, failed, and skipped case counts remain separate from
within-case replication denominators.

## Run The Live Lane

```bash
RUN_LIVE_PARITY=1 bijux-phylogenetics parity \
  --reference-source phytools-live \
  --phytools-shard phytools-fitmk-model-er \
  --summary-out artifacts/parity/phytools-summary.tsv \
  --observations-out artifacts/parity/phytools-observations.tsv \
  --phytools-failure-root artifacts/parity/phytools-failures \
  --json
```

Use `--phytools-case` for an individual fixture. The live lane requires
`Rscript`, phytools, and `RUN_LIVE_PARITY=1`.

## Programmatic Surface

| Public locator | Purpose |
| --- | --- |
| `bijux_phylogenetics.parity.phytools:list_phytools_parity_cases` | enumerate the 65 deterministic and stochastic contracts |
| `bijux_phylogenetics.parity.phytools:list_phytools_parity_shards` | enumerate 21 operation-owned shards |
| `bijux_phylogenetics.parity.phytools:run_phytools_parity_cases` | execute selected live cases under the opt-in boundary |
| `bijux_phylogenetics.parity.phytools:write_phytools_parity_summary_table` | persist case and shard status totals |
| `bijux_phylogenetics.parity.phytools:write_phytools_parity_observation_table` | persist scalar, row, distribution, and mismatch observations |

Discovery defines the eligible case population before random draws occur. A
stochastic case cannot be reselected after viewing its summaries.

## Preserve Both Denominators

A stochastic parity report has two denominators:

1. selected registry cases, partitioned into executed, passed, failed, and
   skipped cases;
2. within each executed case, requested, completed, invalid, excluded, and
   summarized replicates.

A case can execute yet fail because too few valid replicates remain. It can
pass a distribution rule while still retaining individual invalid draws. Do
not substitute the replicate pass fraction for the registry case result, or
the case pass fraction for the within-case sampling record.

## Boundaries That Remain Visible

- ER and SYM rerooting cases are compared only where root-prior semantics can
  be aligned; neighboring models remain native review surfaces without a live
  claim.
- Stochastic-history agreement concerns declared summaries, never exact event
  placement on every branch.
- Density maps and other plotting-adjacent results are compared through their
  underlying numerical ledgers, not pixel identity.
- A passed shard validates its registered fixtures and versioned conventions;
  it does not certify all phytools arguments or future package behavior.
- A passed harness case is not a `matched` Evidence Book verdict; scientific
  promotion requires a governed claim, provenance, checks, and freshness.
