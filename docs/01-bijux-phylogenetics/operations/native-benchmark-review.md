---
title: Native Benchmark Review
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Native Benchmark Review

Use the native benchmark family when the question is not "can I run one owned
inference result" but "what is the current public benchmark posture of the
native maximum-likelihood surface".

The main review lane is the unified suite:

```bash
bijux-phylogenetics benchmark native-maximum-likelihood-suite --json
```

That suite is the public summary surface because it keeps scaling, memory,
truth-recovery, and wrapper-comparison evidence together without flattening
the outcome into one generic pass state.

## Start With The Unified Suite

The suite report keeps these governed outcome classes explicit:

- `native-advantage`
- `native-bug`
- `expected-assumption-difference`
- `unsupported-case`

Read those classes literally:

- `native-advantage` means the owned native surface did something better on the
  governed case
- `native-bug` means the discrepancy is treated as a real owned defect
- `expected-assumption-difference` means the mismatch is explained by a stated
  model or engine boundary
- `unsupported-case` means the case should not be misread as covered

That distinction is what makes the benchmark family release-facing evidence
rather than one synthetic score.

## Use Narrower Lanes For Focused Questions

Use the narrower commands when the review question is specific:

- `benchmark native-maximum-likelihood-speed`
- `benchmark native-maximum-likelihood-memory`
- `benchmark native-maximum-likelihood-accuracy`
- `benchmark maximum-likelihood-wrapper-comparison`

Examples:

```bash
bijux-phylogenetics benchmark native-maximum-likelihood-speed --json
bijux-phylogenetics benchmark native-maximum-likelihood-memory --json
bijux-phylogenetics benchmark native-maximum-likelihood-accuracy --json
bijux-phylogenetics benchmark maximum-likelihood-wrapper-comparison --json
```

Use the suite first for release posture. Use the narrower lanes when one
performance or correspondence question needs its own artifact trail.

## Override Governed Inputs Deliberately

The suite exposes bounded override flags so one reviewer can stress a specific
native lane without editing code:

```bash
bijux-phylogenetics benchmark native-maximum-likelihood-suite \
  --scaling-model jc69 \
  --scaling-search-method nni \
  --accuracy-model k80 \
  --accuracy-search-method nni \
  --size-class medium \
  --taxon-count 16 \
  --site-count 512 \
  --start-tree-count 1 \
  --start-tree-seed 17 \
  --upper-branch-length-bound 1.0 \
  --json
```

Those flags are for bounded review. They do not change the repository's
governed release story by themselves.

## What To Inspect In The JSON Output

At minimum, inspect:

- `metrics.observation_count`
- `metrics.summary_row_count`
- `metrics.native_advantage_case_count`
- `metrics.native_bug_case_count`
- `metrics.expected_assumption_difference_case_count`
- `metrics.unsupported_case_count`
- `data.summary_rows`
- `data.observations`
- `data.limitations`

Use `summary_rows` to see the release-facing status breakdown. Use
`observations` to inspect the actual case ids, dimensions, comparison
surfaces, and notes behind that summary.

## How This Relates To The Rest Of The Runtime Story

- use [native maximum-likelihood workflows](native-maximum-likelihood-workflows.md)
  when the main output is one owned tree-inference result object
- use this benchmark workflow when the main output is one review of native
  speed, memory, recovery, and wrapper-comparison posture
- use [native benchmark evidence](../quality/native-benchmark-evidence.md) when
  the question is how to interpret the benchmark family as one public trust
  surface

## Current Boundaries

- this workflow documents the public benchmark family for native
  maximum-likelihood inference
- this page does not claim that one benchmark run closes wider scientific
  validity questions
- this page does not relabel wrapper-comparison evidence as native-only proof
- this page does not hide bugs or unsupported cases behind aggregate success
  language
