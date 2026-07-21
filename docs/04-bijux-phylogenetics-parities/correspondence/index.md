---
title: Interpreting Correspondence
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Interpreting Correspondence

Correspondence is a claim about named outputs under named assumptions. It is
not a promise that two libraries expose identical APIs, optimize identical
objectives, or will agree for every dataset.

## Comparison Record

| Field | Reviewer question |
| --- | --- |
| reference source | Did the run use a live tool or a checked fixture? |
| case and shard | Which operation, model, and input family were selected? |
| versions | Which Bijux, Python, R, package, or executable versions ran? |
| normalization | How were taxa, states, units, roots, likelihood signs, and support scales aligned? |
| observations | Which exact scalar, vector, matrix, topology, or distribution summaries were compared? |
| tolerance | What absolute, relative, structural, or distribution rule decides agreement? |
| case status | Did the registered run pass, fail, or skip? |
| mismatch reason | Which execution, normalization, field, row, tolerance, or availability condition explains a non-pass? |
| artifacts | Where are inputs, raw outputs, observations, warnings, and failure details retained? |

## Runtime Case Status

- `passed`: the selected case executed and every registered observation rule
  held, including declared tolerances where applicable;
- `failed`: Bijux or reference execution failed, required output was missing,
  or at least one registered comparison rule was violated;
- `skipped`: the live reference environment was unavailable or the case was
  not executable under the selected lane.

Retain separate selected, passed, failed, and skipped counts. “Passed” is a
runtime harness outcome, not the `matched` verdict stored in a scientific
evidence bundle.

## Promotion To An Evidence Verdict

An Evidence Book bundle may interpret retained parity observations only after
it binds them to a claim, source, governed inputs, explicit checks, provenance,
and freshness. Its allowed verdicts are `matched`, `matched_with_tolerance`,
`mismatch_explained`, `mismatch_unexplained`, and `not_comparable`.

Promotion is not a spelling conversion. A failed case needs scientific review
before it can become an explained or unexplained mismatch. A skipped case may
leave evidence unavailable, but does not automatically establish
`not_comparable`. Resolution state—such as open follow-through or a named
blocker—remains a separate axis.

## Reproducing A Comparison

The checked-fixture lane runs without R:

```bash
bijux-phylogenetics parity --reference-source checked-fixture --json
```

Live R lanes are explicit and can be narrowed to registered cases or shards:

```bash
RUN_LIVE_PARITY=1 bijux-phylogenetics parity \
  --reference-source ape-live \
  --ape-shard load_tree \
  --summary-out artifacts/parity/ape-summary.tsv \
  --observations-out artifacts/parity/ape-observations.tsv \
  --json
```

Use the tool-specific registry rather than assuming an R function name has a
Bijux equivalent. Unsupported operations should fail with a named boundary,
not fall through to a neighboring approximation.

## Denominator Check

Before quoting correspondence, record how many cases were registered,
selected, executed, passed, failed, and skipped. Then identify the exact
shards. A statement about “all cases” is invalid when it silently means only
the subset that executed successfully.

## Review Resources

- [Verification matrix](verification-matrix.md) identifies the strongest
  current claim for each external reference.
- [Artifacts and reproduction](artifacts-and-reproduction.md) explains how to
  retain selections, observations, manifests, and failures.
