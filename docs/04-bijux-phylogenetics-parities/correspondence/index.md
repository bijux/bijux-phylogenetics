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
| verdict | Did the case match, mismatch, skip, or remain not comparable? |
| artifacts | Where are inputs, raw outputs, observations, warnings, and failure details retained? |

## Verdict Semantics

- **matched**: every declared observation satisfies its acceptance rule.
- **matched with tolerance**: bounded numerical differences satisfy a stated
  tolerance; this is not exact identity.
- **mismatched**: at least one comparable observation violates the rule.
- **not comparable**: assumptions, outputs, or ownership do not support the
  proposed comparison.
- **skipped**: the live reference could not run or the governed case was not
  selected; a skip is never counted as a pass.

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

## Review Resources

- [Verification matrix](verification-matrix.md) identifies the strongest
  current claim for each external reference.
- [Artifacts and reproduction](artifacts-and-reproduction.md) explains how to
  retain selections, observations, manifests, and failures.
