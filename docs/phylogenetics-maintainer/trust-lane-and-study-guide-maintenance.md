---
title: Trust Lane And Study Guide Maintenance
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Trust Lane And Study Guide Maintenance

The public trust story now depends on more than general docs and release
guides. It also depends on:

- quality pages that explain validation lanes and how to read them
- foundation pages that explain how to interpret runtime, native, wrapper,
  evidence, and release claims
- PCM study guides that explain what the public evidence studies do and do not
  prove

## Review These Together

When trust language or study coverage changes, review:

- `docs/phylogenetics/quality/`
- `docs/phylogenetics/foundation/`
- `docs/phylogenetics-evidence-book/`
- root and package README pointers into those sections

## Questions To Ask

- does the public quality section still explain what each validation lane
  proves and does not prove
- does the foundation section still separate runtime, native, wrapper,
  evidence, and release claims clearly
- do the PCM study guides still describe their analytical sequence and bounded
  scope honestly
- do entry-point pages still route readers toward the strongest trust-reading
  surfaces

## Failure Pattern

The failure here is trust compression:

- many different proof surfaces exist
- the docs stop distinguishing them clearly
- readers start overreading one lane or one study as broader closure than it is

That is the exact pattern these guides are meant to prevent.
