---
title: Python Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# Python Surface

The canonical Python import surface lives under `bijux_phylogenetics`.

Use the canonical package name when you need the durable runtime API. The
`phylogenetic` package exists as a compatibility alias, not as a second
independent runtime.

The public promise is ownership clarity: imports should resolve to one runtime
family, not to competing public surfaces with different meanings.
