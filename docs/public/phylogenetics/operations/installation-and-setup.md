---
title: Installation and Setup
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-13
---

# Installation And Setup

Use the canonical runtime package when you want the main product surface.

## Basic Route

- install the runtime package
- confirm the CLI resolves
- run a small tree inspection or validation workflow
- direct local run products into `artifacts/`

The operational goal is predictable local behavior, not hidden notebook state
or scattered temporary outputs.

## External Engine Runtime

Some owned workflows execute scientific binaries directly instead of using
mocked adapters. Those lanes need the local runtime to be installed before the
real integration or governed scientific-validation commands will pass.

On macOS with Homebrew, use:

```bash
make install-external-engine-runtime
```

That target installs the repository's current local engine set:

- `mafft`
- `trimal`
- `iqtree2`
- `FastTree`
- `mb` from MrBayes
- `beast` from BEAST2

If you install those executables by another route, make sure they resolve on
`PATH` before you run the real local lanes.

## Local Verification Routes

Use the split verification commands when you want a real executable check
instead of the default lightweight test surface.

```bash
make test-external-engines
make test-scientific-validation-slow
```

`make test-external-engines` runs the compact real-engine integration lane plus
the non-slow governed scientific-validation reruns.
`make test-scientific-validation-slow` keeps the slower FASTA-to-tree golden
reruns separate so they can use a longer timeout budget without slowing every
local engine check.

Keep the outputs from those runs under the repository `artifacts/` tree. The
owned make surfaces already do that by default.
