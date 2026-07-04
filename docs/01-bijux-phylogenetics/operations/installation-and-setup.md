---
title: Installation and Setup
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# Installation And Setup

Use the canonical runtime package when you want the main product surface.

## Basic Route

1. Install the runtime package.
2. Confirm the CLI resolves.
3. Run a small validation or inspection command.
4. Keep local run products under `artifacts/`.

The operational goal is predictable local behavior, not hidden notebook state
or scattered temporary outputs.

## External Engine Runtime

Some workflows execute real scientific binaries rather than mock adapters. On
macOS with Homebrew, install the current local engine set with:

```bash
make install-external-engine-runtime
```

That target installs:

- `mafft`
- `trimal`
- `iqtree2`
- `FastTree`
- `mb` from MrBayes
- `beast` from BEAST2

If you install those executables another way, make sure they resolve on
`PATH` before you run the real local lanes.

## Local Verification Routes

Use these commands when you want an executable-backed verification rather than
the default lightweight surface:

```bash
make test-external-engines
make test-scientific-validation-slow
```

- `make test-external-engines`
  Runs the compact real-engine integration lane plus the non-slow governed
  scientific-validation reruns.
- `make test-scientific-validation-slow`
  Keeps the longer FASTA-to-tree reruns separate so they do not slow every
  local engine check.

The repository-owned make targets already route outputs into `artifacts/` by
default.
