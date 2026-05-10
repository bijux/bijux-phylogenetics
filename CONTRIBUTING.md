# Contributing

This repository publishes one runtime package and one maintainer package, but
it should still feel like one coherent system. Good contributions improve one
surface without weakening the shared contracts that keep the workspace aligned.

## Start With Placement

Choose the narrowest correct home for every change:

- package runtime behavior belongs in the owning package under `packages/`
- repository-wide automation belongs in `makes/`
- shared tool policy belongs in `configs/`
- checked-in API contracts, when they exist, belong in versioned package
  directories under `apis/`
- repository handbook and maintainer guidance belong at the root or in `docs/`

If the same rule would need to be copied into multiple packages, that is a sign
it probably belongs at the repository level instead.

## Working Principles

- prefer one shared fix over many local copies
- keep generated outputs under `artifacts/`
- keep names durable and descriptive rather than tied to temporary planning
  language
- write docs for humans first, especially in the root handbook
- keep adapter and runtime contracts explicit instead of hiding them behind one
  convenience bucket

## Local Workflow

1. Read the root entry points in `README.md` and `docs/`.
2. Change the smallest surface that can honestly own the behavior.
3. Run the checks that match the area you touched.
4. Update contracts, docs, or metadata when the change affects shared behavior.
5. Commit in a small, coherent step with a durable message.

## Validation Expectations

Use the generated root help to discover the current command surface:

- `make help`
- `make list`
- `make list-all`

Typical validation flows:

- root docs and handbook changes: `make docs-check`
- root workflow or config changes: `make help`, `make list-all`, and the
  relevant root or package targets
- package-scoped changes: `make -C packages/<package> <target>`
- repository package dispatch changes: `make <target> PACKAGE=<package>`

If you touch checked-in API contracts, schemas, or publish flows, verify the
resulting artifacts rather than assuming the repository state stayed
consistent.

Do not keep empty root-owned surfaces alive with placeholder files. If this
repository does not yet own a checked-in API contract tree, leave `apis/`
absent instead of adding a root `README.md`.

## Commits

Use Conventional Commits with clear, durable intent.

Good examples:

- `build(repo): centralize pytest configuration`
- `docs(governance): clarify repository security reporting`
- `refactor(ci): replace duplicated package workflow wiring`

Avoid commit messages that depend on temporary project language or internal
iteration labels. Someone reading the history two years later should still know
what changed and why.

## Pull Request Quality Bar

- explain the user-facing or maintainer-facing effect
- call out compatibility impact when behavior or naming changed
- mention any follow-up that still remains open
- keep unrelated edits out of the same change

## When Shared Files Move

If you move or rename root-owned files, update the references that point to
them:

- package metadata and URLs
- workflow paths and automation
- docs links and handbook navigation
- contract or release references

## Getting Help

Open an issue for questions that are not security-sensitive. For security
issues, use the private reporting guidance in `SECURITY.md`.
