---
title: Release Review Workflow
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Release Review Workflow

Use this workflow when the question is whether one repository state looks
publishable, not merely whether one command or one test lane passed.

## Start With Governed Repository Reports

The repository already exposes governed readiness surfaces through root make
targets:

```bash
make report-publish-readiness
make check-publish-readiness
make report-release-readiness
make check-release-readiness
```

Those commands write governed report artifacts under `artifacts/root/` and
separate scorecard production from hard-gate enforcement.

## What Each Surface Answers

- `report-publish-readiness`
  Builds the repository publish-readiness scorecard without acting as the final
  release gate.
- `check-publish-readiness`
  Fails when the repository still has unresolved publish-readiness issues.
- `report-release-readiness`
  Builds the broader release-readiness evidence and governance surface.
- `check-release-readiness`
  Enforces the full release-readiness gate.

Use the `report-*` surfaces when you are auditing or preparing. Use the
`check-*` surfaces when the repository should already be held to the gate.

## Add The Runtime Truth Review

When the release question includes full-suite test truth and real-engine test
truth together with owned workflow outputs, run the governed CLI summary:

```bash
bijux-phylogenetics report release-truth \
  --test-report artifacts/pytest/full-suite.xml \
  --real-engine-test-report artifacts/pytest/real-engine.xml \
  --out artifacts/release-truth-report.html \
  --json
```

That surface complements the repository readiness reports by reviewing the
runtime and validation story in one place.

## Read The Release Story In This Order

For a serious release review, read:

1. public package and runtime docs
2. package metadata and package roles
3. publish-readiness and release-readiness reports
4. release-truth workflow output when available
5. evidence-book freshness, integrity, and coverage reports

That order makes it easier to catch contradiction before treating a green gate
as enough.

## Release Questions To Ask Explicitly

- do the README, package READMEs, and public docs describe the same runtime
  truth
- do metadata, maturity labels, and public links describe the same product
  surface
- do native workflow claims stay bounded by documented native contracts
- do wrapper-backed surfaces remain labeled as wrappers
- do evidence freshness and coverage reports support every scientific claim
  made for this release

## What This Workflow Does Not Replace

- it does not replace package build verification
- it does not replace benchmark or parity review
- it does not replace the evidence-book
- it does not turn one green check into full scientific closure
