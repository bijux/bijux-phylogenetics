---
title: Validation Lanes
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Validation Lanes

The repository does not use one generic test lane because different questions
need different proof surfaces.

## Start By Question

- does one owned runtime surface behave correctly in isolation
  Look at runtime tests.
- do repository governance and release contracts still hold
  Look at governance and readiness checks.
- do wrapper contracts still parse and execute as expected
  Look at external-engine contract and real-engine lanes.
- do stored scientific outputs still match the governed workflow contract
  Look at scientific-validation lanes.
- does one owned workload still fit the expected runtime and memory envelope
  Look at benchmark and stress lanes.

## Main Lane Families

- runtime behavior lanes
  Exercise owned Python and CLI behavior, packaging, and stable runtime
  contracts.
- governance lanes
  Exercise repository boundaries, managed standards, and release or publish
  rules.
- external-engine contract lanes
  Exercise fake-engine and parser expectations without requiring real local
  binaries.
- real-engine execution lanes
  Exercise compact real runs against actual local tools when those executables
  are present.
- scientific-validation lanes
  Re-run governed workflows and compare the durable reviewer-facing outputs.
- benchmark and stress lanes
  Measure runtime, memory, output size, and owned workload pressure.

## Reading Rule

Do not collapse those lanes into one generic confidence score.

- a runtime unit or integration test does not prove broader scientific closure
- a fake-engine contract lane does not prove the real executable still works
- a real-engine run does not prove every stored governed output still matches
- a scientific-validation rerun does not replace benchmark or release review
- a benchmark pass does not turn one approximation or wrapper lane into a
  native claim

## Practical Relation To The Rest Of The Handbook

- use [test strategy](test-strategy.md) for the overall repository testing
  approach
- use [release review workflow](release-review-workflow.md) when the question is
  whether one repository state looks publishable
- use [native benchmark evidence](native-benchmark-evidence.md) when the
  question is how to interpret native benchmark outcomes
- use the [evidence-book overview](../../phylogenetics-evidence-book/index.md)
  when the question becomes claim closure rather than lane structure
