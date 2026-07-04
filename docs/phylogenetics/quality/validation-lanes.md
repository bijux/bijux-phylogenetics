---
title: Validation Lanes
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# Validation Lanes

The repository does not use one generic test lane because different questions
need different proof surfaces.

## Start By Question

| Question | Lane family |
| --- | --- |
| Does one owned runtime surface behave correctly in isolation? | Runtime tests |
| Do repository governance and release contracts still hold? | Governance and readiness checks |
| Do wrapper contracts still parse and execute as expected? | External-engine contract and real-engine lanes |
| Do stored scientific outputs still match the governed workflow contract? | Scientific-validation lanes |
| Does an owned workload still fit the expected runtime and memory envelope? | Benchmark and stress lanes |

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
  Re-run governed workflows and compare durable reviewer-facing outputs.
- benchmark and stress lanes
  Measure runtime, memory, output size, and workload pressure.

## Reading Rule

Do not collapse those lanes into one generic confidence score. a runtime unit or integration test does not prove broader scientific closure.
- A fake-engine contract lane does not prove the real executable still works.
- A real-engine run does not prove every stored governed output still matches.
- A scientific-validation rerun does not replace benchmark or release review.
- A benchmark pass does not turn an approximation or wrapper lane into a
  native claim.

## Related Reading

- [Test strategy](test-strategy.md)
- [Native benchmark evidence](native-benchmark-evidence.md)
- [Release review workflow](release-review-workflow.md)
- [Evidence Book](../../phylogenetics-evidence-book/index.md)
