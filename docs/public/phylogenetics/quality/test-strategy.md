---
title: Test Strategy
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-13
---

# Test Strategy

The repository does not rely on one test family.

It uses:

- runtime tests for behavior and packaging
- governance tests for repository contracts
- standards checks for shared sync integrity
- evidence-book checks for governed trust surfaces
- external-engine contract tests for adapter and parser expectations
- real local engine execution lanes for MAFFT, trimAl, IQ-TREE2, FastTree,
  MrBayes, and BEAST when those executables are present
- governed scientific-validation lanes that rerun checked real workflows
  against stored reviewer-facing outputs

That layered strategy matters because a green runtime test suite alone does not
prove the public repository is honest or publishable.

## External Engine Lanes

The repository keeps three distinct boundaries for external tools.

- `engine_contract` tests exercise fake-executable and parser contracts. They
  prove the owned Python surfaces call the engines and parse their artifacts the
  way the repository expects, without requiring local scientific binaries.
- `engine_real` tests under `packages/bijux-phylogenetics/tests/real_local/`
  execute the real binaries on compact inputs. They prove MAFFT alignment,
  trimAl trimming, IQ-TREE2 inference and bootstrap support, FastTree
  approximate inference, MrBayes acceptance, and BEAST parsing against the
  actual executables.
- `scientific_validation` tests rerun governed reference workflows and compare
  the stored reviewer-facing outputs. They exist to catch silent drift in the
  durable public contract, not just API regressions.

Those boundaries are intentional. A fake-engine unit test should not be treated
as proof that the scientific binary still works locally, and a real-engine
rerun should not be treated as proof that every stored reference artifact still
matches the governed contract.

## Local Commands

Use the repository make surfaces when you want the separated lanes explicitly.

- `make test-external-engines` runs the real executable integration lane and the
  non-slow scientific-validation lane.
- `make test-scientific-validation-slow` runs the slower governed FASTA-to-tree
  reference reruns separately, with a longer timeout budget.

That split keeps the routine local engine check reviewable while preserving a
durable place for slower golden-workflow validation.
