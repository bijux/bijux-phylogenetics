---
title: Test Strategy
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# Test Strategy

The repository does not rely on one test family.

It uses:

- runtime tests for behavior and packaging
- governance tests for repository contracts
- standards checks for shared sync integrity
- evidence-book checks for governed trust surfaces

That layered strategy matters because a green runtime test suite alone does not
prove the public repository is honest or publishable.
