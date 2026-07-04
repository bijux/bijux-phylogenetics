---
title: Runtime Operations
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# Runtime Operations

This section covers installation, local execution, and practical workflow use.

Use it when you want to get work done with the runtime rather than study the
interface map or release posture.

The operations docs should make the package feel runnable at real scale. That
means they need to cover workflow-first usage, native-first usage, and
review-first usage instead of pretending every reader follows one path.

## Operational Reading Order

1. [Installation and setup](installation-and-setup.md)
2. [Common workflows](common-workflows.md)
3. [Native maximum-likelihood workflows](native-maximum-likelihood-workflows.md)
4. [Native Bayesian workflows](native-bayesian-workflows.md)
5. [Native benchmark review](native-benchmark-review.md)
6. [Operational boundaries](operational-boundaries.md)

This section explicitly includes native maximum-likelihood workflows, native
Bayesian workflows, and native benchmark review so the flagship runtime
surfaces are discoverable from the main operations landing page.

## What This Section Helps You Do

- install the canonical runtime surface
- run documented workflows without scattering local outputs
- choose between workflow-oriented and native-oriented operating modes
- understand when a workflow is operationally convenient but scientifically
  review-heavy

## Operational Modes

| Mode | Best for | Representative guides |
| --- | --- | --- |
| Workflow-first | end-to-end runs that should emit reviewable artifacts | installation, common workflows, artifact guides |
| Native-first | direct use of owned maximum-likelihood or Bayesian contracts | native maximum-likelihood workflows, native Bayesian workflows |
| Review-first | understanding performance posture and benchmark limits | native benchmark review, operational boundaries |

## What Readers Should Be Able To Do After This Section

- install the canonical runtime without relying on local folklore
- run an end-to-end artifact-producing workflow
- switch from workflow use into lower-level native use when that is the right
  contract
- understand when an operationally valid run still needs careful scientific
  reading

## Why This Section Has To Be Better Than A Command Dump

The runtime now publishes enough depth that operations cannot stop at
`pip install` plus one toy command. These pages should help a reader operate:

- packaged examples
- workflow APIs that write governed artifacts
- native maximum-likelihood surfaces
- supported native Bayesian surfaces
- benchmark review paths that explain what a performance result really means
