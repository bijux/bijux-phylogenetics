---
title: Product Guide
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-10
---

# Bijux Phylogenetics Product Guide

`bijux-phylogenetics` is the public product guide for the repository runtime.
It explains what the phylogenetics surface does today, where the product
boundary stops, how the runtime is organized, and which interfaces are meant to
be stable enough for serious reuse.

This guide is for readers who want the product answer first, not the source
tree answer first. It should let someone understand the repository without
digging through package names, maintainer scripts, or evidence internals.

<div class="bijux-quicklinks">
  <a class="md-button md-button--primary" href="../../index.md">Open the documentation home</a>
  <a class="md-button" href="foundation/">What this runtime is for</a>
  <a class="md-button" href="architecture/">How the runtime is organized</a>
  <a class="md-button" href="interfaces/">Commands and public contracts</a>
  <a class="md-button" href="operations/">Install and operate locally</a>
  <a class="md-button" href="quality/">Quality gates and limits</a>
</div>

## What This Guide Covers

- what the runtime publishes and what it refuses to claim
- how the public runtime boundary differs from the evidence-book boundary
- which packages, commands, and artifacts are public-facing
- how to install, inspect, and rerun common workflows responsibly
- which quality and publication gates protect the runtime surface

## Read By Question

- what is this repository actually trying to do:
  [foundation](foundation/index.md)
- how are runtime responsibilities split across packages and layers:
  [architecture](architecture/index.md)
- which command-line and Python surfaces should I treat as stable:
  [interfaces](interfaces/index.md)
- how do I get a working local setup and run common workflows:
  [operations](operations/index.md)
- how do I judge trust, limits, and release readiness:
  [quality](quality/index.md)

## Public Versus Internal

The public handbook explains the product and runtime surface. Internal
maintainer rules, managed standards behavior, release operations, and repo
health contracts are intentionally kept separate under the internal docs area.

The evidence-book is also separate. It is the scientific trust surface for
governed parity and coverage boundaries, not the general runtime handbook.

## Start Here

- [foundation](foundation/index.md) for scope, ownership, and non-goals
- [architecture](architecture/index.md) for runtime layout and boundary design
- [interfaces](interfaces/index.md) for CLI, Python, and artifact surfaces
- [operations](operations/index.md) for install, setup, and practical usage
- [quality](quality/index.md) for tests, release gates, and honesty rules
