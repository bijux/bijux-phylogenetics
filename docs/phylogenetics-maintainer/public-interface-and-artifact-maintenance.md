---
title: Public Interface And Artifact Maintenance
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Public Interface And Artifact Maintenance

The repository now documents several distinct public interface layers:

- surface-selection guidance
- CLI and Python reference pages
- native inference reference pages
- artifact contract and artifact consumption guides
- package metadata links to those guides

Those layers should change together when one public workflow, native contract,
or output family changes materially.

## Review Surfaces Together

When one public interface changes, review:

- `docs/phylogenetics/interfaces/index.md`
- `docs/phylogenetics/interfaces/surface-selection.md`
- the affected CLI, API, Python, or native reference pages
- `docs/phylogenetics/interfaces/artifact-contracts.md`
- `docs/phylogenetics/interfaces/artifact-consumption-guide.md`
- root and package README links
- package `pyproject.toml` URLs

## Questions To Ask

- can one reader still choose the right public surface without reading the full
  handbook
- do the reference pages still point to the right practical guides
- do artifact docs still name the stable output families honestly
- do package metadata links still reinforce the same public interface story
- has one new output family become stable enough that it now belongs in the
  public artifact-contract surface

## Failure Pattern

The usual failure is partial interface drift:

- one new runtime workflow exists
- the workflow guide gets updated
- interface selection still points readers elsewhere
- artifact docs never mention the new stable outputs
- package metadata links still describe the older surface map

That leaves the repository looking larger in code than it does in its public
contract story.
