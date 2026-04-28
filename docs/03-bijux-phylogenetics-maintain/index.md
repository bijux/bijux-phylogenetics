---
title: Maintainer Overview
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-04-28
---

# Maintainer Overview

This repository consumes the shared Bijux standards surface through
`.bijux/shared/`, `makes/`, `configs/`, and the repository-owned
`bijux-phylogenetics-dev` package.

## Maintainer Rules


- treat `makes/`, `configs/`, and `.bijux/shared/` as managed standard surfaces
- run `make sync-badges` after editing badge templates or managed README blocks
- run `make sync-license-assets` after adding a new package or changing root legal assets
