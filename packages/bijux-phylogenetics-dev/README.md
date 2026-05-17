# bijux-phylogenetics-dev

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://github.com/bijux/bijux-phylogenetics)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)
<!-- bijux-phylogenetics-badges:generated:end -->

Repository-owned developer tooling and automation helpers for
bijux-phylogenetics.

This package owns repository checks for docs synchronization, badge rendering,
license asset management, quality gates, and release support. It is intended
for maintainers working from the repository workspace, not as the runtime entry
point for end users.

Its packaging surface now also owns the installed-artifact smoke proof that
builds the runtime wheel and sdist, installs each into a clean virtual
environment, copies packaged example inputs through the installed runtime API,
and verifies packaged example and dataset resources through the public CLI.

## Read this next

- maintainer handbook: [Maintainer handbook](https://bijux.io/bijux-phylogenetics/internal/maintain/)
- source directory: [Developer tooling source directory](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics-dev)
- changelog: [Developer tooling changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/bijux-phylogenetics-dev/CHANGELOG.md)
- security policy: [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
