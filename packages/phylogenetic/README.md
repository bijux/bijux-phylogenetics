# phylogenetic

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/phylogenetic/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)
[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)

[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)
[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)

[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/phylogenetics/)
[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

Compatibility alias package for `bijux-phylogenetics`.

This distribution installs the canonical phylogenetics runtime under the
shorter `phylogenetic` package and CLI names. The runtime behavior, report
contracts, and public API continue to live in `bijux-phylogenetics`; this
package exists so installations can use the shorter distribution name without
forking the implementation surface.

That includes the same documented native likelihood, native inference,
benchmark, workflow, and reporting surfaces that the canonical package
publishes.

Choose this package when you want the shorter install and command names while
keeping the same public runtime contract.

## Install

`phylogenetic` supports Python 3.11 and newer.

```bash
python3.11 -m pip install phylogenetic
phylogenetic --help
```

## Read this next

- [Runtime product guide](https://bijux.io/bijux-phylogenetics/phylogenetics/)
- [Runtime architecture](https://bijux.io/bijux-phylogenetics/phylogenetics/architecture/)
- [Surface selection guide](https://bijux.io/bijux-phylogenetics/phylogenetics/interfaces/surface-selection/)
- [Artifact consumption guide](https://bijux.io/bijux-phylogenetics/phylogenetics/interfaces/artifact-consumption-guide/)
- [Native inference and benchmarks](https://bijux.io/bijux-phylogenetics/phylogenetics/interfaces/native-inference-and-benchmarks/)
- [Native maximum-likelihood workflows](https://bijux.io/bijux-phylogenetics/phylogenetics/operations/native-maximum-likelihood-workflows/)
- [Native Bayesian workflows](https://bijux.io/bijux-phylogenetics/phylogenetics/operations/native-bayesian-workflows/)
- [Native benchmark review](https://bijux.io/bijux-phylogenetics/phylogenetics/operations/native-benchmark-review/)
- [Release review workflow](https://bijux.io/bijux-phylogenetics/phylogenetics/quality/release-review-workflow/)
- [Evidence-book guide](https://bijux.io/bijux-phylogenetics/phylogenetics-evidence-book/)
- [Canonical runtime package on GitHub](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
- [Alias package changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/phylogenetic/CHANGELOG.md)
- [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
