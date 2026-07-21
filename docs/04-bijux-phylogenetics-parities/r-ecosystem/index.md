---
title: R Ecosystem Correspondence
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# R Ecosystem Correspondence

Three native parity harnesses execute checked-in R runners against live APE,
GEIGER, and phytools installations. Additional study fixtures retain narrower
R-derived comparisons, including `ape` plus `nlme` PGLS outputs and course
workflows that use packages such as `MCMCglmm` without claiming a live
drop-in-equivalence surface.

| Reference | Strongest governed relationship | Main scientific territory |
| --- | --- | --- |
| [APE](ape.md) | live case registry plus checked fixtures | trees, topology, distances, DNA utilities, covariance, contrasts, ancestral states |
| [GEIGER](geiger.md) | live model registry plus generated review report | continuous and discrete evolutionary-model fitting |
| [phytools](phytools.md) | live case registry plus stochastic-summary policies | signal, Mk models, simulation, stochastic maps, ancestral estimates |
| [`nlme`](study-references.md#ape-plus-nlme) | checked reference outputs with APE covariance | GLS/PGLS coefficients and model comparison |
| [`treeio`, `tidytree`, `ggtree`](study-references.md#treeio-tidytree-and-ggtree) | governed structural evidence or plot-only provenance | study-specific tree data and visualization surfaces |
| [`MCMCglmm`](study-references.md#mcmcglmm) | external study boundary | Bayesian phylogenetic mixed models in PCM4 and PCM5 |

## Execution Boundary

Live parity is opt-in through `RUN_LIVE_PARITY=1` and requires `Rscript` plus
the named package. The command emits structured summary and observation rows;
failed and skipped cases can be preserved under a failure-artifact directory.
The default parity command uses checked fixtures and does not imply that R ran
on the current machine.

## Why Several Reference Packages Matter

No single R package owns every convention used in phylogenetics. APE is the
structural baseline, GEIGER is a model-fitting reference, and phytools provides
specialized comparative and stochastic-history behavior. Agreement with one
does not silently grant agreement with another.

[Study-specific R references](study-references.md) explains the narrower
relationships that enter through Evidence Book workflows rather than a live
package registry.
