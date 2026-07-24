---
title: Study-Specific R References
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Study-Specific R References

Some R packages enter the repository through a named scientific study rather
than a general-purpose live parity registry. Their role is narrower: they may
own a reference calculation, establish a structural mapping, or identify an
open computational boundary. The Evidence Book preserves that role without
turning package presence in an R script into a blanket compatibility claim.

| Reference | Governed role | Current claim | Explicit limit |
| --- | --- | --- | --- |
| `nlme` | GLS optimization over covariance structures created with APE | checked PGLS coefficients, likelihoods, and model comparisons for named primate fixtures | no standalone live `nlme` registry |
| `treeio` | node lookup and tree–trait joining in the primate longevity study | verified node identifiers and joined tip count for the governed fixture | no general tree-object API equivalence claim |
| `tidytree` | tabular tree representation in the primate reference workflow | package version and reference provenance retained with the structural comparison | no broad verb-by-verb parity registry |
| `ggtree` | course-study visualization surface | source fragments and plot intent retained | figures are plot-only; no pixel or aesthetic parity verdict |
| `MCMCglmm` | posterior computation in the cichlid mixed-model studies | source-backed external computation boundary | no governed reproduction of the PCM4 or PCM5 posterior core |

## Record Origin Matrix

| Origin | Execution identity | Permitted statement |
| --- | --- | --- |
| live R registry | current R process, package version, registered case and fixture | the selected cases passed, failed, or skipped under their observation rules |
| checked reference fixture | frozen output identity, producing environment, input identity, and comparison code | the named fixture corresponds under its checked rules |
| structural study record | source fragment, dataset, tree, and structural observations | the named join, node, or table contract was observed |
| plot-only source | source code and intended figure semantics | the visualization provenance is known; analytical or pixel correspondence is not claimed |
| external computation boundary | source and model identity without comparable execution | the intended method is located, but its numerical claim is `not_comparable` |

These origins are not interchangeable. Re-running an adjacent native method
does not turn an external boundary into a live registry result.

## APE Plus `nlme`

The PGLS reference is a composition, not a claim about either package in
isolation. APE constructs Brownian or Pagel-lambda correlation structures and
`nlme::gls` fits the corresponding regression. Bijux compares the scientific
outputs that the checked fixture declares: coefficients, log-likelihood,
information criteria, estimated lambda where applicable, and selected
diagnostics.

That comparison is reproducible from repository-owned inputs and frozen
reference outputs. It is not labelled `nlme-live`, and a successful checked
fixture does not establish coverage for every `gls` correlation structure,
optimizer, missing-data policy, or variance model.

## `treeio`, `tidytree`, And `ggtree`

The primate longevity dossier separates structural data correspondence from
rendered presentation:

- `treeio` supplies reference node lookup and tree–trait joins.
- `tidytree` supplies the reference tabular tree representation and recorded
  package context.
- `ggtree` fragments describe the intended visual analyses, but remain
  plot-only because no governed image-comparison contract exists.

For the governed primate fixture, the structural comparison records a joined
tip count of 75 and exact reference node identifiers for `Pan_paniscus`,
`Hylobates_lar`, and `Node32`. Those observations support the named study
claim; they do not imply comprehensive compatibility with Bioconductor tree
objects.

## `MCMCglmm`

PCM4 and PCM5 use `MCMCglmm` for Gaussian and generalized phylogenetic mixed
models, including random-effect covariance, posterior summaries, and chain
diagnostics. Bijux has adjacent native Bayesian, comparative-regression, and
external MrBayes/BEAST capabilities, but none is presented as a substitute for
that model family.

The Evidence Book therefore retains the source programs, datasets, governed
claim records, and `not_comparable` verdicts while identifying `MCMCglmm` as
the unresolved computational owner. Closing that boundary would require an
independent governed runner, aligned priors and parameterization, posterior
summary rules, convergence criteria, and reproducible comparison artifacts.

## Promotion Into A Study Verdict

Mentioning or loading an R package is provenance, not verification. A
study-specific result becomes a parity claim only when the repository can name
the input, executed implementation, versions, compared observations,
acceptance rule, artifacts, and verdict. Until then, the package remains a
documented dependency or open boundary.

Reference-package evidence enters an Evidence Book verdict only when a bundle
binds it to a claim ID, exact inputs, source and runtime identities, checks,
eligible observations, limitations, freshness dependencies, and one governed
verdict. The bundle must account for failed, skipped, excluded, and absent
outputs rather than copying only successful rows.

A live harness case can remain passed while the study verdict is
`not_comparable` because the study uses a different model, dataset, or claim.
A checked fixture can support `matched_with_tolerance` without implying a live
R execution on the reviewer's machine. State both records independently.
