# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

Canonical runtime package for the `bijux-phylogenetics` repository.

This package is the public Python API and CLI for tree validation, alignment
handling, comparative trait analysis, ancestral-state reconstruction,
discrete-state evolution analysis, external-engine orchestration, native
likelihood and inference, Bayesian posterior summarization, diversification
review, evidence bundles, explicit parsimony scoring, and HTML report
generation.

Choose this package when you want the full runtime contract rather than the
shorter compatibility alias.

## Package Structure At A Glance

- runtime package
- CLI workflows
- workflow Python API
- native runtime contracts
- likelihood foundations
- maximum-likelihood inference
- Bayesian inference
- benchmark contracts

## Install

`bijux-phylogenetics` supports Python 3.11 and newer.

```bash
python3.11 -m pip install bijux-phylogenetics
bijux-phylogenetics --help
```

The installed runtime also ships packaged example inputs:

```python
from pathlib import Path

from bijux_phylogenetics.core import copy_example_inputs

copy_example_inputs(Path("artifacts/example-inputs"))
```

## What This Package Already Covers

| Capability family | Examples |
| --- | --- |
| Trees | validation, rootedness review, support normalization, clade extraction, MRCA lookup, comparison, rendering, tree-set inspection |
| Alignments | FASTA validation, trimming, coding checks, translation, partition-aware input handling |
| Native likelihood | nucleotide, protein, codon, and discrete Mk finite-state likelihood foundations |
| Native inference | native maximum-likelihood tree inference results and native Bayesian public inference entry points |
| Comparative analysis | PGLS, signal, Brownian and OU modeling, discrete-state evolution, comparative diagnostics |
| Ancestral analysis | continuous and discrete reconstruction, uncertainty ledgers, transition review, report artifacts |
| Parsimony | Fitch, Wagner, Sankoff, Dollo, Camin-Sokal, ACCTRAN, DELTRAN, bootstrap, jackknife, NNI, SPR, ratchet |
| Reports and artifacts | reviewer-facing TSV, JSON, HTML, manifest, benchmark, and evidence-linked bundle surfaces |

## Why The Package Is Deeper Than A Typical Wrapper Toolkit

The package surface is broad because it joins several families that often live
in separate projects:

- owned tree and topology semantics
- workflow APIs that emit typed, reviewable results
- native maximum-likelihood and supported native Bayesian public surfaces
- comparative, ancestral, and parsimony families
- packaged datasets, example inputs, reports, and benchmark contracts

## Why This Package Is Not Just A Wrapper

- it owns the `PhyloTree` runtime and the topology, support, and tree-distance
  semantics that sit underneath multiple workflows
- it exposes typed Python workflow results with stable artifact writers instead
  of only printing shell output
- it includes native maximum-likelihood tree inference results and native
  Bayesian public inference entry points as documented public surfaces
- it keeps wrapper-backed orchestration honest by separating native ownership
  from external-engine execution

## Flagship User Paths

- use the CLI when you want one governed workflow with durable files
- use `bijux_phylogenetics.api` when you want typed workflow results in Python
- use native likelihood or Bayesian modules when you want lower-level owned
  runtime contracts
- use the benchmark and evidence-book guides when you need trust-oriented
  review rather than only execution

## Public Runtime Families

The public runtime families are:

- CLI workflows
- workflow Python APIs
- native runtime modules
- comparative analysis
- reports and artifacts

## Python Workflow Surface

The stable notebook-and-pipeline surface lives under `bijux_phylogenetics.api`.

It exposes typed workflow results for:

- FASTA validation
- multiple-sequence alignment
- alignment trimming
- full FASTA-to-tree execution
- maximum-likelihood tree inference
- branch-support estimation
- topology comparison
- PGLS comparative modeling
- discrete ancestral reconstruction
- reviewer-facing report generation
- config-driven workflow execution

```python
from pathlib import Path

from bijux_phylogenetics.api import (
    render_report_workflow,
    run_comparative_model_workflow,
    run_sequence_to_tree_workflow,
)

workflow = run_sequence_to_tree_workflow(
    Path("dataset/sequences.fasta"),
    out_dir=Path("artifacts/sequence-to-tree"),
    sequence_type="dna",
)

comparative = run_comparative_model_workflow(
    Path("dataset/tree.nwk"),
    Path("dataset/traits.tsv"),
    response="response",
    predictors=["predictor_one"],
    lambda_value=1.0,
)

report = render_report_workflow(
    tree_path=workflow.output_paths["tree"],
    alignment_path=workflow.output_paths["trimmed_alignment"],
    traits_path=Path("dataset/traits.tsv"),
    metadata_path=Path("dataset/metadata.tsv"),
    out_path=Path("artifacts/sequence-to-tree/report.html"),
)

workflow.write_json(Path("artifacts/sequence-to-tree/workflow.json"))
workflow.write_tsv(Path("artifacts/sequence-to-tree/workflow.tsv"))
comparative.write_tsv(Path("artifacts/comparative-model.tsv"))
```

## Native Inference And Benchmark Surfaces

Native Inference And Benchmark Surfaces are part of the public package
contract.

Public native surfaces include:

- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result(...)`
- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result_from_alignment(...)`
- `bijux_phylogenetics.phylo.likelihood.NucleotideMaximumLikelihoodResult`
- `bijux_phylogenetics.bayesian.run_bayesian_inference(...)`
- `bijux_phylogenetics.benchmark.benchmark_native_maximum_likelihood_suite(...)`

These surfaces matter because the package now owns native maximum-likelihood
tree inference results, native Bayesian public inference entry points, and
benchmark review contracts in addition to the higher-level workflow wrappers.

Representative public runtime surfaces include:

- `bijux_phylogenetics.api.run_sequence_to_tree_workflow(...)`
- `bijux_phylogenetics.api.run_comparative_model_workflow(...)`
- `bijux_phylogenetics.compare.topology`
- `bijux_phylogenetics.phylo.likelihood`
- `bijux_phylogenetics.bayesian`
- `bijux_phylogenetics.benchmark`

## Public Reading Rule

Runtime breadth and evidence closure are related, but they are not the same
claim.

- A documented runtime surface is usable.
- A native surface is locally implemented rather than only wrapped.
- A study or benchmark surface explains how much trust that claim currently
  carries.

That separation is deliberate. The package can expose substantial runtime depth
without pretending that every surrounding scientific claim is already closed by
the evidence-book.

## Read This Next

- [Surface selection guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/surface-selection/)
- [Artifact consumption guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/artifact-consumption-guide/)
- [Native Inference And Benchmark Surfaces](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/native-inference-and-benchmarks/)
- [Native maximum-likelihood workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-maximum-likelihood-workflows/)
- [Native Bayesian workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-bayesian-workflows/)
- [Native benchmark review](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-benchmark-review/)
- [Runtime architecture](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/architecture/)
- [Release review workflow](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/quality/release-review-workflow/)
- [Evidence-book guide](https://bijux.io/bijux-phylogenetics/02-bijux-phylogenetics-evidence-book/)
