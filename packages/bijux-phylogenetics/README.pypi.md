# bijux-phylogenetics

Canonical runtime package for the `bijux-phylogenetics` repository.

`bijux-phylogenetics` is the full public Python and CLI surface for tree
validation, alignment handling, native likelihood and inference, comparative
analysis, ancestral reconstruction, explicit parsimony, and evidence-linked
report generation.

Choose this package when you want the full runtime contract rather than the
shorter compatibility alias.

The full handbook includes diagrams, architecture maps, and study navigation
on the documentation site. This PyPI page stays text-first so it renders
cleanly in the package index.

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

## What This Package Covers

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

## Why This Package Is More Than A Thin Wrapper

- it owns the `PhyloTree` runtime and the topology, support, and tree-distance
  semantics that sit underneath multiple workflows
- it exposes typed Python workflow results with stable artifact writers instead
  of only printing shell output
- it includes native maximum-likelihood tree inference results and native
  Bayesian public inference entry points as documented public surfaces
- it keeps wrapper-backed orchestration honest by separating native ownership
  from external-engine execution

## Trust And Evidence

The runtime surface is broader than the closed evidence surface, and the
evidence surface is narrower than the full runtime surface. That distinction
is deliberate:

- a documented runtime surface is usable
- a native surface is locally implemented rather than only wrapped
- an evidence-book study explains how much trust a scientific claim currently
  carries

Use the package for execution. Use the docs and evidence-book when you need
quality boundaries, parity review, or study-specific provenance.

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

## Choose This Package When

- you want the canonical package name for notebooks, papers, or downstream
  integration
- you want both CLI workflows and typed Python workflow results
- you want the full public runtime contract rather than the shorter alias
  surface

## Read Next

- [Runtime guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- [Surface selection guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/surface-selection/)
- [Artifact consumption guide](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/artifact-consumption-guide/)
- [Native Inference And Benchmark Surfaces](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/native-inference-and-benchmarks/)
- [Native maximum-likelihood workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-maximum-likelihood-workflows/)
- [Native Bayesian workflows](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-bayesian-workflows/)
- [Native benchmark review](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-benchmark-review/)
- [Evidence-book guide](https://bijux.io/bijux-phylogenetics/02-bijux-phylogenetics-evidence-book/)
