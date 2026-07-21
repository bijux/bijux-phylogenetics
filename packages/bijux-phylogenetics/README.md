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

The canonical Python distribution for the Bijux phylogenetics runtime. It
provides owned tree, likelihood, inference, comparative, ancestral, parsimony,
simulation, dataset, and reporting contracts together with explicit adapters
for selected external engines.

## Installation

```bash
python -m pip install bijux-phylogenetics
bijux-phylogenetics --help
```

Python 3.11 or newer is required. External-engine workflows also require their
named executable—such as MAFFT, trimAl, IQ-TREE2, FastTree, MrBayes, or
BEAST—to be installed separately.

## Choose The First Scientific Route

| Starting material | Public route | First durable record |
| --- | --- | --- |
| unaligned FASTA | sequence-to-tree workflow | admitted-sequence report and stage-linked workflow result |
| tree plus continuous predictors | comparative model workflow | reconciled population, model specification and coefficient result |
| tree plus discrete trait | ancestral reconstruction workflow | state model, node identities and uncertainty-aware estimates |
| two candidate trees | topology comparison workflow | shared-taxon population, comparison mode and split observations |
| existing tree, alignment, traits or metadata | report workflow | machine manifest linked to the source scientific records |

The route is selected by the scientific object and intended result, not by the
shortest command. Native methods are computed by this distribution. Adapter
workflows retain Bijux orchestration while naming the external executable that
owns alignment, trimming, or inference.

| Decision boundary | Record that must survive |
| --- | --- |
| scientific object | admitted input identity and exclusions |
| method owner | native implementation or named external executable |
| computation | typed result, status, warnings, and diagnostics |
| persistence | structured outputs and complete manifest inventory |
| consumption | explicit acceptance or refusal for one bounded claim |

Do not skip admission because an example runs successfully. Packaged examples
prove installation and resource wiring; they do not approve a project dataset
or analysis design.

## Confirm Package Identity

Check the installed distribution and import before diagnosing scientific
behavior:

```python
from importlib.metadata import version

import bijux_phylogenetics

print(version("bijux-phylogenetics"))
print(bijux_phylogenetics.__version__)
print(bijux_phylogenetics.__all__)
```

The distribution version identifies the installed artifact. The import version
identifies the runtime loaded by Python. Investigate a disagreement before
trusting results; editable installs, environment shadowing, and mixed
environments can otherwise make a correct command execute the wrong code.

## Verify The Installed Artifact

Run the packaged demonstration before connecting project data:

```bash
bijux-phylogenetics demo run \
  --out artifacts/capability-demo \
  --manifest artifacts/capability-demo/manifest.json \
  --json
```

This command reads examples from the installed distribution, copies them into
a writable run directory, and produces tree, dataset, phylogenetic-input, and
topology-comparison reports. Inspect the manifest and the four reports under
`artifacts/capability-demo/reports/`. A successful run verifies package data,
runtime imports, parsing, comparison, reporting, and artifact linkage in the
current environment; it is not an inference-accuracy result.

## Runtime Map

| Layer | Responsibility |
| --- | --- |
| trees, alignments, traits, and models | validated scientific inputs |
| `bijux_phylogenetics.api` | curated workflow composition |
| native methods | Bijux-owned numerical computation |
| external-engine adapters | invocation, capture, parsing, and normalized results while preserving engine ownership |
| typed results | scientific fields, status, diagnostics, warnings, and refusal state |
| JSON, TSV, manifests, figures, and reports | persistent machine and review projections |

## Trace The Installed Result Lineage

A runtime result should remain traceable after the Python process exits:

| Lineage component | Required identity |
| --- | --- |
| runtime | installed distribution version, import version, and implementation owner |
| inputs | paths, content hashes, schemas, taxon/partition identity, and exclusions |
| configuration | model, method, transforms, seed policy, and execution settings |
| result | type/schema, status, diagnostics, warnings, and scientific fields |
| inventory | manifest identity and every accepted output checksum |
| projections | source result and transform behind each JSON, TSV, tree, figure, or report |

The manifest identifies the accepted projections; the typed result owns status,
warnings, and scientific fields. A projection that cannot be linked back to
the installed distribution, input identity, configuration, and typed result is
not a reconstructable runtime result.

For an external-engine route, insert executable path and version, invocation,
captured streams, native files, and parser state between configuration and the
normalized result. Normalization preserves ownership; it does not make the
external computation native.

## Capability Families

| Family | Examples |
| --- | --- |
| tree runtime | rootedness, support, clades, MRCA, splits, topology distance, tree sets, SVG rendering |
| alignment runtime | FASTA validation, trimming, translation, coding diagnostics, partitions |
| native likelihood | nucleotide, protein, codon, and discrete Mk finite-state models |
| native inference | maximum-likelihood tree results and supported Bayesian DNA inference |
| comparative methods | PGLS, Brownian and OU fits, signal, regression diagnostics, discrete evolution |
| ancestral methods | continuous/discrete reconstruction, uncertainty, stochastic mapping, tree-set review |
| parsimony | Fitch, Wagner, Sankoff, Dollo, Camin–Sokal, optimization, resampling, topology search |
| context and review | datasets, simulation, benchmarks, biogeography, diversification, reports, evidence bundles |

## What Installation Does And Does Not Supply

| Surface | Included in the Python distribution | Separate requirement |
| --- | --- | --- |
| native tree, likelihood, comparative, ancestral, parsimony, simulation, and supported Bayesian methods | yes | scientifically appropriate inputs and analysis design |
| packaged examples and demonstration inputs | yes | none for the installed-wheel proof |
| external-engine adapters and parsers | yes | the named executable and its runtime dependencies |
| R parity harness integration | Python-side contract | R, required packages, and explicit live-execution opt-in |
| governed study evidence | checked in with the repository and published through the Evidence Book | do not infer that every installed method has study evidence |
| HTML documentation and reports | generated presentation surfaces | structured artifacts for machine consumption |

Installing an adapter does not install or reimplement its external engine.
Installing the runtime does not convert Evidence Book boundaries into favorable
verdicts.

## Select A Public Surface

| Consumer | Use | Preserve |
| --- | --- | --- |
| Python application or notebook | curated `bijux_phylogenetics.api` workflows and documented domain imports | typed result, configuration, warnings, serialized artifacts |
| shell pipeline or scheduled execution | `bijux-phylogenetics` commands | command, exit status, JSON/TSV output, manifest, stderr |
| reviewer or downstream tool | versioned JSON, TSV, manifests, and linked reports | schema/version identity and complete output inventory |

Human-readable HTML is an inspection surface, not a stable data-exchange API.
Use structured output when another program needs to consume the result.

## Accept Or Refuse A Result

Do not equate file creation or a zero process exit with a usable scientific
result. Acceptance depends on the contract that produced it:

| Producer | Accept only when | Preserve on refusal |
| --- | --- | --- |
| Python workflow | result status, required output inventory, and method diagnostics satisfy the caller's policy | typed failure or partial state, warnings, configuration, and any valid diagnostics |
| CLI workflow | exit status agrees with structured status and every required artifact is present | stderr, JSON error record, manifest, and produced subset |
| native inference | model/search or chain state is complete enough for the intended interpretation | start/chain diagnostics, boundary estimates, rejected moves, and seed |
| external adapter | executable identity, invocation, native output, parser state, and normalized result are mutually consistent | captured streams, exit state, native partial output, and parse failure |

The caller owns the final suitability decision. The runtime exposes the state
needed to make it without treating warnings, skipped components, or partial
output as success.

## Publish A Reconstructable Result

Before a result leaves the analysis environment, retain the input identities,
runtime version, complete configuration, random seeds, typed or structured
result, warnings and diagnostics, required output inventory, and the acceptance
rule used for the scientific claim. For external execution, add executable
version, invocation, captured streams, native files, and parser state.

A Newick tree, coefficient table, or figure is a projection. Publish it beside
the structured record that explains how it was produced and where its
interpretation stops.

## Separate Reuse From Reinterpretation

Reusing a valid result does not permit changing its scientific meaning. A
consumer may create a new table, figure, or report when the projection remains
linked to the original status, model, units, ordering, exclusions, and
diagnostics. A changed transform, denominator, node mapping, factor reference,
or support interpretation creates a new analytical projection and must record
that rule.

When the original typed state is unavailable, consume the manifest and
versioned structured artifacts. Do not reconstruct a result from HTML, console
text, filenames, or row position. If the retained contract cannot express the
needed interpretation, refuse the reuse or rerun from the identified inputs.

## Workflow API

Use the curated `bijux_phylogenetics.api` namespace for application and
notebook workflows:

```python
from pathlib import Path

from bijux_phylogenetics.api import (
    render_report_workflow,
    run_comparative_model_workflow,
    run_sequence_to_tree_workflow,
)

sequence_run = run_sequence_to_tree_workflow(
    Path("dataset/sequences.fasta"),
    out_dir=Path("artifacts/sequence-to-tree"),
    sequence_type="dna",
)

comparative_run = run_comparative_model_workflow(
    Path("dataset/tree.nwk"),
    Path("dataset/traits.tsv"),
    response="response",
    predictors=["predictor"],
    lambda_value=1.0,
)

render_report_workflow(
    tree_path=sequence_run.output_paths["tree"],
    alignment_path=sequence_run.output_paths["trimmed_alignment"],
    traits_path=Path("dataset/traits.tsv"),
    metadata_path=Path("dataset/metadata.tsv"),
    out_path=Path("artifacts/sequence-to-tree/report.html"),
)

sequence_run.write_json(Path("artifacts/sequence-to-tree/workflow.json"))
comparative_run.write_tsv(Path("artifacts/comparative-model.tsv"))
```

Structured results retain information that console output cannot: diagnostics,
warnings, output identities, and serialization methods. Persist those records
with model configuration and input identity when a run must be reviewed.

External-engine results must also retain the executable name and version,
captured invocation, native output, and parser status. Normalization into a
Bijux result never transfers ownership of the external algorithm.

## Native Inference Entry Points

- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result(...)`
- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result_from_alignment(...)`
- `bijux_phylogenetics.phylo.likelihood.NucleotideMaximumLikelihoodResult`
- `bijux_phylogenetics.bayesian.run_bayesian_inference(...)`
- `bijux_phylogenetics.benchmark.benchmark_native_maximum_likelihood_suite(...)`

Native means that this package owns the computation. Tests, benchmarks, and
Evidence Book bundles separately define which numerical or study claims are
supported.

## Packaged Example Inputs

```python
from pathlib import Path

from bijux_phylogenetics.core import copy_example_inputs

copy_example_inputs(Path("artifacts/example-inputs"))
```

Packaged-resource smoke tests exercise this API from built wheels and source
distributions, not only from the repository checkout.

## Documentation

- [Product handbook](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- [Scientific domains](https://bijux.io/bijux-phylogenetics/02-bijux-phylogenetics-domains/)
- [Python surface](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/python-surface/)
- [Sequence-to-tree study](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/sequence-to-tree-study/)
- [Comparative and ancestral study](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/comparative-and-ancestral-study/)
- [Result and refusal semantics](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/result-and-refusal-semantics/)
- [Reproducibility records](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/reproducibility-records/)
- [CLI reference](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/cli-surface/)
- [Native maximum-likelihood workflow](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-maximum-likelihood-workflows/)
- [Native Bayesian workflow](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/native-bayesian-workflows/)
- [Artifact contracts](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/artifact-contracts/)
- [Evidence Book](https://bijux.io/bijux-phylogenetics/03-bijux-phylogenetics-evidence-book/)

Apache License 2.0.
