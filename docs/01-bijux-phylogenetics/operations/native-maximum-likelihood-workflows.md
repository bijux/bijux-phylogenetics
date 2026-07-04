---
title: Native Maximum-Likelihood Workflows
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Native Maximum-Likelihood Workflows

Use the native maximum-likelihood surface when you want one owned Bijux result
object for nucleotide tree inference instead of one broader wrapper-oriented
workflow bundle.

This public surface is Python-first today. The stable entry point is
`infer_nucleotide_maximum_likelihood_result_from_alignment(...)`, which builds
one `NucleotideMaximumLikelihoodResult` with the selected model, fitted
parameters, final tree, accepted search trace, optional support reports, and
optional benchmark metadata.

## Start From One Alignment Path

```python
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    infer_nucleotide_maximum_likelihood_result_from_alignment,
    write_nucleotide_maximum_likelihood_result_json,
)

result = infer_nucleotide_maximum_likelihood_result_from_alignment(
    Path("dataset/alignment.fasta"),
    model_name="auto",
    model_selection_criterion="aic",
    search_method="nni",
    start_tree_count=4,
    start_tree_seed=17,
)

print(result.model_name)
print(result.final_tree_newick)
print(result.final_log_likelihood)

write_nucleotide_maximum_likelihood_result_json(
    Path("artifacts/native-maximum-likelihood/result.json"),
    result,
)
```

Use `model_name="auto"` when the goal is governed model selection across the
supported nucleotide family. Use one explicit model such as `jc69`, `k80`,
`f81`, `hky85`, or `gtr` when the review question already fixes that choice.

## Add Support Evidence When It Matters

Bootstrap and SH-like support remain optional because they are more expensive
than the core tree search.

```python
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    infer_nucleotide_maximum_likelihood_result_from_alignment,
)

result = infer_nucleotide_maximum_likelihood_result_from_alignment(
    Path("dataset/alignment.fasta"),
    model_name="hky85",
    search_method="nni",
    start_tree_count=3,
    start_tree_seed=17,
    bootstrap_replicate_count=100,
    bootstrap_seed=5,
    sh_like_resampling_replicate_count=100,
    sh_like_resampling_seed=11,
)

assert result.support_reports.bootstrap is not None
assert result.support_reports.sh_like is not None
```

Use those support lanes when the output will be reviewed as one inferred tree,
not just as one intermediate optimization result.

## Attach Wrapper-Comparison Metadata Deliberately

The native ML result can also carry wrapper-correspondence benchmark metadata.
That is useful when the output needs one attached trust summary, but it should
not be confused with full evidence-book closure.

```python
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    infer_nucleotide_maximum_likelihood_result_from_alignment,
)

result = infer_nucleotide_maximum_likelihood_result_from_alignment(
    Path("dataset/alignment.fasta"),
    model_name="jc69",
    search_method="nni",
    include_wrapper_correspondence_benchmark=True,
)

assert result.benchmark_metadata is not None
print(result.benchmark_metadata.wrapper_correspondence_case_count)
```

When the real question is the repository's current native benchmark posture,
use the public evidence guide in
[native benchmark evidence](../quality/native-benchmark-evidence.md) instead of
inferring that answer from one attached result payload.

## Read The Result As One Stable Contract

The main fields to treat as stable are:

- `model_name`
- `model_selection_strategy`
- `model_selection_criterion`
- `final_tree_newick`
- `final_topology_fingerprint`
- `final_log_likelihood`
- `parameter_values`
- `search_method`
- `run_summaries`
- `search_trace_rows`
- `support_reports`
- `warning_messages`
- `benchmark_metadata`

That result object also supports `to_dict()`, `from_dict(...)`,
`write_json(...)`, and `load_json(...)` so one notebook or pipeline does not
need to reverse-engineer ad hoc output.

## Current Boundaries

- this public contract is strongest for native nucleotide maximum-likelihood
  tree inference
- this page does not claim that every lower-level likelihood helper has the
  same stability or public review weight
- this page does not turn wrapper-backed workflow orchestration into a native
  ownership claim
- this page does not replace the benchmark family when the review question is
  speed, memory, recovery accuracy, or wrapper comparison
