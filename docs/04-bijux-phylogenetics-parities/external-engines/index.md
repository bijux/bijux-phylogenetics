---
title: External Engine Correspondence
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# External Engine Correspondence

MAFFT, trimAl, IQ-TREE2, FastTree, MrBayes, and BEAST remain external. They are
the computational owners of their algorithms. Bijux provides the governed
boundary around them:
input checks, executable discovery, explicit commands, failure capture,
required-output validation, typed parsing, checksums, manifests, diagnostics,
and evidence packages.

That is correspondence, not reimplementation. A successful adapter run means
that a named external engine produced an accepted artifact under a recorded
contract. It does not mean that Bijux independently reproduced the algorithm
or that another engine would produce the same scientific result.

```mermaid
flowchart LR
    inputs["Governed inputs<br/>checksums and semantics"]
    preflight["Executable preflight<br/>path and version"]
    execution["External engine<br/>command and exit state"]
    native["Native files<br/>stdout and stderr"]
    validation["Structural validation<br/>and typed parsing"]
    record["Manifest · checksums<br/>review artifacts"]

    inputs --> preflight --> execution --> native --> validation --> record
```

## Choose The Surface By Scientific Operation

| Operation | Computational owner | Bijux correspondence | Detailed guide |
| --- | --- | --- | --- |
| multiple-sequence alignment | MAFFT | ordinary and codon-aware execution, alignment validation, run manifest | [alignment and trimming](alignment-and-trimming.md) |
| alignment trimming | trimAl | named trimming modes, retained-site diagnostics, validated handoff | [alignment and trimming](alignment-and-trimming.md) |
| model selection and maximum likelihood | IQ-TREE2 | model/report parsing, topology, likelihood, support, reproducibility | [tree inference](tree-inference.md) |
| approximate large-tree inference | FastTree | streamed preflight, topology, resource record, SH-like support parsing | [tree inference](tree-inference.md) |
| Bayesian tree inference | MrBayes | NEXUS preparation and execution, trace/tree parsing, posterior review | [Bayesian inference](bayesian-inference.md) |
| Bayesian time-tree inference | BEAST | XML preparation and execution, calibration, log/tree and posterior review | [Bayesian inference](bayesian-inference.md) |

## Three Kinds Of Record

The same parser can be useful without proving the same execution history.
Every result must identify which record was created:

| Record | What happened | What can be claimed |
| --- | --- | --- |
| prepared analysis | Bijux wrote or validated NEXUS/XML/configuration; no posterior engine ran | the configuration satisfies the supported structural contract |
| parsed artifact | Bijux read governed native output created earlier or elsewhere | the supplied artifact parses and yields the reported observations |
| live workflow run | Bijux resolved the executable, recorded its version and command, observed the process, and validated current outputs | the named engine executed under the retained run identity |

File presence alone cannot promote a prepared or parsed record into a live-run
claim. Conversely, a live exit code of zero does not replace output validation.

## Preserve The Artifact Lineage

A normalized result is a child of native engine output, not a replacement for
it. The accepted record must make every transformation traversable in both
directions: a reviewer can move from a reported value to the native file and
run that produced it, and from a declared run to every accepted projection.

```mermaid
flowchart LR
    source["Source input<br/>checksum"]
    prepared["Prepared engine input<br/>checksum"]
    run["Engine execution<br/>identity"]
    native["Native output<br/>inventory"]
    parsed["Parsed observations<br/>and diagnostics"]
    normalized["Normalized result<br/>and downstream handoff"]

    source --> prepared --> run --> native --> parsed --> normalized
```

Each edge needs an explicit parent identity and transformation description.
For alignment workflows that includes translation, back-translation, and site
selection. For inference it includes model/search selection and tree parsing.
For posterior workflows it includes chain or state selection, burn-in, and
subsampling. If a link is absent, only the independently identifiable segment
can be reviewed; the end-to-end run is incomplete.

### Native Files Remain Authoritative

The native file is the authority for what the external engine emitted. A
parsed observation is an addressable reading of that file, and a normalized
result is a deliberately narrower projection for comparison or downstream
use. Neither projection may silently replace engine-specific structure.

Every projection records the native parent checksum, parser and schema
identity, extracted fields, unit or scale conversions, and values that could
not be represented. When a parser changes, the native artifact is parsed into
a new child record; the earlier projection remains attributable to the parser
that created it. This makes parser corrections reviewable without rewriting
the historical engine run.

## The Accepted Run Identity

A reviewable external-engine record joins four identities:

```text
run identity = inputs + engine + configuration + accepted outputs
```

- **Inputs:** absolute source identity, content checksums, taxon set, sequence
  type, partitions, and any tree or calibration inputs.
- **Engine:** requested and resolved executable, native version text, command,
  start/end state, exit code, timeout, and captured streams.
- **Configuration:** algorithm mode, model, seed, threads, replicate count,
  burn-in, and other engine-specific choices.
- **Outputs:** required native files, structural validation, checksums, parsed
  observations, and warnings.

`--resume` is accepted only when these identities still agree with a completed
manifest. An active marker, timeout, missing output, changed checksum, or
malformed artifact remains incomplete. The caller chooses whether to reject or
clean that state; it is never silently reclassified as success.

## Correspondence Does Not Merge Estimands

Some outputs can be placed beside one another only after their meanings are
retained. IQ-TREE2 ultrafast-bootstrap percentages, IQ-TREE2 SH-aLRT support,
FastTree SH-like local support, and posterior clade probability are different
statistical quantities. A display may normalize them to a shared range, but
every row must keep its engine, method, scale, and denominator.

The same rule applies to models and posteriors. A common tree topology does not
make likelihood functions, priors, search procedures, or support estimands
identical.

## Know What Was Actually Run

Start with the adapter inventory and the command-specific help:

```bash
bijux-phylogenetics adapter --help
bijux-phylogenetics adapter align --help
bijux-phylogenetics adapter compare-engines --help
bijux-phylogenetics adapter beast-run --help
```

For a completed run, inspect the manifest rather than reconstructing identity
from filenames. The durable review packet includes the manifest, all required
native outputs, their checksums, parser summaries, diagnostics, and the exact
artifact passed downstream.

## Claim Ladder

External-engine evidence becomes stronger one rung at a time:

1. the input and configuration are structurally valid;
2. the executable identity and command are known;
3. the process completed and required outputs passed validation;
4. parsed scientific quantities are internally coherent;
5. repeated or cross-engine observations satisfy a declared comparison rule;
6. an Evidence Book claim governs the observations, denominator, tolerance,
   and verdict.

Do not jump from executable availability to scientific validation. The
[correspondence verification matrix](../correspondence/verification-matrix.md)
explains how adapter records relate to live parity cases and Evidence Book
verdicts.
