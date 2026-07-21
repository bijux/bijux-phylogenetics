---
title: MAFFT And trimAl Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# MAFFT And trimAl Correspondence

MAFFT owns multiple-sequence alignment and trimAl owns external alignment
trimming. Bijux owns the surrounding contract: input validation, executable
resolution, argument construction, captured output, stage fingerprints,
alignment validation, diagnostics, and workflow artifacts.

```mermaid
flowchart LR
    fasta["Unaligned FASTA"]
    preflight["Sequence and taxon<br/>preflight"]
    mafft["MAFFT"]
    aligned["Validated alignment"]
    trimal["trimAl"]
    trimmed["Validated trimmed<br/>alignment"]
    review["Occupancy · length · taxa<br/>manifest · logs"]

    fasta --> preflight --> mafft --> aligned --> trimal --> trimmed --> review
```

## What Is Verified

| Boundary | Verified behavior |
| --- | --- |
| executable | resolved path, native version text, invocation and exit state are retained |
| input | FASTA parsing, taxon identity, sequence type, non-empty content and coding checks where requested |
| MAFFT output | a non-empty parseable alignment is produced with the expected taxon set |
| trimAl output | a non-empty parseable alignment is produced and before/after dimensions remain reviewable |
| workflow state | stage inputs, configuration, output paths, fingerprints and failures are recorded |
| downstream handoff | the exact validated alignment used by inference is identified in the manifest |

Success does not mean the alignment is biologically optimal. Strategy choice,
gap treatment, trimming threshold, coding frame, and downstream sensitivity
remain study decisions.

## Use The Adapters

```bash
bijux-phylogenetics adapter align \
  dataset/sequences.fasta \
  --out artifacts/engines/aligned.fasta

bijux-phylogenetics adapter trim \
  artifacts/engines/aligned.fasta \
  --out artifacts/engines/trimmed.fasta
```

Use the command help for version-specific options and executable overrides.
`adapter fasta-to-tree` composes alignment, optional trimming, and inference
when one governed workflow bundle is preferable to separate stages.

## Correspondence, Not Reimplementation

Bijux does not reproduce MAFFT's alignment algorithm or trimAl's trimming
algorithm and does not relabel their results as native. The correspondence
claim is that a recorded engine run produced a validated artifact through a
known adapter contract. A numerical parity claim would require an independent
reference output and a separate acceptance rule.
