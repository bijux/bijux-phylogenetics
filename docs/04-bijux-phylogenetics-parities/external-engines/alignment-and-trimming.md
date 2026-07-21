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

## Stage Acceptance And Handoff

| Stage | Accept when | Preserve for the next stage |
| --- | --- | --- |
| sequence preflight | identifiers are unique, declared type is coherent, and exclusions are explicit | input hash, accepted/excluded taxa, alphabet and coding diagnostics |
| MAFFT execution | executable/version and command are known, exit state is accepted, and output parses | native streams, mode, output alignment, taxon and dimension checks |
| codon-aware alignment | translation policy, genetic code, exclusions, guide alignment, and back-translation are coherent | sequence type, codon decisions, protein guide identity, final nucleotide alignment |
| trimAl execution | input alignment identity matches the handoff and trimmed output parses | mode, threshold, before/after dimensions, removed-site summary |
| inference handoff | selected post-alignment artifact is unambiguous and complete | exact alignment path/hash, partition/coding state, manifest lineage |

The input to inference is the accepted handoff artifact, not whichever FASTA
file is newest in the output directory. A composed `fasta-to-tree` run must
retain the same stage identities and failures as the separate commands.

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

`--resume` is reuse under an identity check, not permission to trust arbitrary
existing output. Input, configuration, executable, manifest, and recorded
outputs must still match. Incomplete output is rejected or cleaned according
to the explicit incomplete-run policy; it is never silently promoted.

## Correspondence, Not Reimplementation

Bijux does not reproduce MAFFT's alignment algorithm or trimAl's trimming
algorithm and does not relabel their results as native. The correspondence
claim is that a recorded engine run produced a validated artifact through a
known adapter contract. A numerical parity claim would require an independent
reference output and a separate acceptance rule.
