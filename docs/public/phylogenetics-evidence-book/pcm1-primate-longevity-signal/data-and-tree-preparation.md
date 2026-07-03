---
title: PCM1 Data and Tree Preparation
audience: public
type: study-guide
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM1 Data And Tree Preparation

The first half of PCM1 is about getting the study state into a form that later
comparative work can honestly use.

That includes:

- loading the primate table into a governed shape
- making type repair and missing-data decisions explicit
- importing and checking the trimmed tree
- proving that tree tips and trait rows correspond before model fitting

This matters because later signal results are only meaningful if those earlier
joins and exclusions are already visible.

## What This Part Answers

- are the primate traits in one governed shape before model work starts
- does the trimmed tree still match the trait table
- are missing-value and exclusion decisions visible before any signal claim is
  made

## Why This Is A Trust Surface

This part of PCM1 is not only data cleaning. It is the proof that later
comparative outputs are anchored to one reviewable prepared study state rather
than one hidden preprocessing step.

## Read Next

- [signal and ancestral results](signal-and-ancestral-results.md)
