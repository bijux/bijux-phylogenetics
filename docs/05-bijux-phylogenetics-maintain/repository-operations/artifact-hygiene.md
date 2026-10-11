---
title: Artifact Hygiene
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-21
---

# Artifact Hygiene

Repository execution writes transient environments, caches, reports, build
products, and run results under `artifacts/`. Checked-in evidence and governed
generated documentation remain in their declared repository-owned locations;
artifact cleanup must never treat them as disposable run output.

```mermaid
flowchart LR
    run["Repository run"] --> classify{"Classify output"}
    classify -->|transient| local["artifacts/root or<br/>artifacts/package"]
    classify -->|governed evidence| evidence["Checked-in bundle<br/>contract"]
    classify -->|governed generated| generated["Declared generated<br/>destination"]
    aliases["Disposable root aliases"] -.-> local
```

## Classify Before Writing

| Output class | Destination | Retention rule |
| --- | --- | --- |
| tool environment or cache | `artifacts/root/` | disposable; recreate from governed configuration |
| repository report | `artifacts/root/` | retain for the review interval; do not commit by default |
| package build, test, or benchmark output | `artifacts/<package>/` | disposable unless a named publication process promotes it |
| cross-package scratch output | `artifacts/tmp/` | disposable and safe to regenerate |
| checked-in Evidence Book artifact | governed bundle path | preserve provenance and validate with evidence checks |
| generated site/reference source | declared `docs/`, API, or schema destination | update only through its owning generator and integrity check |

Location follows ownership, not filename extension. A JSON report from a local
audit is transient; a JSON schema or evidence manifest can be governed source.

## Artifact Lifecycle Record

For any report retained beyond the command that produced it, record:

| Field | Why it is required |
| --- | --- |
| producer | identifies the owning command and package rather than inferring ownership from the path |
| revision and configuration | ties the report to the source and policy it evaluated |
| inputs | distinguishes a fresh run from reuse of an unrelated prior artifact |
| completion state | separates complete reports from diagnostic partial output |
| output inventory | makes expected, present, missing, failed, and skipped records reviewable |
| integrity | carries hashes or equivalent identity when another gate consumes the report |
| retention purpose | states which active review or release decision still depends on it |

A report without revision and producer identity is useful as local diagnostic
output but cannot serve as reconstructable release evidence.

## Promote, Retain, And Expire Deliberately

An artifact changes class only through an owning workflow. Copying a transient
report into a governed directory does not give it provenance, and leaving a
file under `artifacts/` does not make it eligible for indefinite reuse.

| Lifecycle decision | Required record | Filesystem action |
| --- | --- | --- |
| retain for an active review | producer, revision, inputs, completion, inventory, integrity, consumer, and review purpose | keep an immutable review packet under its run identity |
| promote into governed output | destination owner, schema, generator, provenance, validation, freshness, and replacement policy | write through the governing producer, then inspect and commit the declared destination |
| supersede a retained report | old and new identities, reason, consumers, and whether conclusions changed | preserve the old packet until dependent decisions move or expire |
| expire transient output | verified absence of active consumers and governed references | remove only the explicit artifact-owned path |
| reject contaminated output | conflicting parent identities, stale files, or incomplete inventory | preserve diagnostic evidence if needed, then recreate in a clean run root |

Retention is consumer-aware. A release-readiness report cannot be deleted
while it is the cited basis of an active decision, and it cannot be reused for
a later revision because the filename is stable. A governed evidence artifact
is changed through its dossier workflow, never through transient cleanup.

Before deletion, resolve symlinks and parent identities, enumerate active
consumers, and verify that no manifest or decision record points to the target.
After promotion or expiration, run the owner-specific integrity check and
confirm the worktree exposes every intended governed change.

## Canonical Paths

Prefer the canonical artifact-owned paths directly when a maintainer command needs a concrete executable or output root:

- root environment: `artifacts/root/check-venv/`
- tox work directory: `artifacts/root/tox/`
- hypothesis state: `artifacts/root/hypothesis/`
- benchmark scratch space: `artifacts/root/benchmarks/`
- documentation site: `artifacts/root/docs/site/`

Those are the durable contract paths. Root aliases are compatibility links;
their targets remain the canonical locations.

Commands may replace transient reports at canonical paths. A reviewer who
needs to retain one decision must copy or archive the complete report packet
through a governed process, including its revision and input identities; a
screenshot or isolated JSON fragment is not sufficient.

## Root Compatibility Aliases

The root `.venv`, `.tox`, `.hypothesis`, and `.benchmarks` paths should be
symlink aliases into `artifacts/root/`.

- `make setup` creates missing aliases and never deletes a real directory
- `make clean-root-artifacts` preserves aliases and their canonical targets
- maintainers should not depend on them as the only route to governed
  environments or caches

If one alias is a real directory, inspect its contents and size before using
the explicit artifact-alias migration command. Ordinary setup and cleanup must
not replace or delete it.

The governed audit also rejects root `.pytest_cache`, `.ruff_cache`, and
`.mypy_cache` directories. Run `make check-artifact-governance` to verify root
aliases, forbidden caches, tox work location, and Make target output paths.

## Safe Cleanup

Use `make clean-root-artifacts` for ungoverned root caches. Treat cache and
run-product removal as a separate, explicitly scoped operation.
Before broader cleanup:

1. resolve every target to an explicit path beneath this repository's
   `artifacts/` tree;
2. distinguish a symlink from a real directory and inspect its target;
3. preserve reports required for an active review or release decision;
4. exclude checked-in evidence, schemas, documentation sources, and generated
   destinations governed by another command;
5. run the narrow producing command again when reproducibility is being
   verified;
6. confirm `git status` and `make check-artifact-governance` after cleanup.

Never infer a destructive cleanup target from an unresolved environment
variable, broad glob, workspace root, or parent directory.

## Review Checklist

- Does a new command write into `artifacts/` rather than the root?
- Does a new guide point to the canonical artifact path instead of only to
  `./.venv`?
- Does cleanup remove disposable aliases and stray caches without touching
  checked-in evidence?
- Can every retained report be tied to the revision and command that produced
  it?
- Does a producer create its destination rather than depending on unrelated
  prior runs?
- Does the worktree remain free of untracked run products after verification?

## Failure Conditions

Treat a real directory at a governed alias path, a forbidden root cache,
output outside the owned artifact tree, cleanup that can reach checked-in
records, an untraceable release report, or a command that relies on stale
prior output as an operational defect. Do not hide the defect with ignore
rules; correct the producer, alias, or cleanup boundary.

Also treat cross-run contamination as a defect: a command must not accept a
prior output merely because it occupies the expected path. Resume or reuse
requires the producer's identity checks over inputs, configuration, executable
where relevant, completion state, and recorded output inventory.
