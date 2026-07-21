---
title: Distribution And Package Map
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Distribution And Package Map

The repository publishes two runtime distributions and keeps its maintainer
automation in a third package. The split is about naming and responsibility,
not separate scientific editions.

```mermaid
flowchart TD
    source["Canonical runtime source<br/>bijux_phylogenetics"]
    canonical["bijux-phylogenetics<br/>canonical distribution"]
    alias["phylogenetic<br/>compatibility distribution"]
    dev["bijux-phylogenetics-dev<br/>repository tooling"]
    users["Python and CLI users"]
    maintainers["Docs, evidence, API,<br/>release maintainers"]

    source --> canonical --> users
    source --> alias --> users
    dev --> maintainers
    dev -. validates .-> canonical
    dev -. validates .-> alias
```

## Distribution Contracts

| Distribution | Import or command role | Intended consumer | Ownership boundary |
| --- | --- | --- | --- |
| `bijux-phylogenetics` | canonical `bijux_phylogenetics` runtime and CLI | analysts, library users, workflow authors | scientific runtime and public execution surfaces |
| `phylogenetic` | shorter compatibility installation and command story | users choosing the concise distribution name | same runtime behavior; no reduced feature tier |
| `bijux-phylogenetics-dev` | repository checks and generation utilities | maintainers and CI | documentation, evidence, schema, packaging, and release truth |

## Installed Identity Chain

```mermaid
flowchart LR
    request["Requested distribution"]
    metadata["Installed metadata<br/>and resolved version"]
    import["Import namespace"]
    command["Console command"]
    runtime["Canonical runtime<br/>contract"]

    request --> metadata --> import --> runtime
    metadata --> command --> runtime
```

| Requested distribution | Import | Command | Runtime owner |
| --- | --- | --- | --- |
| `bijux-phylogenetics` | `bijux_phylogenetics` | `bijux-phylogenetics` | `bijux_phylogenetics` modules in the canonical wheel |
| `phylogenetic` | `phylogenetic` facade, including `phylogenetic.api` | `phylogenetic` | the installed `bijux-phylogenetics` dependency reached through the facade |
| `bijux-phylogenetics-dev` | `bijux_phylogenetics_dev` | repository Make targets and maintainer modules | repository verification and release policy, never scientific execution |

For the alias, record both distribution versions. The requested
`phylogenetic` version identifies the facade artifact; the resolved
`bijux-phylogenetics` version identifies the computation it delegates to.

## Inspect Installed Identity

Canonical installation:

```bash
python -c "from importlib.metadata import version; print(version('bijux-phylogenetics'))"
python -c "import bijux_phylogenetics as p; print(p.__version__, p.__name__)"
bijux-phylogenetics --version
```

Compatibility installation:

```bash
python -c "from importlib.metadata import version; print(version('phylogenetic')); print(version('bijux-phylogenetics'))"
python -c "import phylogenetic as p; print(p.__version__, p.__name__)"
phylogenetic --version
```

The two metadata values in the alias environment answer different questions:
which facade was requested and which canonical runtime performed the work.
Retain both in a reproducibility record.

## Choice Rules

Install `bijux-phylogenetics` when you want the canonical name in an
environment specification, publication method, or support request. Install
`phylogenetic` when compatibility with the shorter distribution name is the
primary concern. Do not install both expecting independent runtimes.

The development package is not required for normal scientific use. Its
commands protect repository integrity; they are not an alternate analysis API.

## Source Ownership

The runtime implementation lives under
`packages/bijux-phylogenetics/src/bijux_phylogenetics/`. The alias package is a
distribution boundary around the canonical runtime contract. The development
package lives under `packages/bijux-phylogenetics-dev/` and may inspect runtime
metadata or artifacts without becoming their scientific owner.

```mermaid
flowchart LR
    alias_code["phylogenetic facade"] -->|declared dependency| runtime_dist["bijux-phylogenetics distribution"]
    runtime_dist --> runtime_code["bijux_phylogenetics runtime"]
    dev_code["bijux_phylogenetics_dev"] -. inspects and validates .-> runtime_dist
    dev_code -. does not execute analyses .-> runtime_code
```

The facade forwards curated top-level attributes and resolves compatible
submodules to the canonical runtime. Its local code owns alias import and
command behavior only; scientific modules remain owned by
`bijux_phylogenetics`.

## Compatibility Meaning

Compatibility requires more than both distributions installing successfully.
Their advertised imports, commands, version story, metadata, package data, and
documented behavior must remain aligned. Release checks should fail when the
alias silently exposes a smaller or stale contract.

Compatibility is evaluated at observable boundaries: import resolution,
public attribute discovery, parser construction, command dispatch, resource
access, structured serialization, and failure behavior. Matching a single
headline value is not an alias compatibility test.

## Where To Report A Problem

| Symptom | Likely owner |
| --- | --- |
| numerical result or runtime behavior | canonical runtime package |
| alias install/import/command mismatch | compatibility package boundary |
| generated docs, evidence index, or schema drift | maintainer package |
| external executable missing or incompatible | adapter plus the named engine environment |

Always include the installed distribution name and version. For adapter
failures, include the external executable version and captured diagnostic
without removing the provenance needed to reproduce the failure.

## Package Boundary Invariants

- scientific runtime modules do not depend on the maintainer package;
- the alias depends on and delegates to the canonical runtime;
- the maintainer package may validate public runtime locators without
  exporting repository-only evidence helpers through the runtime;
- wheels and source distributions contain their declared package data and do
  not absorb repository docs, tests, evidence records, or local artifacts;
- package versions, command identity, and public links remain consistent at
  release time.
