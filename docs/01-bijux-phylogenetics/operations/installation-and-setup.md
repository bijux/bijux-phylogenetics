---
title: Installation And Setup
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Installation And Setup

Install the canonical runtime unless compatibility with the shorter
distribution name is a requirement.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install bijux-phylogenetics
```

The supported Python floor is 3.11. Record both Python and distribution
versions in reproducible work.

## Verify The Environment

```bash
python -c "import importlib.metadata as m; print(m.version('bijux-phylogenetics'))"
python -c "import bijux_phylogenetics; print(bijux_phylogenetics.__name__)"
bijux-phylogenetics --help
```

Verification separates three failures that are otherwise easy to confuse:
distribution installation, Python import, and command discovery.

## Canonical And Alias Distributions

```bash
python -m pip install phylogenetic
```

`phylogenetic` provides the shorter compatibility distribution for the same
runtime story. It is not a reduced edition and should not be installed beside
the canonical distribution to create two independent runtimes. See the
[package map](../architecture/package-map.md).

## External Engines

Some workflows invoke MAFFT, trimAl, IQ-TREE2, FastTree, MrBayes, or BEAST.
Installing the Python package does not install those executables. For an
adapter workflow:

1. install the engine through a channel appropriate to your platform;
2. confirm its executable is discoverable;
3. record the engine version;
4. run the adapter's smallest validation or help path;
5. preserve captured diagnostics with the analysis artifacts.

The runtime should report missing executables as adapter-environment failures,
not as native numerical failures.

## Pinning For Reproducibility

Exploration can use a compatible release range. A governed analysis should
record an exact runtime version and external-engine versions, plus input and
configuration identities. A lock file or environment specification is useful
only when it is retained with enough platform context to recreate the
environment.

```text
Python:                  3.11.x
bijux-phylogenetics:     exact installed version
external engine:         name and exact version, when used
platform:                OS and architecture
inputs:                  paths plus durable identity/checksum
configuration:           persisted parameters and seed policy
```

## Source Checkout For Contributors

Repository development uses the tracked lock and repository automation:

```bash
make install
make check
make docs-serve
```

Transient outputs belong under `artifacts/`. The root environment and docs
site are managed beneath `artifacts/root/`; do not treat generated output as
source documentation.

## Common Failure Diagnosis

| Symptom | Check first |
| --- | --- |
| package installs but import fails | active interpreter and environment |
| import works but command is missing | environment `bin` directory and command name |
| native method imports but adapter fails | external executable discovery and version |
| run starts but rejects data | taxon, format, rootedness, partition, or model validation |
| result differs after replay | versions, input identity, nondeterminism, seeds, and defaults |

After installation, use [surface selection](../interfaces/surface-selection.md)
and then a method-specific operations guide.
