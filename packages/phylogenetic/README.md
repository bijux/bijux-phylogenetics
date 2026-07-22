# phylogenetic

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/phylogenetic/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)
[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)

[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)
[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)

[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`phylogenetic` is the short-name compatibility distribution for the complete
Bijux phylogenetics runtime. It depends on `bijux-phylogenetics` and exposes
the same workflow contract through `phylogenetic` import and command names.

It is not a fork, legacy edition, or reduced feature set.

```text
install/import phylogenetic
            |
            v
compatibility distribution and short-name facade
            |
            v
bijux-phylogenetics canonical runtime
            |
            v
canonical result, diagnostics, and artifacts

The result also retains the requested compatibility name.
```

The short name selects an entry route, not a second scientific computation.
Artifacts should preserve both the requested alias identity and the resolved
canonical runtime identity so another person can reconstruct what was
installed and what actually performed the work.

## Preserve Identity For Reproduction And Citation

The compatibility name and scientific owner serve different purposes. Keep
both where the route matters; use the canonical name where the implementation
must be unambiguous.

| Record | Identity to retain | Reason |
| --- | --- | --- |
| environment lock | `phylogenetic` request plus resolved `bijux-phylogenetics` version | reproduces dependency resolution |
| command transcript | invoked alias and canonical runtime reported by the result | reproduces the exercised facade and computation owner |
| manifest or support bundle | both distribution identities, loaded module path, and schema | distinguishes alias failure from runtime failure |
| methods or software citation | canonical `bijux-phylogenetics` identity, with alias invocation when operationally relevant | credits the scientific implementation without hiding the entry route |
| compatibility test | alias surface paired with the exact canonical surface | proves delegation rather than independent numerical agreement |

Do not cite `phylogenetic` as an independent implementation or count an alias
run as replication. Conversely, do not erase the alias from a defect report:
the invoked facade is the surface whose compatibility is being tested.

## Installation

```bash
python -m pip install phylogenetic
phylogenetic --help
```

Python 3.11 or newer is required.

## Choose The Compatibility Name Intentionally

| Situation | Recommended name | Reason |
| --- | --- | --- |
| new environment, publication or long-lived manifest | `bijux-phylogenetics` | names the scientific owner directly |
| existing code already imports `phylogenetic` | `phylogenetic` | preserves a supported consumer-facing route |
| concise interactive use | either | both resolve the same canonical runtime contract |
| independent implementation or evidence confirmation | neither interpretation | the alias is delegation, not independent corroboration |

Do not install both names to obtain two votes on a result. The compatibility
distribution depends on the canonical distribution and delegates scientific
work to it.

## Resolve Both Distribution Identities

The compatibility name and canonical runtime should be inspected together:

```python
from importlib.metadata import version

import bijux_phylogenetics
import phylogenetic

print("requested distribution:", version("phylogenetic"))
print("canonical distribution:", version("bijux-phylogenetics"))
print("alias import:", phylogenetic.__version__)
print("canonical import:", bijux_phylogenetics.__version__)
print("shared API facade:", phylogenetic.api is bijux_phylogenetics.api)
```

Record all four versions when reporting an installation problem. Local dirty
or editable builds can carry different local-version markers even when the
alias resolves to the intended canonical API; compare the installed artifacts
and runtime path before deciding whether the difference is a defect.

## Verify The Alias End To End

The compatibility command can run the same packaged demonstration:

```bash
phylogenetic demo run \
  --out artifacts/phylogenetic-capability-demo \
  --manifest artifacts/phylogenetic-capability-demo/manifest.json \
  --json
```

The command should report `status: ok` and write the same classes of input and
report artifacts as the canonical command. Record the resolved canonical
runtime version when reporting a defect; alias behavior is meaningful only in
relation to that runtime.

## Prove Delegation, Not Similarity

Alias correspondence is an identity chain:

| Link | Required correspondence |
| --- | --- |
| requested distribution | `phylogenetic` resolves the declared `bijux-phylogenetics` dependency |
| import and command facade | short names delegate to canonical public surfaces |
| canonical runtime | `bijux_phylogenetics` remains the scientific implementation owner |
| result | schemas, diagnostics, artifacts, resources, and refusal semantics remain canonical |

Verification therefore checks more than similar output values. It establishes
that the compatibility distribution resolves the intended canonical artifact,
the facade delegates to canonical public surfaces, packaged resources come
from the canonical runtime, commands preserve exit and stream behavior, and
serialized results retain the same schema and ownership.

If two commands produce similar trees while loading different canonical
artifacts, alias identity has failed. If they load the same runtime but differ
only in the invoked command name, there is one scientific execution, not two
independent confirmations.

### Follow One Alias Result

```text
environment lock (phylogenetic version)
        -> resolved dependency (bijux-phylogenetics version)
        -> alias import or command
        -> canonical public surface
        -> result manifest (both distribution identities)
```

A support report that records only `phylogenetic` loses the implementation
identity. A report that records only `bijux-phylogenetics` loses the route the
consumer exercised. Preserve both; compare behavior at the canonical surface
and compatibility at the facade boundary.

## Name Mapping

| Canonical | Compatibility name |
| --- | --- |
| distribution `bijux-phylogenetics` | distribution `phylogenetic` |
| `import bijux_phylogenetics` | `import phylogenetic` |
| `from bijux_phylogenetics.api import …` | `from phylogenetic.api import …` |
| command `bijux-phylogenetics` | command `phylogenetic` |

The scientific runtime, external-engine requirements, artifact semantics,
quality boundaries, and Evidence Book are shared.

## Preserve Alias Identity In Artifacts

An artifact created through `phylogenetic` has one scientific owner and two
relevant distribution identities. Retain the invoked compatibility name for
reproduction and the resolved `bijux-phylogenetics` version for computation
ownership.

| Artifact fact | Why it matters |
| --- | --- |
| requested distribution and version | reconstructs how the environment was installed |
| invoked import or command name | reproduces the consumer-facing route |
| resolved canonical distribution and module path | identifies the implementation that executed |
| result and schema version | establishes the shared serialization contract |
| packaged-resource identity | proves the alias consumed canonical resources rather than a local copy |

Do not rewrite canonical owner fields to the alias name, and do not drop the
alias invocation from an operational record. The first would invent a second
scientific owner; the second would hide the route whose compatibility is under
review.

## Hand Off Across Names Without Forking Identity

A downstream consumer may load a result through either import name. Preserve
the original invoked name, both distribution versions, canonical module owner,
result or schema version, and artifact checksums. Do not regenerate a result
solely to replace `phylogenetic` with `bijux_phylogenetics` in presentation.

When moving source code to the canonical import, compare structured results
and manifests under the same environment and inputs. A successful name change
establishes interface compatibility; it is not a new scientific execution and
must not enlarge the evidence denominator.

## Workflow Example

```python
from pathlib import Path

from phylogenetic.api import (
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

sequence_run.write_json(Path("artifacts/sequence-to-tree/workflow.json"))
comparative_run.write_tsv(Path("artifacts/comparative-model.tsv"))
```

## Choosing A Distribution Name

Use `phylogenetic` when concise commands and imports are valuable. Use the
canonical `bijux-phylogenetics` name in publications, environment manifests,
release notes, and support reports when package identity should be explicit.

Installing both does not create independent runtimes. When diagnosing an
environment, record both the requested compatibility distribution and the
resolved canonical runtime version.

## Compatibility Contract

The alias is verified for import behavior, commands, version alignment,
package resources, examples, and built-artifact installation. Report a defect
when a documented canonical workflow cannot be expressed through the alias or
when metadata implies different scientific capability.

| Compatibility layer | Contract |
| --- | --- |
| distribution | installs the intended canonical runtime version |
| import facade | exports the supported canonical public API |
| command facade | preserves parsing, exit status, streams, and structured output |
| packaged resources | resolve from the canonical runtime without copied scientific ownership |
| result contract | preserves model identity, diagnostics, serialization, and artifact meaning |

Compatibility covers naming, imports, commands, resources, and version
alignment. It does not create a second implementation, independent release
line, or separate evidence base. External-engine availability and scientific
limits remain identical under both names.

## Compatibility Failure Classes

| Failure | Meaning | Evidence to retain |
| --- | --- | --- |
| alias distribution missing | packaging or dependency resolution failed | installer output and environment package list |
| `import phylogenetic` fails while canonical import succeeds | facade packaging or import contract failed | traceback and both distribution versions |
| documented alias symbol is absent | public export correspondence failed | symbol path and canonical equivalent |
| alias command differs in parsing or exit state | CLI facade contract failed | commands, streams, exit codes, and structured outputs |
| packaged resource differs | distribution assembly contract failed | wheel hashes and resource inventories |
| numerical result differs after identical dispatch | canonical-runtime identity or serialization path must be investigated | full configurations, manifests, results, and loaded module paths |

Scientific disagreement is not an expected alias feature. The short name may
change how the runtime is requested, never the model, algorithm, or evidence
contract.

## Diagnose A Difference

1. Capture `phylogenetic --version` and `bijux-phylogenetics --version`.
2. Compare the resolved installed distributions in the same environment.
3. Run the smallest documented workflow through both command names.
4. Compare structured outputs and manifests before comparing HTML rendering.
5. Report divergent import, command, resource, metadata, or serialization
   behavior as an alias-contract defect.

Do not repair a difference by copying source into the alias package. The alias
must continue to delegate to the canonical implementation; duplicated runtime
code would create two scientific owners and make future evidence ambiguous.

## Documentation

- [Product handbook](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- [Scientific domains](https://bijux.io/bijux-phylogenetics/02-bijux-phylogenetics-domains/)
- [Distribution map](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/architecture/package-map/)
- [Surface selection](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/surface-selection/)
- [Workflow operations](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/operations/)
- [Result and refusal semantics](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/result-and-refusal-semantics/)
- [Reproducibility records](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/reproducibility-records/)
- [Artifact contracts](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/interfaces/artifact-contracts/)
- [Evidence Book](https://bijux.io/bijux-phylogenetics/03-bijux-phylogenetics-evidence-book/)
- [Canonical package source](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)

Apache License 2.0.
