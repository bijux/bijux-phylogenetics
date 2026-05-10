# Analytical Surface Coverage

This map records which major runtime analysis surfaces are already backed
by governed evidence-book bundles and which remain bounded or uncovered.

- analytical surfaces: `6`

## Coverage Status Counts

- `bounded`: `4`
- `uncovered`: `2`

## Surfaces

### Tree diagnostics and tree-data correspondence

- surface id: `tree-diagnostics`
- coverage: `bounded`
- linked evidence bundles: `3`
- summary: The evidence-book proves governed lecture-derived tree import, diagnostics, and tree-data correspondence behavior, but it does not yet close the broader generic tree validation surface.

Runtime surfaces:
- `validate`
- `inspect`
- `topology`
- `dataset crosswalk`

Linked evidence:
- `primate-longevity-signal/evidence-006` — `matched`
- `primate-longevity-signal/evidence-007` — `matched`
- `primate-longevity-signal/evidence-008` — `matched`

Known gaps:
- No dedicated governed evidence family yet covers generic malformed-tree fixtures outside the Lund-derived studies.
- Topology and branch-length comparison workflows remain runtime-tested but not yet mapped into first-class evidence families.

### Taxon identity and normalization

- surface id: `taxon-identity`
- coverage: `bounded`
- linked evidence bundles: `1`
- summary: The taxon-trust study makes taxon identity, synonym handling, and exclusion reasoning reviewer-visible, but the broader taxonomy surface still relies on one evidence family.

Runtime surfaces:
- `taxonomy`
- `normalize-taxa`
- `taxon crosswalk`
- `tree accepted-name export`

Linked evidence:
- `taxon-trust/evidence-001` — `matched`

Known gaps:
- Additional evidence families are still needed for rank inference, accepted-name export, and cross-dataset reconciliation edge cases.

### Alignment diagnostics and trimming

- surface id: `alignment-diagnostics`
- coverage: `uncovered`
- linked evidence bundles: `0`
- summary: The runtime exposes a broad alignment surface, but the evidence-book does not yet contain governed parity or trust bundles for it.

Runtime surfaces:
- `alignment trim`
- `alignment coding`
- `alignment translate`
- `alignment identity-matrix`

Linked evidence:
- none yet

Known gaps:
- No evidence family yet covers trimming, coding-quality, translation, or identity-matrix behavior.

### Comparative models, signal testing, and regression

- surface id: `comparative-analysis`
- coverage: `bounded`
- linked evidence bundles: `15`
- summary: Comparative workflows are the strongest evidence-grounded surface in the repository, but EB, transformed-tree, and ancestral parity remain open coverage boundaries and one diagnostic scalar mismatch is still unresolved.

Runtime surfaces:
- `comparative readiness`
- `comparative signal`
- `comparative pgls`
- `ancestral compare`

Linked evidence:
- `primate-longevity-signal/evidence-001` — `matched_with_tolerance`
- `primate-longevity-signal/evidence-002` — `matched`
- `primate-longevity-signal/evidence-003` — `matched`
- `primate-longevity-signal/evidence-004` — `matched`
- `primate-longevity-signal/evidence-005` — `matched`
- `primate-longevity-signal/evidence-009` — `matched`
- `primate-pgls-and-signal/evidence-001` — `matched`
- `primate-pgls-and-signal/evidence-002` — `matched`
- `primate-pgls-and-signal/evidence-003` — `matched_with_tolerance`
- `primate-pgls-and-signal/evidence-004` — `matched_with_tolerance`
- `primate-pgls-and-signal/evidence-005` — `matched_with_tolerance`
- `primate-pgls-and-signal/evidence-006` — `not_comparable`
- `comparative-trust-boundaries/evidence-001` — `matched`
- `comparative-trust-boundaries/evidence-002` — `matched_with_tolerance`
- `comparative-trust-boundaries/evidence-003` — `matched`

Known gaps:
- Transformed-tree EB parity is still explicit debt rather than closed evidence.
- Ancestral-mode parity from the Lund PCM2 lecture is still intentionally unclaimed.

### Distance analysis and distance maturity

- surface id: `distance-analysis`
- coverage: `uncovered`
- linked evidence bundles: `0`
- summary: Distance workflows have runtime and maturity tests, but they are not yet mapped into governed evidence-book families.

Runtime surfaces:
- `distance matrix workflows`
- `distance maturity reports`

Linked evidence:
- none yet

Known gaps:
- No checked-in evidence bundle yet proves governed distance outputs or cross-tool parity.

### Reviewer-facing reporting and evidence rendering

- surface id: `reporting-surfaces`
- coverage: `bounded`
- linked evidence bundles: `2`
- summary: Reviewer summaries, catalog indexes, and governed study reports are checked in and portable, but the broader runtime report families are not yet all represented by evidence families.

Runtime surfaces:
- `report tree`
- `report dataset`
- `figure packaging`
- `reviewer summaries`

Linked evidence:
- `primate-longevity-signal/evidence-001` — `matched_with_tolerance`
- `taxon-trust/evidence-001` — `matched`

Known gaps:
- Generic HTML reporting surfaces still depend more on runtime tests than on evidence-book-backed coverage.
