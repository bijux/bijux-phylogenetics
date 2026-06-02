# Closure Criteria

These criteria define what the repository must prove before it can claim
hard closure for foundational trust, visualization competition, performance competition,
migration usefulness, or reviewer readiness.

- criteria: `5`
- current freshness-tracked bundles: `10`
- integrity-tracked bundles: `19`

## foundational-numerical-trust

- current status: `bounded`

Pass when:
- Major runtime analytical surfaces no longer include uncovered foundational categories.
- Unresolved scalar mismatches are closed or explicitly relegated to bounded, non-foundational surfaces.
- Coverage-gap debt for core comparative and tree workflows has been reduced to bounded, explained exceptions.

Blocking factors:
- 2 analytical surfaces are still uncovered by governed evidence families.
- 12 unresolved scalar mismatch rows remain visible in the mismatch archive.
- 4 coverage-gap debt entries remain reviewer-visible.

Supporting evidence:
- `analytical-surface-coverage.json`
- `mismatch-archive.json`
- `coverage-gaps.json`

## visualization-competition

- current status: `blocked_by_foundational_numerical_trust`

Pass when:
- Foundational numerical trust is already closed for the compared workflow family.
- Figure-comparison evidence bundles exist for the specific visual surface being claimed competitive.

Blocking factors:
- Foundational numerical trust is still bounded rather than closed.
- Plot-only scientific debt remains visible for the PCM1 visual surfaces.

Supporting evidence:
- `scientific-debt-register.json`
- `claim-reaudit.json`

## performance-competition

- current status: `blocked_by_foundational_numerical_trust`

Pass when:
- Foundational numerical trust is already closed for the workflow family.
- Representative benchmark evidence bundles compare Bijux against the reference toolchain under governed inputs.

Blocking factors:
- No governed benchmark evidence families yet exist in the evidence-book.
- Foundational numerical trust is still bounded rather than closed.

Supporting evidence:
- `analytical-surface-coverage.json`
- `completion-gates.json`

## migration-usefulness

- current status: `bounded`

Pass when:
- Migration studies cover more than one major analytical family for each public migration claim.
- Side-by-side R and Bijux mappings exist for the claimed workflow family and link to governed Evidence IDs.

Blocking factors:
- Migration evidence is currently strongest for Lund-derived comparative workflows only.
- Alignment and distance migration claims would currently outrun the evidence-book.

Supporting evidence:
- `teaching-and-migration.json`
- `analytical-surface-coverage.json`

## reviewer-readiness

- current status: `bounded`

Pass when:
- A skeptical reviewer can find linked evidence status for every major analytical surface.
- Coverage gaps, unresolved mismatches, and bounded claims remain explicit in one navigation path.
- High-level docs link directly to the evidence-book closure surfaces.

Blocking factors:
- 3 surfaces still rely on bounded rather than evidence-grounded coverage.
- 0 false-confidence surface findings would block reviewer readiness if nonzero.

Supporting evidence:
- `claim-reaudit.json`
- `analytical-surface-coverage.json`
- `coverage-gaps.json`
