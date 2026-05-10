# Primate PGLS and signal evidence study Student-Safe Reproducibility

- study id: `primate-pgls-and-signal`
- categories: `teaching-study, migration-study`
- supported scope: regenerate the governed PCM2 study outputs from the checked-in study sources and repository-managed reference script
- entrypoint: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python evidence-book/studies/primate-pgls-and-signal/build_evidence.py`

## Portable Prerequisites

- Python 3.11
- UV project environment at artifacts/root/venv
- Rscript on PATH
- optional R_LIBS_USER rooted at artifacts/root/r-library

## Forbidden Assumptions

- no workstation-local /Users paths
- no sibling repository checkout assumptions
- no hidden manual patching of reference outputs

## Expected Paths

- `evidence-book/studies/primate-pgls-and-signal/evidence-001/scalar-parity-table.json`
- `evidence-book/studies/primate-pgls-and-signal/teaching-guide.json`
- `evidence-book/studies/primate-pgls-and-signal/migration-guide.json`
- `evidence-book/studies/primate-pgls-and-signal/student-safe-reproducibility.json`

## Not Claimed

- EB and ancestral-mode parity remains an explicit coverage boundary
- plot rendering is summarized as diagnostics rather than figure equivalence
