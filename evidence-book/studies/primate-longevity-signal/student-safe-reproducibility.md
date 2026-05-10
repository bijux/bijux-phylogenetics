# Primate Longevity Signal Student-Safe Reproducibility

- study id: `primate-longevity-signal`
- categories: `teaching-study, migration-study`
- supported scope: regenerate the governed teaching and migration summaries from the checked-in evidence bundle and study indexes
- entrypoint: `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c "from pathlib import Path; from bijux_phylogenetics.evidence.book import write_evidence_book_index; write_evidence_book_index(Path('.'))"`

## Portable Prerequisites

- Python 3.11
- UV project environment at artifacts/root/venv

## Forbidden Assumptions

- no workstation-local /Users paths
- no sibling repository checkout assumptions
- no hidden lecture data paths outside the repository

## Expected Paths

- `evidence-book/studies/primate-longevity-signal/teaching-guide.json`
- `evidence-book/studies/primate-longevity-signal/migration-guide.json`
- `evidence-book/studies/primate-longevity-signal/student-safe-reproducibility.json`
- `evidence-book/index/teaching-and-migration.json`

## Not Claimed

- raw Lund workbook reconstruction is not re-executed from inside this repository
- rendered-figure equivalence is still outside the current teaching trust claim
