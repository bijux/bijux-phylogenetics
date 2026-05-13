# Evidence Review Ritual

This ritual keeps parity gaps, stale bundles, and broad confidence claims
under recurring maintainers review instead of treating the evidence-book as one-off work.

- rituals: `3`
- current maturity tier snapshot: `reviewable_but_incomplete`

## weekly-evidence-triage

- cadence: `weekly`
- goal: Keep freshness, integrity, and open gaps visible before they drift into release-time surprises.

Required inputs:
- `freshness-report.json`
- `integrity-report.json`
- `coverage-gaps.json`
- `mismatch-archive.json`

Commands:
- `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c "from pathlib import Path; from bijux_phylogenetics.evidence.workbench import refresh_evidence_book; refresh_evidence_book(Path('.'))"`
- `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c "from pathlib import Path; from bijux_phylogenetics.evidence.book import validate_evidence_book; report = validate_evidence_book(Path('.')); print({'valid': report.valid, 'issue_count': len(report.issues)})"`

## pre-release-closure-review

- cadence: `before release tags`
- goal: Re-check broad confidence claims against the current closure state before public release messaging goes out.

Required inputs:
- `claim-reaudit.json`
- `closure-criteria.json`
- `evidence-maturity-scorecard.json`
- `completion-gates.json`

Commands:
- `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 pytest packages/bijux-phylogenetics/tests/test_evidence_book_contract.py packages/bijux-phylogenetics/tests/test_evidence_book_repository.py packages/bijux-phylogenetics/tests/test_evidence_book_outputs.py packages/bijux-phylogenetics/tests/test_evidence_workbench.py`
- `UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c "from pathlib import Path; from bijux_phylogenetics.evidence.closure import build_claim_reaudit, build_evidence_maturity_scorecard; print(build_claim_reaudit(Path('.'))['downgraded_claim_count']); print(build_evidence_maturity_scorecard(Path('.'))['maturity_tier'])"`

## new-study-intake-review

- cadence: `for every new evidence study family`
- goal: Force new study intake to declare evidence ownership, closure effects, and unresolved boundaries up front.

Required inputs:
- `analytical-surface-coverage.json`
- `completion-gates.json`
- `scientific-debt-register.json`

Commands:
- `bijux-phylogenetics evidence book studies --json`
- `bijux-phylogenetics evidence book build --json`
- `bijux-phylogenetics evidence book validate --json`
