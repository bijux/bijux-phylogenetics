# Real-Local Validation

The `real_local/` suite is the governed installed-environment validation lane.

It is intentionally separate from the fast `engine_contract` tests that use
fixture executables and parser-only contract coverage.

Use this lane when you need evidence that Bijux still works with installed
artifacts and installed external engines, not only with fake runners:

- `test_installability_smoke.py` builds wheel and sdist artifacts, installs
  each into a clean virtual environment, runs `bijux-phylogenetics --help`,
  validates packaged FASTA input, renders a tree report from packaged example
  data, runs a comparative signal command from packaged example data, and
  verifies that the built distributions include governed packaged resources.

- `test_engine_preflight_lane.py` exercises `phylo preflight` against the local
  MAFFT, trimAl, and IQ-TREE environment and verifies that the selected-workflow
  gate honestly blocks `beast-posterior` when no BEAST executable is available.
- `test_alignment_engine_validation_matrix.py` exercises MAFFT, trimAl,
  IQ-TREE, and FastTree on small real inputs and writes one matrix JSON artifact
  under the pytest temporary directory through the product-owned
  `run_alignment_engine_validation_matrix()` surface.
- `test_bayesian_engine_validation_matrix.py` exercises a real MrBayes run and
  records BEAST either from a live executable run when available or from the
  governed checked-in BEAST posterior corpus when the executable is absent.
  That governed corpus now includes the prepared XML, posterior log, posterior
  trees, consensus tree, maximum clade credibility tree, and checked burn-in or
  ESS reference summaries under one durable shared fixture owner, all routed
  through `run_bayesian_engine_validation_matrix()`.
- `test_external_engine_validation_matrix.py` merges those governed alignment
  and Bayesian cases into one reviewer-facing matrix over MAFFT, trimAl,
  IQ-TREE, FastTree, MrBayes, and BEAST through
  `run_external_engine_validation_matrix()`.
- `test_bayesian_execution_lane.py` also proves that `adapter beast-run
  --resume` and `adapter mrbayes-run --resume` only reuse one verified
  completed manifest when the installed executable, checked inputs, and
  recorded outputs still match.
- `test_bayesian_timeout_resume_safety.py` is the fixture-backed contract lane
  for the shared Bayesian posterior safety runtime. It keeps BEAST and
  MrBayes identical on timeout markers, killed-process markers, malformed
  output rejection, missing-executable handling, clean restart, and verified
  resume behavior.

Each matrix case records the reviewer-facing engine name, validation name,
validation mode, executable path when one was used, version text, command,
exit code, runtime, output paths, and output hashes.

Run the full real-local lane with:

```bash
./.venv/bin/pytest tests/real_local -m real_local
```

Run only the external-engine subset with:

```bash
./.venv/bin/pytest tests/real_local -m engine_real
```
