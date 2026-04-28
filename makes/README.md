# Make Architecture

`makes/` is the operational control plane for this repository. It should keep
the root entrypoints small, keep shared implementation synchronized from
`.bijux/shared/bijux-makes-py/` into `makes/bijux-py/`, and keep package-specific policy in
leaf package profiles.

## Design Goals

- expose a small, durable root command surface
- keep shared implementation byte-identical across sibling repositories
- let package profiles declare identity and local policy without copying gate
  logic
- keep repository-specific behavior local when it depends on package catalog,
  paths, release policy, or API freeze policy

## Layout

```text
shared/
└── bijux-py/
```

```text
makes/
├── bijux-py/
│   ├── ci/
│   ├── repository/
│   ├── root/
│   ├── api-contract.mk
│   ├── api-freeze.mk
│   ├── api-live-contract.mk
│   ├── api.mk
│   ├── bijux.mk
│   ├── package-catalog.mk
│   └── package.mk
├── packages/
├── api-freeze.mk
├── env.mk
├── packages.mk
├── publish.mk
└── root.mk
```

`.bijux/shared/bijux-makes-py/` is the system-level shared make layer that must stay
byte-identical across sibling repositories. `makes/bijux-py/` is the local
execution mirror used by root include paths.

Repositories may add another top-level makefile when a repository-specific
archetype earns a durable home, but that file should extend the same layout
instead of copying shared logic out of `.bijux/shared/bijux-makes-py/`.

## Neighbor Contract

Shared gates in `.bijux/shared/bijux-makes-py/` and its local mirror consume the repository
`configs/` tree rather than hardcoded tool flags. Every repository is expected
to expose the same config surface:

```
configs/
├── coveragerc.ini
├── deptry.toml
├── mypy.ini
├── package-lock.json
├── package.json
├── pytest.ini
├── ruff.toml
└── schemathesis.toml
```

## Layer Boundaries

- `.bijux/shared/bijux-makes-py/` is the system-level shared make source of truth.
- `makes/bijux-py/` is the local include mirror that should stay identical to
  `.bijux/shared/bijux-makes-py/`.
- `bijux-py/bijux.mk` verifies both local mirror integrity and cross-repository
  shared integrity.
- `bijux-py/ci/` owns shared gates such as `lint`, `test`, `quality`,
  `security`, `build`, `docs`, and `sbom`.
- `bijux-py/api*.mk` owns shared API contract generation, live contract checks,
  freeze-mode behavior, and API dispatch helpers.
- `bijux-py/package.mk` owns the shared package gate surface and the common
  package bootstrap.
- `bijux-py/package-catalog.mk` owns shared package catalog traversal and root
  dispatch helpers.
- `bijux-py/repository/` owns shared repository-level checks for environment,
  publication, config layout, and make layout.
- `bijux-py/root/` owns shared root help, docs, lifecycle, and package dispatch
  orchestration.
- `packages/*.mk` are leaf package profiles. They should declare package
  identity, include the right archetypes, and define only package-local
  overrides.
- `env.mk`, `packages.mk`, `publish.mk`, and `root.mk` are repository-owned
  entrypoints and policy layers.
- `api-freeze.mk` stays repository-local because freeze and drift policy is tied
  to repository-owned contracts and release intent.

## Placement Rules

1. If the behavior must stay identical in every repository, put it in
   `.bijux/shared/bijux-makes-py/` and mirror it into `makes/bijux-py/`.
2. If the behavior is shared inside one repository but depends on that
   repository's package catalog, paths, or publication policy, keep it in a
   clearly named top-level makefile under `makes/`.
3. If the behavior belongs to exactly one package, keep it in
   `makes/packages/<package>.mk`.
4. If a package profile starts carrying real workflow logic instead of
   declarations, move that logic up into a shared module or a repository-local
   archetype.

## What Does Not Belong Here

- duplicated gate recipes across package profiles
- repository identity or path assumptions embedded in `.bijux/shared/bijux-makes-py/`
- long-lived migration wrappers that only preserve old names
- file names based on temporary planning language instead of durable intent

## Working Rules

- root makefiles are stable entrypoints, not dumping grounds
- package profiles should stay declarative and short
- shared gates should read from `configs/`, not from hand-written per-recipe
  tool flags
- API freeze behavior stays local, while shared API workflow logic stays in
  `.bijux/shared/bijux-makes-py/`
- every new makefile should earn a durable domain name that will still make
  sense years later

## Verification

- `make help`, `make list`, and `make list-all` expose the root command surface
- `make check-shared-bijux-py` verifies local mirror integrity and shared
  byte-identity across repositories
- `make check-config-layout` verifies that the repository exports the expected
  shared config surface
- `make check-make-layout` verifies that the repository keeps the expected make
  directories, entrypoints, and package profiles

## Refactoring Test

Before adding a rule, ask three questions:

1. Does every repository need the exact same implementation?
2. Is this repository policy rather than shared framework behavior?
3. Is this owned by one package only?

Choose the highest honest layer. The best `makes/` tree is the one where each
rule has one obvious home and no file exists just to repeat another one.
