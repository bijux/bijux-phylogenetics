# Security Policy

This repository uses coordinated vulnerability disclosure. Security reports are
handled privately until we understand impact and have a remediation path.

## What This Policy Covers

Security reports may cover:

- the monorepo layout and shared tooling
- package source code under `packages/`
- shared configuration under `configs/`
- shared automation under `makes/`
- published release artifacts generated from this repository

This policy also covers checked-in API contracts when the issue could mislead
operators, weaken validation, or permit unsafe behavior at an integration
boundary.

## What To Report

Examples of in-scope reports include:

- authentication or authorization bypass
- unsafe defaults that expose data or execution surfaces
- contract validation gaps that allow malformed or dangerous inputs through
- supply-chain weaknesses in build, publish, or artifact handling
- secrets exposure in tracked files or generated release artifacts

Low-quality bug reports and routine support questions should go through normal
issue channels instead.

## Reporting a Vulnerability

Please report security issues privately through one of these channels:

- GitHub private vulnerability reporting for the affected package or repository
- Email: [`bijan@bijux.io`](mailto:bijan@bijux.io)

Please include:

- affected package or shared subsystem
- version or commit reference
- impact summary
- reproduction steps or proof of concept
- any known mitigations

Do not open a public issue for unpatched vulnerabilities.

## What Happens Next

After a report is received, we will:

1. acknowledge receipt
2. confirm whether the report is in scope
3. assess severity and affected surfaces
4. prepare or plan remediation
5. coordinate disclosure when it is responsible to do so

## Response Targets

Best-effort targets:

- acknowledgement within 48 hours
- initial assessment within 5 business days
- fix windows based on severity and exploitability

These are targets, not guarantees, but they are the bar we intend to meet.

## Safe Harbor

Good-faith research is welcome if it:

- avoids privacy violations and service disruption
- stays within systems and data you control
- stops immediately after proving impact
- is reported privately without public disclosure before remediation

Please do not use social engineering, denial-of-service, destructive payloads,
or unauthorized data access as part of research.
